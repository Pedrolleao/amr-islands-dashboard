"""31 - the two cluster payloads for the dashboard.

`clusters.json` is the deep layer: the 806 clusters that have been analysed in full (411 resistance,
328 stress, 67 virulence), each with its variation track, domain architecture, fold matches and the
islands it sits in. This is what the cluster detail view reads.

`cluster_index.json` is the shallow layer: all 208,757 clusters with product, prevalence, pangenome
class and Pfam architecture, so a gene seen in an island always has a name and can be searched for.
It is a structure of arrays, one slot per cluster in id order, with no id stored at all - the ids are
a dense range, so position IS the id.

Sizes: the deep layer is small because there are only 806 of them; the shallow layer measured 3.9 MB,
which is small enough that no prefix index or island-adjacency restriction is needed.

Usage: python3 scripts/31_cluster_payload.py [--group G]
"""
import importlib.util
import os
import time

import pandas as pd

from amr_common import (DASHBOARD, DOMAINS, ISLANDS, PANGENOME, STRUCTURES, VARIATION, read_tsv)

_spec = importlib.util.spec_from_file_location(
    "dash", os.path.join(os.path.dirname(os.path.abspath(__file__)), "30_dashboard_common.py"))
dash = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(dash)

SETS = [("amr", ""), ("stress", "_stress"), ("virulence", "_virulence")]
MAX_ISLANDS_LISTED = 60          # the detail view shows this many, with the true total beside it
MAX_FOLD_HITS = 5
KEEP_TIERS = ("tier1_multi_acquired", "tier2_single_acquired")
TIER3 = "tier3_core_only"          # backbone-only windows, counted separately

FIELDS = ["set", "gene", "product", "n_members", "n_genomes", "prevalence_permille",
          "pangenome_class", "anchor_tier", "amr_classes", "architecture", "n_domains",
          "covered_permille", "aa_len", "is_fragment", "mean_plddt_x10", "min_plddt_x10",
          "ptm_x1000", "max_variation_x10", "source", "status", "catalog_symbol",
          "catalog_identity_x10", "n_lineages", "n_copies", "n_distinct_seqs",
          "mean_identity_x10000", "n_positions_ge_5pct", "has_variation", "n_islands",
          "amr_subclasses", "n_islands_t3"]   # appended last so earlier positions never move


def suffixed(base, suffix):
    return base if suffix == "" else base.with_name(base.name + suffix)


def num(v, scale=1, default=-1):
    try:
        f = float(v)
    except (TypeError, ValueError):
        return default
    if f != f:
        return default
    return int(round(f * scale))


def load_sets():
    """One frame per analysed cluster, carrying which set it came from."""
    frames = []
    for name, suffix in SETS:
        s = read_tsv(suffixed(STRUCTURES, suffix) / "structure_summary.tsv")
        s["set"] = name
        # the amr summary predates the status column; everything in it was predicted
        if "status" not in s.columns:
            s["status"] = "predicted"
        for col in ("source", "parent_cluster_id", "module", "module_start", "module_end"):
            if col not in s.columns:
                s[col] = ""
        frames.append(s)
    s = pd.concat(frames, ignore_index=True)
    s["module"] = s.cluster_id.map(dash.cluster_module)
    return s[s.module == 0].copy()          # module rows belong to structures.json, not here


def cluster_to_islands():
    """Which islands each cluster appears in: (tier-1/2 count, tier-3 count, ids with tier 1/2 first)."""
    isl = read_tsv(ISLANDS / "islands.tsv", usecols=["island_id", "island_tier"])
    tier_of = dict(zip(isl.island_id, isl.island_tier))
    keep = {i for i, t in tier_of.items() if t in KEEP_TIERS or t == TIER3}
    g = read_tsv(ISLANDS / "gene_island.tsv", usecols=["island_id", "cluster_id"])
    g = g[g.island_id.isin(keep) & g.cluster_id.notna() & (g.cluster_id != "")]
    out = {}
    for cid, grp in g.groupby("cluster_id", sort=False):
        ids = set(grp.island_id)
        i12 = sorted(dash.island_num(i) for i in ids if tier_of[i] in KEEP_TIERS)
        i3 = sorted(dash.island_num(i) for i in ids if tier_of[i] == TIER3)
        out[dash.cluster_num(cid)] = (len(i12), len(i3), (i12 + i3)[:MAX_ISLANDS_LISTED])
    return out


def build_deep(ann, structures, islands_of, d):
    """The 806 analysed clusters, with every layer attached."""
    var_summary, var_pos, domains, pfam_rows, folds = {}, {}, {}, {}, {}
    for name, suffix in SETS:
        vf = suffixed(VARIATION, suffix) / "cluster_variation_summary.tsv"
        if vf.exists():
            for r in read_tsv(vf).itertuples():
                var_summary[r.cluster_id] = r
        pf = suffixed(VARIATION, suffix) / "position_variation.tsv"
        if pf.exists():
            p = pd.read_csv(pf, sep="\t", usecols=["cluster_id", "position", "ref_residue",
                                                   "var_frac_lineages", "var_frac_copies",
                                                   "catalog_site"])
            for cid, grp in p.groupby("cluster_id", sort=False):
                var_pos[cid] = grp.sort_values("position")
        for r in read_tsv(suffixed(DOMAINS, suffix) / "cluster_domains.tsv").itertuples():
            domains[r.cluster_id] = r
        pd_file = suffixed(DOMAINS, suffix) / "pfam_domains.tsv"
        for cid, grp in read_tsv(pd_file).groupby("cluster_id", sort=False):
            pfam_rows[cid] = grp.sort_values("start")
        hf = suffixed(DOMAINS, suffix) / "foldseek_hits.tsv"
        if hf.exists():
            for cid, grp in read_tsv(hf).groupby("cluster_id", sort=False):
                folds[cid] = grp.head(MAX_FOLD_HITS)

    pfam_table, pfam_index = [], {}
    recs, var_out, dom_out, fold_out, isl_out = {}, {}, {}, {}, {}
    ann_by_id = ann.set_index("cluster_id")

    for s in structures.itertuples():
        cid = s.cluster_id
        n = dash.cluster_num(cid)
        a = ann_by_id.loc[cid] if cid in ann_by_id.index else None
        dm = domains.get(cid)
        vs = var_summary.get(cid)
        n_isl, n_isl3, ids = islands_of.get(n, (0, 0, []))

        recs[n] = [
            {"amr": 0, "stress": 1, "virulence": 2}[s.set],
            d.enc("gene", s.family), d.enc("product", s.product_top),
            num(s.n_members), num(a.n_genomes if a is not None else -1),
            num(a.prevalence if a is not None else 0, 1000),
            d.enc("class", a.pangenome_class if a is not None else ""),
            d.enc("tier", s.anchor_tier), d.enc("amrclass", s.amr_classes),
            d.enc("arch", dm.architecture if dm is not None else ""),
            num(dm.n_domains if dm is not None else 0, default=0),
            num(dm.fraction_covered if dm is not None else 0, 1000, default=0),
            num(s.aa_len), 1 if str(s.is_fragment) == "True" else 0,
            num(s.mean_plddt, 10), num(s.min_plddt, 10), num(s.ptm, 1000),
            num(s.max_variation_pct, 10), d.enc("source", s.source or "esmfold"),
            d.enc("status", s.status), d.enc("catalog", s.catalog_symbol),
            num(vs.catalog_identity, 10) if vs is not None else -1,
            num(vs.n_lineages) if vs is not None else -1,
            num(vs.n_copies) if vs is not None else -1,
            num(s.n_distinct_seqs),
            num(vs.mean_identity_to_reference_lineages, 10000) if vs is not None else -1,
            num(vs.n_positions_ge_5pct) if vs is not None else -1,
            1 if cid in var_pos else 0, n_isl,
            d.enc("amrsubclass", a.amr_subclasses if a is not None else ""),
            n_isl3,
        ]

        if cid in var_pos:
            g = var_pos[cid]
            sites = [int(p) for p, c in zip(g.position, g.catalog_site.fillna(""))
                     if str(c) not in ("", "nan")]
            top = g.nlargest(5, "var_frac_lineages")[["position", "var_frac_lineages"]]
            var_out[n] = {
                "ref": "".join(str(x) for x in g.ref_residue),
                "lin": dash.qual_string(g.var_frac_lineages * 100),
                "cop": dash.qual_string(g.var_frac_copies * 100),
                "sites": sites,
                "top": [[int(p), num(v, 1000)] for p, v in top.itertuples(index=False)],
            }
        if cid in pfam_rows:
            out = []
            for r in pfam_rows[cid].itertuples():
                key = str(r.pfam_acc)
                if key not in pfam_index:
                    pfam_index[key] = len(pfam_table)
                    pfam_table.append([key, str(r.pfam_name), str(r.pfam_desc)])
                out.append([pfam_index[key], int(r.start), int(r.end)])
            dom_out[n] = out
        if cid in folds:
            fold_out[n] = [[d.enc("db", r.database), d.enc("target", r.target),
                            d.enc("tdesc", r.theader), num(r.alntmscore, 1000),
                            num(r.fident, 1000), 1 if str(r.same_fold) == "True" else 0]
                           for r in folds[cid].itertuples()]
        if ids:
            isl_out[n] = ids

    return recs, var_out, dom_out, fold_out, isl_out, pfam_table


def build_index(ann, arch_all, analysed):
    """All 208,757 clusters as parallel arrays. Position is the id, so no id is stored."""
    ann = ann.sort_values("cluster_id").reset_index(drop=True)
    nums = ann.cluster_id.map(dash.cluster_num).to_numpy()
    if not (nums[0] == 1 and nums[-1] == len(nums) and (nums[1:] - nums[:-1] == 1).all()):
        raise SystemExit("cluster ids are not a dense 1..N range; the index needs an id column")

    arch = ann.cluster_id.map(arch_all).fillna("")
    d = dash.Dicts()
    prod_idx = [d.enc("product", p) for p in ann.product_top.fillna("")]
    arch_idx = [d.enc("arch", a) for a in arch]
    cls_idx = [d.enc("class", c) for c in ann.pangenome_class.fillna("")]
    setcode = {"AMR": 0, "STRESS": 1, "VIRULENCE": 2}
    st = [setcode.get(e, -1) for e in ann.element_type.fillna("")]
    return {
        "n": len(ann),
        "dicts": d.as_dict(),
        "nm": [int(v) for v in ann.n_members.fillna(0)],
        "pv": [int(round(float(v) * 1000)) for v in ann.prevalence.fillna(0)],
        "pi": prod_idx, "ai": arch_idx, "ci": cls_idx, "st": st,
    }


def main():
    t0 = time.time()
    ann = read_tsv(PANGENOME / "cluster_annotation.tsv",
                   usecols=["cluster_id", "element_type", "pangenome_class", "product_top",
                            "prevalence", "n_members", "n_genomes", "gene_symbol_top",
                            "amr_subclasses"])
    structures = load_sets()
    print(f"{len(structures)} analysed clusters; {len(ann):,} in the pangenome", flush=True)

    islands_of = cluster_to_islands()
    print(f"cluster -> island index: {len(islands_of):,} clusters appear in a kept island",
          flush=True)

    d = dash.Dicts()
    recs, var_out, dom_out, fold_out, isl_out, pfam_table = build_deep(ann, structures, islands_of, d)
    deep = {"built": time.strftime("%Y-%m-%d"), "fields": FIELDS, "dicts": d.as_dict(),
            "pfam": pfam_table, "recs": recs, "var": var_out, "dom": dom_out,
            "fold": fold_out, "isl": isl_out}
    mb1 = dash.write_payload(deep, "clusters.json")

    arch_all = read_tsv(DOMAINS.with_name("Domains_all") / "cluster_domains.tsv",
                        usecols=["cluster_id", "architecture"])
    arch_all = dict(zip(arch_all.cluster_id, arch_all.architecture.fillna("")))
    idx = build_index(ann, arch_all, set(structures.cluster_id))
    mb2 = dash.write_payload(idx, "cluster_index.json")

    print("\nchecks")
    print(f"  deep records: {len(recs)} (expected 806); with variation: {len(var_out)}; "
          f"with domains: {len(dom_out)}; with fold hits: {len(fold_out)}")
    print(f"  distinct Pfam families referenced: {len(pfam_table)}")
    print(f"  index: {idx['n']:,} clusters, {len(idx['dicts']['product']):,} distinct products, "
          f"{len(idx['dicts']['arch']):,} distinct architectures")
    print(f"  analysed clusters flagged in the index: {sum(1 for v in idx['st'] if v >= 0)}")
    print(f"  sizes: clusters.json {mb1:.2f} MB, cluster_index.json {mb2:.2f} MB")
    print(f"\ndone in {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
