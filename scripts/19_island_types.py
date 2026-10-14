"""19 - group islands into recurrent TYPES, so the same cassette can be recognised across genomes.

An island's architecture is written as an ordered, strand-aware string of protein-cluster codes.
Two islands with the same string are the same arrangement of the same gene families, wherever they
were found.

ORIENTATION IS NORMALISED. The same cassette inserts in either direction, so the string is stored as
min(forward, reverse) with every strand flipped in the reverse form. Without this, every inverted
copy would look like a second, unrelated type.

THREE KEYS, BECAUSE EXACT MATCHING ALONE IS TOO BRITTLE HERE. 74% of islands are truncated by a
contig edge, so demanding an identical gene list would split one real cassette into dozens of
fragments-of-a-type.

    arch_exact_key    every gene, strand-aware                 strictest
    arch_core_key     only INFORMATIVE genes - AMR, stress,    the default type definition
                      virulence, mobile elements, and rare
                      clusters (prevalence < 0.15). Backbone
                      and unclustered genes are dropped.
    amr_content_key   the set of AMR gene families, no order   "what does it carry"

RECURRENCE IS COUNTED IN SNP CLUSTERS, NOT GENOMES. This dataset is full of clonal outbreak
expansions: 400 genomes of one outbreak carrying one plasmid is ONE acquisition event, and counting
genomes would make it the most "recurrent" type in the study. `n_pds` counts distinct Pathogen
Detection SNP clusters and `clonality_ratio` = genomes / SNP clusters makes the inflation visible.

No p-values. Every anchor is inside an island by construction, so any naive enrichment test against
all islands is circular. The counts are reported plainly instead.

Usage: python3 scripts/19_island_types.py [--group G]
"""
import hashlib
import time

import pandas as pd

from amr_common import GROUP, GROUPS, ISLANDS, read_tsv, write_tsv

RARE_PREV = 0.15          # below this a cluster counts as informative for the core key
MIN_PDS = 3               # distinct SNP clusters before a type is "recurrent across lineages"


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def canonical(tokens):
    """min(forward, reverse-with-strands-flipped) -> orientation-independent architecture."""
    fwd = ",".join(tokens)
    rev = ",".join(("-" + t[1:] if t[0] == "+" else "+" + t[1:]) if t[0] in "+-" else t
                   for t in reversed(tokens))
    return (fwd, "fwd") if fwd <= rev else (rev, "rev")


def main():
    t0 = time.time()
    prefix = GROUPS[GROUP].get("island_prefix", "ISL") + "T"
    gt = read_tsv(ISLANDS / "island_gene_table.tsv")
    isl = read_tsv(ISLANDS / "islands.tsv")
    gt["pos"] = pd.to_numeric(gt["pos_in_island"])
    gt["prev"] = pd.to_numeric(gt["pc_prevalence_eff"], errors="coerce").fillna(0)
    log(f"{len(isl):,} islands, {len(gt):,} genes")

    informative = ((gt["draw_role"] != "other") & (gt["draw_role"] != "hypothetical")) | \
                  ((gt["cluster_id"] != "") & (gt["prev"] < RARE_PREV))
    gt["_tok"] = [("+" if s == "+" else "-") + (c if c else "X")
                  for s, c in zip(gt["strand"], gt["cluster_id"])]

    rows = []
    for iid, grp in gt.sort_values(["island_id", "pos"]).groupby("island_id", sort=False):
        toks = list(grp["_tok"])
        core_toks = list(grp.loc[informative.loc[grp.index], "_tok"])
        exact, orient = canonical(toks)
        core, _ = canonical(core_toks) if core_toks else ("", "fwd")
        rows.append((iid, exact, core, orient, len(core_toks)))
    arch = pd.DataFrame(rows, columns=["island_id", "arch_exact_key", "arch_core_key",
                                       "arch_orientation", "n_informative"])
    arch["arch_hash"] = [hashlib.md5(k.encode()).hexdigest()[:12] for k in arch["arch_core_key"]]
    m = isl.merge(arch, on="island_id")
    m["amr_content_key"] = [";".join(sorted(set(filter(None, v.split(";")))))
                            for v in m["amr_gene_families"]]

    # ---- types, grouped on the core architecture ------------------------------------------------
    types = []
    for key, grp in m[m["arch_core_key"] != ""].groupby("arch_core_key", sort=False):
        pds = set(filter(None, grp["PDS_acc"]))
        n_g = grp["download_acc"].nunique()
        fams = sorted({f for v in grp["amr_gene_families"] for f in v.split(";") if f})
        types.append({
            "arch_core_key": key, "n_islands": len(grp), "n_genomes": n_g, "n_pds": len(pds),
            "n_species": grp["scientific_name"].str.split().str[:2].str.join(" ").nunique(),
            "key_len": grp["n_informative"].iloc[0],
            "amr_gene_families": ";".join(fams), "n_amr_families": len(fams),
            "median_cargo_clusters": pd.to_numeric(grp["n_cargo_clusters"]).median(),
            "amr_classes": ";".join(sorted({c for v in grp["amr_classes"]
                                            for c in v.split(";") if c})),
            "median_n_genes": pd.to_numeric(grp["n_genes"]).median(),
            "frac_truncated": round((grp["truncation"] != "none").mean(), 3),
            "frac_with_mge": round((pd.to_numeric(grp["n_mge_genes"]) > 0).mean(), 3),
            "frac_integron": round(pd.to_numeric(grp["has_integron_integrase"]).mean(), 3),
            "tier_composition": ";".join(f"{k}:{v}" for k, v in
                                         grp["island_tier"].value_counts().items()),
            "clonality_ratio": round(n_g / len(pds), 2) if pds else "",
            "rep_island_id": grp.sort_values("window_complete", ascending=False)["island_id"].iloc[0],
        })
    tdf = pd.DataFrame(types)
    # total sort order, so ids never permute between runs
    tdf = tdf.sort_values(["n_pds", "n_islands", "arch_core_key"],
                          ascending=[False, False, True]).reset_index(drop=True)
    tdf["island_type_id"] = [f"{prefix}{i:06d}" for i in range(1, len(tdf) + 1)]
    tdf["recurrence_tier"] = [
        "recurrent_cross_lineage" if p >= MIN_PDS else
        "recurrent_within_lineage" if g >= 3 else "rare" if g == 2 else "singleton"
        for p, g in zip(tdf["n_pds"], tdf["n_genomes"])]
    cols = ["island_type_id", "recurrence_tier", "n_islands", "n_genomes", "n_pds",
            "clonality_ratio", "n_species", "key_len", "n_amr_families", "amr_gene_families",
            "amr_classes", "median_n_genes", "frac_truncated", "frac_with_mge", "frac_integron",
            "tier_composition", "rep_island_id", "arch_core_key"]
    write_tsv(tdf[cols], "island_types.tsv", folder=ISLANDS)

    member = m.merge(tdf[["arch_core_key", "island_type_id"]], on="arch_core_key", how="left")
    write_tsv(member[["island_id", "island_type_id", "arch_hash", "arch_orientation",
                      "arch_core_key", "arch_exact_key", "amr_content_key", "n_informative"]],
              "island_type_membership.tsv", folder=ISLANDS)

    # ---- what travels with what ------------------------------------------------------------------
    # Pairs of distinct CARGO CLUSTERS, not of family labels: one cluster carries several labels for
    # the same gene (aac(6')-Ib / -Ib' / -Ib-cr / -Ib11), which would otherwise fill the table with
    # a gene "co-occurring" with itself.
    anchors = read_tsv(ISLANDS / "anchor_clusters.tsv")
    label_of = {}
    for r in anchors.itertuples():
        fams = sorted(set(filter(None, r.amr_gene_families.split(";"))))
        label_of[r.cluster_id] = (fams[0] if fams else r.cluster_id)
    pairs = {}
    for r in m.itertuples():
        cl = sorted(set(filter(None, r.cargo_cluster_ids.split(";"))))
        pds = r.PDS_acc
        for i in range(len(cl)):
            for j in range(i + 1, len(cl)):
                pairs.setdefault((cl[i], cl[j]), set()).add(pds or r.download_acc)
    co = pd.DataFrame([{"cluster_a": a, "cluster_b": b, "family_a": label_of.get(a, a),
                        "family_b": label_of.get(b, b), "n_pds": len(v)}
                       for (a, b), v in pairs.items()])
    if len(co):
        co = co.sort_values("n_pds", ascending=False)
    write_tsv(co, "amr_family_cooccurrence.tsv", folder=ISLANDS)

    # ---- QC -----------------------------------------------------------------------------------------
    print("\n===== QC 19_island_types =====")
    print(f"islands with an informative architecture: {len(m[m['arch_core_key']!='']):,} "
          f"of {len(m):,}")
    print(f"island types: {len(tdf):,}")
    print("\nrecurrence:")
    for t, n in tdf["recurrence_tier"].value_counts().items():
        print(f"  {t:26s} {n:7,} types, {int(tdf.loc[tdf.recurrence_tier==t,'n_islands'].sum()):8,} islands")
    print(f"\ntypes seen in >= {MIN_PDS} SNP clusters: {(tdf['n_pds']>=MIN_PDS).sum():,}")
    print("\ntop 15 types by distinct SNP clusters (n_genomes shown beside, to expose clonality):")
    for r in tdf.head(15).itertuples():
        print(f"  {r.island_type_id} pds={r.n_pds:4,} genomes={r.n_genomes:5,} "
              f"clonality={r.clonality_ratio:5} isl={r.n_islands:5,} "
              f"mge={r.frac_with_mge:.2f} {r.amr_gene_families[:52]}")
    if len(co):
        print("\ntop 12 AMR family pairs travelling together (distinct SNP clusters):")
        for r in co.head(12).itertuples():
            print(f"  {r.n_pds:5,}  {r.family_a} + {r.family_b}")
    log(f"done in {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
