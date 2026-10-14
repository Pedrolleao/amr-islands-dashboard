"""32 - the island browse payload.

One record per kept island - tier 1 (two or more acquired resistance clusters) and tier 2 (one).
The 19,084 tier-3 backbone-only controls are deliberately dropped: they exist to calibrate the
tiering, not to be browsed.

Records are positional arrays against a field list declared once. Written as objects with their key
names repeated, 46,921 islands would spend about a megabyte on the word "scientific_name" alone.

The gene lists are NOT here. They are an order of magnitude bigger and go into genemaps.json
(step 33), which the page decodes lazily for whichever island is selected.

Usage: python3 scripts/32_island_payload.py [--group G]
"""
import importlib.util
import os
import time

from amr_common import ISLANDS, read_tsv

_spec = importlib.util.spec_from_file_location(
    "dash", os.path.join(os.path.dirname(os.path.abspath(__file__)), "30_dashboard_common.py"))
dash = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(dash)

KEEP_TIERS = ["tier1_multi_acquired", "tier2_single_acquired", "tier3_core_only"]
# tier 3 (backbone-only windows around a core anchor) is carried since 24 September; the page shows
# tiers 1 and 2 by default and tier 3 on request

FIELDS = ["tier", "target", "pds", "species", "contig", "ord_lo", "n_genes", "bp_lo",
          "span_kb_x10", "window_complete", "frac_window_x1000", "truncation", "n_anchor",
          "n_anchor_cargo", "n_anchor_core", "n_anchor_point", "anchor_profile",
          "amr_gene_families", "n_amr_families", "amr_classes", "n_stress_genes",
          "n_virulence_genes", "n_mge_genes", "frac_mge_x1000", "mge_kinds", "integron",
          "n_pseudo", "frac_hypothetical_x1000", "n_clusters", "type"]

TYPE_FIELDS = ["recurrence_tier", "n_islands", "n_genomes", "n_pds", "clonality_x100", "n_species",
               "n_amr_families", "amr_gene_families", "median_n_genes", "frac_truncated_x1000",
               "frac_with_mge_x1000", "frac_integron_x1000", "rep_island"]


def num(v, scale=1, default=0):
    try:
        f = float(v)
    except (TypeError, ValueError):
        return default
    return default if f != f else int(round(f * scale))


def main():
    t0 = time.time()
    isl = read_tsv(ISLANDS / "islands.tsv")
    total = len(isl)
    isl = isl[isl.island_tier.isin(KEEP_TIERS)].copy()
    print(f"islands: {len(isl):,} kept of {total:,} ({' + '.join(KEEP_TIERS)})", flush=True)

    member = read_tsv(ISLANDS / "island_type_membership.tsv",
                      usecols=["island_id", "island_type_id"])
    type_of = dict(zip(member.island_id, member.island_type_id))

    d = dash.Dicts()
    ids, recs = [], []
    for r in isl.itertuples():
        tid = type_of.get(r.island_id, "")
        ids.append(dash.island_num(r.island_id))
        recs.append([
            KEEP_TIERS.index(r.island_tier),
            d.enc("target", r.target_acc), d.enc("pds", r.PDS_acc),
            d.enc("species", r.scientific_name), d.enc("contig", r.contig),
            num(r.ord_lo), num(r.n_genes), num(r.bp_lo), num(r.span_kb, 10),
            num(r.window_complete), num(r.frac_window, 1000),
            d.enc("truncation", r.truncation),
            num(r.n_anchor), num(r.n_anchor_cargo), num(r.n_anchor_core), num(r.n_anchor_point),
            d.enc("profile", r.anchor_profile),
            d.enc("families", r.amr_gene_families), num(r.n_amr_families),
            d.enc("amrclass", r.amr_classes),
            num(r.n_stress_genes), num(r.n_virulence_genes), num(r.n_mge_genes),
            num(r.frac_mge, 1000), d.enc("mge", r.mge_kinds), num(r.has_integron_integrase),
            num(r.n_pseudo), num(r.frac_hypothetical, 1000), num(r.n_clusters),
            dash.type_num(tid) if tid else -1,
        ])

    # only the types that a kept island actually belongs to
    referenced = {rec[FIELDS.index("type")] for rec in recs} - {-1}
    types_all = read_tsv(ISLANDS / "island_types.tsv")
    types_all["num"] = types_all.island_type_id.map(dash.type_num)
    types_all = types_all[types_all.num.isin(referenced)]
    types = {}
    for r in types_all.itertuples():
        types[int(r.num)] = [
            d.enc("rtier", r.recurrence_tier), num(r.n_islands), num(r.n_genomes), num(r.n_pds),
            num(r.clonality_ratio, 100), num(r.n_species), num(r.n_amr_families),
            d.enc("families", r.amr_gene_families), num(r.median_n_genes),
            num(r.frac_truncated, 1000), num(r.frac_with_mge, 1000), num(r.frac_integron, 1000),
            dash.island_num(r.rep_island_id),
        ]

    cooc = read_tsv(ISLANDS / "amr_family_cooccurrence.tsv")
    cooc_out = [[dash.cluster_num(r.cluster_a), dash.cluster_num(r.cluster_b),
                 d.enc("gene", r.family_a), d.enc("gene", r.family_b), num(r.n_pds)]
                for r in cooc.itertuples()]

    payload = {"built": time.strftime("%Y-%m-%d"), "fields": FIELDS, "type_fields": TYPE_FIELDS,
               "dicts": d.as_dict(), "ids": ids, "isl": recs, "types": types, "cooc": cooc_out}
    mb = dash.write_payload(payload, "islands.json")

    print("\nchecks")
    print(f"  islands: {len(recs):,}; tier1 {sum(1 for r in recs if r[0] == 0):,}, "
          f"tier2 {sum(1 for r in recs if r[0] == 1):,}, tier3 {sum(1 for r in recs if r[0] == 2):,}")
    print(f"  island types referenced: {len(types):,} of {len(types_all):,} read "
          f"({len(referenced):,} distinct ids in the kept set)")
    print(f"  co-occurrence edges: {len(cooc_out):,}")
    print(f"  dictionaries: " + ", ".join(f"{k} {d.size(k):,}" for k in sorted(d.as_dict())))
    print(f"  ids are unique: {len(set(ids)) == len(ids)}")
    print(f"  size {mb:.2f} MB")
    print(f"\ndone in {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
