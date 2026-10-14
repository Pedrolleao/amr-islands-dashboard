"""22 - the hard checks, and the generated README for the Islands/ folder.

The checks are the ones that would catch a silent, plausible-looking failure. Each prints PASS/FAIL
and the number it was judged on, so a re-run after a new release can be compared line by line.

Usage: python3 scripts/22_island_qc.py [--group G]
"""
import time

import pandas as pd

from amr_common import GROUP, ISLANDS, PANGENOME, ROOT, read_tsv

CONTROL_ISOLATE = "PDT000050772.2"
CONTROL_CONTIG = "NZ_JYFY02000004.1"
results = []


def check(name, ok, detail):
    results.append((name, bool(ok), detail))
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}: {detail}")


def main():
    t0 = time.time()
    isl = read_tsv(ISLANDS / "islands.tsv")
    gt = read_tsv(ISLANDS / "island_gene_table.tsv")
    gi = read_tsv(ISLANDS / "gene_island.tsv")
    ann = read_tsv(PANGENOME / "cluster_annotation.tsv")
    anchors = read_tsv(ISLANDS / "anchor_clusters.tsv")
    link = read_tsv(ISLANDS / "hit_gene_link.tsv")
    genome_index = read_tsv(PANGENOME / "genome_index.tsv")
    types = read_tsv(ISLANDS / "island_types.tsv")
    for c in ("n_genes", "ord_lo", "ord_hi", "window_complete", "n_cargo_clusters", "n_mge_genes"):
        isl[c] = pd.to_numeric(isl[c])
    n_genomes = len(genome_index)

    print("===== QC 22_island_qc =====")

    # 1. geometry
    check("ordinal contiguity",
          ((isl["ord_hi"] - isl["ord_lo"] + 1) != isl["n_genes"]).sum() == 0,
          f"{((isl['ord_hi'] - isl['ord_lo'] + 1) != isl['n_genes']).sum()} violations of "
          f"ord_hi - ord_lo + 1 == n_genes")
    check("gene table matches island sizes", isl["n_genes"].sum() == len(gt),
          f"sum(n_genes)={isl['n_genes'].sum():,} vs {len(gt):,} rows")
    dup = gi.duplicated(["gene_id"]).sum()
    check("no gene in two islands", dup == 0, f"{dup} genes appear in more than one island")

    # 2. the propagation test - the one number that proves anchoring is at CLUSTER level
    gyr = anchors[anchors["anchor_tier"] == "POINT_ONLY"]
    n_anchor_gyr = pd.to_numeric(gyr["n_members"]).sum()
    hits = link[link["symbol"].str.startswith("gyrA_")]
    n_hit_genomes = hits["target_acc"].nunique()
    check("cluster-level propagation (gyrA)",
          n_anchor_gyr > 1.5 * n_hit_genomes,
          f"{n_anchor_gyr:,} anchors from {n_hit_genomes:,} genomes carrying the mutation "
          f"(hit-level anchoring would give ~{n_hit_genomes:,}; ~{n_genomes:,} is correct)")

    # 3. controls, verified by hand before the pipeline existed
    ctl = isl[(isl["target_acc"] == CONTROL_ISOLATE) & (isl["contig"] == CONTROL_CONTIG)]
    cat = ctl[ctl["direct_hit_symbols"].str.contains("catA2", na=False)]
    if len(cat):
        iid = cat.iloc[0]["island_id"]
        genes = gt[gt["island_id"] == iid]
        tn3 = genes[(genes["product"].str.contains("Tn3 family transposase", na=False))
                    & (genes["is_pseudo"] == "1")]
        check("control 1: pseudogene Tn3 transposase is in the island", len(tn3) > 0,
              f"{iid}: {len(tn3)} pseudogene Tn3 transposase(s) - this gene has NO protein "
              f"accession, so it only exists if pseudogenes were kept in the gene index")
        check("control 1: blaSHV-12 in the same island",
              genes["hit_symbol"].eq("blaSHV-12").any(),
              f"{iid} carries {cat.iloc[0]['direct_hit_symbols']}")
        check("control 1: called tier1", cat.iloc[0]["island_tier"] == "tier1_multi_acquired",
              cat.iloc[0]["island_tier"])
    else:
        check("control 1 present", False, "catA2 island not found")

    cas = ctl[pd.to_numeric(ctl["n_cargo_clusters"]) >= 5]
    check("control 2: the cassette array merged into one island", len(cas) >= 1,
          f"{len(cas)} island(s) with >=5 distinct cargo clusters on {CONTROL_CONTIG}; "
          f"symbols: {cas.iloc[0]['direct_hit_symbols'][:80] if len(cas) else '-'}")
    if len(cas):
        check("control 2: the DNA-level hit ere(A) was placed",
              "ere(A)" in cas.iloc[0]["direct_hit_symbols"],
              "ere(A) has no protein accession, so this tests the coordinate path")

    # 4. hit placement
    tiers = link["match_tier"].value_counts()
    t1 = tiers.get("T1_exact", 0) / max(len(link), 1)
    check("hit placement dominated by exact coordinates", t1 > 0.85,
          f"T1_exact {100*t1:.1f}%; T5_none {tiers.get('T5_none', 0)}")
    strong = link[link["match_tier"] == "T1_exact"]
    agree = pd.to_numeric(strong["protein_agrees"], errors="coerce").dropna()
    check("protein accession agreement at T1", agree.mean() > 0.999,
          f"{100*agree.mean():.3f}% of {len(agree):,}")

    # 5. tiering sanity
    spot = {"blaKPC-2": ("cloud", {"MOBILE_AMR", "SHELL_AMR"}),
            "blaNDM-1": (None, {"MOBILE_AMR", "SHELL_AMR"}),
            "blaCMY": (None, {"CORE_AMR", "SOFTCORE_AMR"})}
    for sym, (pan, want) in spot.items():
        rows = anchors[anchors["amr_symbols"].str.contains(sym, regex=False, na=False)]
        got = set(rows["anchor_tier"])
        check(f"anchor tier of {sym}", bool(got & want), f"{sorted(got)} (expected one of {sorted(want)})")

    # 6. truncation is labelled and plausible
    trunc = (isl["truncation"] != "none").mean()
    check("truncation recorded and in a plausible range", 0.2 < trunc < 0.95,
          f"{100*trunc:.1f}% of islands truncated; "
          f"{100*isl['window_complete'].mean():.1f}% have a full window")

    # 7. the clean subset actually exists
    clean = isl[(isl["island_tier"] == "tier1_multi_acquired") & (isl["window_complete"] == 1)]
    check("a clean island set exists", len(clean) > 200,
          f"{len(clean):,} tier1 islands with a complete window, "
          f"{clean['PDS_acc'].replace('', pd.NA).nunique()} SNP clusters")

    # 8. types are deterministic in shape
    check("island types built", len(types) > 0 and types["island_type_id"].is_unique,
          f"{len(types):,} types, ids unique")

    n_fail = sum(1 for _, ok, _ in results if not ok)
    print(f"\n{len(results) - n_fail}/{len(results)} checks passed")

    # ---- README ------------------------------------------------------------------------------------
    amr_cl = ann[pd.to_numeric(ann["n_amr_hits"]) > 0]
    L = []
    a = L.append
    a(f"# AMR resistance islands - {GROUP.replace('_', ' ')}\n")
    a("Built by `scripts/run_islands.py --group " + GROUP + "` (steps 11-22). Every number here is "
      "computed from the files in this folder; re-run `scripts/22_island_qc.py` after a rebuild.\n")
    a("## What an island is\n")
    a("An **island** is the union of +/-10 gene windows around anchors on one contig, overlapping "
      "windows merged, never spanning a contig. An **anchor** is any gene whose protein cluster "
      "(90% identity, 80% coverage, both directions) contains at least one known resistance gene. "
      "Anchoring is deliberately at cluster level, so the annotation carries over to copies the "
      "annotation missed.\n")
    a("## The numbers\n")
    a("| | |")
    a("|---|---|")
    a(f"| genomes | {n_genomes:,} |")
    a(f"| protein clusters | {len(ann):,} |")
    a(f"| clusters with resistance evidence | {len(amr_cl):,} |")
    a(f"| islands | {len(isl):,} |")
    a(f"| genes inside an island | {len(gt):,} |")
    a(f"| island types | {len(types):,} |")
    a(f"| islands with a full +/-10 window | {int(isl['window_complete'].sum()):,} "
      f"({100*isl['window_complete'].mean():.1f}%) |")
    a(f"| islands carrying a mobile element | {int((isl['n_mge_genes']>0).sum()):,} |")
    a("")
    a("| island tier | islands | meaning |")
    a("|---|---|---|")
    for t, lbl in [("tier1_multi_acquired", "two or more distinct acquired resistance clusters"),
                   ("tier2_single_acquired", "one acquired resistance cluster"),
                   ("tier3_core_only", "only chromosomal or point-mutation anchors - the backbone control")]:
        a(f"| `{t}` | {int((isl['island_tier']==t).sum()):,} | {lbl} |")
    a("")
    a("## Read this before interpreting\n")
    a(f"- **Most islands are truncated.** {100*trunc:.0f}% run off a contig edge, and the two "
      "groups differ sharply: backbone islands mostly have a complete window, acquired-resistance "
      "islands mostly do not. That is not a bug, it is the biology of the assemblies - resistance "
      "genes sit in repeat-rich, insertion-sequence-rich regions that short reads assemble badly. "
      f"For any claim about island size or composition use the clean set: tier1 with a complete "
      f"window, {len(clean):,} islands.")
    a("- **A cluster is the unit of 'a different gene', not a gene-family label.** One cluster "
      "carries several labels for the same gene (`aac(6')-Ib`, `-Ib'`, `-Ib-cr`, `-Ib11`), and "
      "tiering on labels inflated the multi-resistance count by 1.4x before this was fixed.")
    a("- **Recurrence is counted in SNP clusters, not genomes.** Clonal outbreaks repeat the same "
      "island hundreds of times; `clonality_ratio` in `island_types.tsv` shows the inflation.")
    a("- **The qnrB cluster is mixed.** At 90% identity the plasmid-borne qnrB alleles cluster with "
      "the chromosomal qnrB intrinsic to *C. werkmanii* and *C. freundii*. No clustering at this "
      "threshold separates them; it needs allele-level assignment.")
    a("- **No enrichment statistics are reported.** Every anchor is inside an island by "
      "construction, so a naive test against all islands would be circular.\n")
    a("## Files\n")
    a("```")
    a("islands.tsv                  one row per island: geometry, truncation, tier, AMR and mobile content")
    a("island_gene_table.tsv        one row per gene in an island - the drawing layer, 55 columns")
    a("gene_island.tsv              compact gene -> island membership")
    a("island_types.tsv             recurrent architectures, counted in SNP clusters")
    a("island_type_membership.tsv   island -> type, with the architecture strings")
    a("amr_family_cooccurrence.tsv  which resistance genes travel together")
    a("anchor_clusters.tsv          the clusters that anchor, and their tier")
    a("anchor_prevalence_sweep.tsv  the calibration behind the core/mobile threshold")
    a("hit_gene_link.tsv            every MicroBIGG-E hit mapped to its coding sequence")
    a("figures/islands/             per-island gene diagrams (PNG + SVG)")
    a("figures/summary/             F01-F04, the figures that justify the tiers")
    a("figures/figure_manifest.tsv  which island each figure shows, and why it was chosen")
    a("../Pangenome/                gene index, protein clusters, cluster annotation")
    a("```\n")
    a("## Colours in the figures\n")
    a("| colour | meaning |")
    a("|---|---|")
    a("| red | acquired resistance gene |")
    a("| rose | resistance-annotated but chromosomal backbone |")
    a("| dark red | resistance by point mutation in an essential gene |")
    a("| blue | metal / biocide resistance |")
    a("| purple | virulence |")
    a("| yellow | mobile element (pale = pseudogene, a disrupted transposase) |")
    a("| grey / white | other / hypothetical |")
    a("")
    a("The grey strip under each diagram is how common each gene's cluster is across all genomes: "
      "dark means core. A red arrow over a pale strip is cargo; a rose arrow over a dark strip is "
      "backbone.\n")
    a("## Checks\n")
    a(f"{len(results) - n_fail}/{len(results)} passed at the last run.\n")
    for name, ok, detail in results:
        a(f"- {'PASS' if ok else 'FAIL'} **{name}** - {detail}")
    (ISLANDS / "README.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"\nwrote {(ISLANDS / 'README.md').relative_to(ROOT)}")
    print(f"done in {time.time()-t0:.0f}s")
    raise SystemExit(1 if n_fail else 0)


if __name__ == "__main__":
    main()
