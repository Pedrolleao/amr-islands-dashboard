"""20 - draw the islands: gene-arrow diagrams, annotated per protein cluster.

Uses dna_features_viewer (installed in the base python; R has no graphics libraries on this machine).

WHAT A FIGURE SHOWS
  - one arrow per gene, in real contig coordinates, pointing along the strand
  - colour = `draw_role` from 18, the single colour dimension shared by every view:
      red        acquired AMR gene (the point of the figure)
      rose       AMR-annotated but chromosomal backbone (ampC, oqxAB ...)
      dark red   resistance by point mutation in an essential gene (gyrA)
      blue       metal / biocide resistance (stress)        purple  virulence
      yellow     mobile-element machinery (pale = pseudogene, i.e. a disrupted transposase)
      grey/white other / hypothetical
  - a prevalence strip under the arrows: how common each gene's protein cluster is across the 6,061
    genomes, dark = core. This is what makes the tier label visible by eye - a red arrow sitting in
    a pale (rare) neighbourhood is cargo; a rose arrow in a dark band is backbone.
  - contig edges drawn as a dotted line with the flank actually available, because most of these
    islands are truncated and a 3-gene contig fragment must not read as a real neighbourhood.

WHICH ISLANDS. There are 66,005; rendering all of them would be a way of not choosing. The selection
is deterministic and recorded in `figure_manifest.tsv` with the reason each island was picked.

Usage: python3 scripts/20_island_figures.py [--group G] [--only ISLAND_ID] [--per-class N]
"""
import sys
import time

import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from dna_features_viewer import GraphicFeature, GraphicRecord

from amr_common import DRAW_ROLE_COLOURS, ISLANDS, read_tsv, write_tsv

SEED = 20260921
MAX_GENES = 110           # a 60-gene multi-resistance region is the most interesting figure
                          # in the set, so the cap has to be well above the typical island
CONTROLS = [("PDT000050772.2", "NZ_JYFY02000004.1")]


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def colour_of(role, is_pseudo):
    if role == "mge" and is_pseudo:
        return DRAW_ROLE_COLOURS["mge_pseudo"]
    return DRAW_ROLE_COLOURS.get(role, DRAW_ROLE_COLOURS["other"])


def draw(island, genes, path, cluster_labels=True):
    genes = genes.sort_values("start_i")
    x0 = genes["start_i"].min()
    span = genes["end_i"].max() - x0 + 1
    feats = []
    for g in genes.itertuples():
        # label what carries meaning: never "hypothetical protein" (it is noise repeated ten
        # times per island), always AMR / stress / virulence / mobile elements, and ordinary genes
        # only when the island is small enough for the labels not to collide.
        if g.draw_role == "hypothetical":
            label = None
        elif g.draw_role == "other":
            label = g.label_text if len(genes) <= 22 else None
        else:
            label = g.label_text
        feats.append(GraphicFeature(
            start=g.start_i - x0, end=g.end_i - x0, strand=1 if g.strand == "+" else -1,
            color=colour_of(g.draw_role, str(g.is_pseudo) == "1"), label=label,
            linewidth=0.6, linecolor="#5d6d7e"))

    rec = GraphicRecord(sequence_length=span, features=feats)
    fig, (ax, ax2) = plt.subplots(
        2, 1, figsize=(max(9, min(0.75 * len(genes), 34)), 4.2),
        gridspec_kw={"height_ratios": [5, 1]})
    rec.plot(ax=ax, figure_width=None, annotate_inline=False)

    # prevalence strip: one cell per gene, dark = the cluster is near-universal
    for g in genes.itertuples():
        prev = float(g.pc_prevalence_eff) if g.pc_prevalence_eff not in ("", None) else 0.0
        ax2.add_patch(Rectangle((g.start_i - x0, 0), max(g.end_i - g.start_i, 1), 1,
                                facecolor=str(1 - 0.85 * prev), edgecolor="white", linewidth=0.3))
    ax2.set_xlim(0, span)
    ax2.set_ylim(0, 1)
    ax2.set_yticks([])
    ax2.set_xticks([])
    ax2.set_ylabel("cluster\nprevalence", fontsize=6, rotation=0, ha="right", va="center")
    for s in ax2.spines.values():
        s.set_visible(False)

    # contig edges
    if island["truncation"] in ("lo", "both"):
        ax.axvline(0, color="k", ls=":", lw=1)
        ax.text(0, ax.get_ylim()[1] * 0.92, f" contig start (flank {island['n_flank_lo']}/10)",
                fontsize=6, va="top")
    if island["truncation"] in ("hi", "both"):
        ax.axvline(span, color="k", ls=":", lw=1)
        ax.text(span, ax.get_ylim()[1] * 0.92, f"contig end (flank {island['n_flank_hi']}/10) ",
                fontsize=6, va="top", ha="right")

    mol = island.get("contig_molecule") or "contig"
    ax.set_title(
        f"{island['island_id']} · {island['target_acc']} · {island['scientific_name']}\n"
        f"{island['contig']} ({mol}) · {island['island_tier']} · "
        f"{island['n_genes']} genes / {island['span_kb']} kb · "
        f"trunc={island['truncation']} frac={island['frac_window']} · "
        f"AMR: {island['amr_gene_families'][:70] or '-'}",
        fontsize=8, loc="left")
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    if str(path).endswith(".png"):
        fig.savefig(str(path)[:-4] + ".svg", bbox_inches="tight")
    plt.close(fig)


def main():
    argv = sys.argv[1:]
    only = argv[argv.index("--only") + 1] if "--only" in argv else None
    per_class = int(argv[argv.index("--per-class") + 1]) if "--per-class" in argv else 20

    isl = read_tsv(ISLANDS / "islands.tsv")
    for c in ("n_genes", "n_amr_families", "n_mge_genes", "window_complete", "n_flank_lo",
              "n_flank_hi", "has_integron_integrase", "n_stress_genes"):
        isl[c] = pd.to_numeric(isl[c])
    gt = read_tsv(ISLANDS / "island_gene_table.tsv")
    gt["start_i"] = pd.to_numeric(gt["start"])
    gt["end_i"] = pd.to_numeric(gt["end"])
    log(f"{len(isl):,} islands, {len(gt):,} genes")

    picks = []                                        # (island_id, reason, rank)
    if only:
        picks = [(only, "requested", 0)]
    else:
        for target, contig in CONTROLS:
            sub = isl[(isl["target_acc"] == target) & (isl["contig"] == contig)]
            for n, r in enumerate(sub.itertuples(), 1):
                picks.append((r.island_id, "control", n))
        complete = isl[isl["window_complete"] == 1]
        t1 = complete[complete["island_tier"] == "tier1_multi_acquired"]
        for n, r in enumerate(t1.nlargest(per_class, "n_amr_families").itertuples(), 1):
            picks.append((r.island_id, "most_amr_families", n))
        comp = t1[t1["n_amr_families"] >= 2].nlargest(per_class, "n_mge_genes")
        for n, r in enumerate(comp.itertuples(), 1):
            picks.append((r.island_id, "most_mobile_elements", n))
        integron = t1[(t1["has_integron_integrase"] == 1)].nlargest(per_class, "n_amr_families")
        for n, r in enumerate(integron.itertuples(), 1):
            picks.append((r.island_id, "class1_integron", n))
        stress = t1[t1["n_stress_genes"] >= 2].nlargest(10, "n_stress_genes")
        for n, r in enumerate(stress.itertuples(), 1):
            picks.append((r.island_id, "amr_plus_stress", n))
        core = complete[complete["island_tier"] == "tier3_core_only"]
        for n, r in enumerate(core.sample(min(10, len(core)), random_state=SEED).itertuples(), 1):
            picks.append((r.island_id, "backbone_control", n))   # drawn on purpose, to compare
        rnd = t1.sample(min(15, len(t1)), random_state=SEED)
        for n, r in enumerate(rnd.itertuples(), 1):
            picks.append((r.island_id, "random_tier1", n))

    seen, ordered = set(), []
    for iid, reason, rank in picks:
        if iid not in seen:
            seen.add(iid)
            ordered.append((iid, reason, rank))

    out_dir = ISLANDS / "figures" / "islands"
    out_dir.mkdir(parents=True, exist_ok=True)
    isl_by_id = isl.set_index("island_id", drop=False)
    rows, skipped = [], 0
    t0 = time.time()
    for n, (iid, reason, rank) in enumerate(ordered, 1):
        if iid not in isl_by_id.index:
            continue
        row = isl_by_id.loc[iid]
        genes = gt[gt["island_id"] == iid]
        if len(genes) == 0 or len(genes) > MAX_GENES:
            skipped += 1
            continue
        path = out_dir / f"{iid}__{row['target_acc']}__{row['island_tier']}.png"
        draw(row, genes, path)
        rows.append({"island_id": iid, "png": path.name, "svg": path.name[:-4] + ".svg",
                     "reason": reason, "rank": rank, "target_acc": row["target_acc"],
                     "scientific_name": row["scientific_name"], "contig": row["contig"],
                     "island_tier": row["island_tier"], "n_genes": row["n_genes"],
                     "n_amr_families": row["n_amr_families"],
                     "amr_gene_families": row["amr_gene_families"],
                     "n_mge_genes": row["n_mge_genes"], "truncation": row["truncation"],
                     "PDS_acc": row["PDS_acc"]})
        if n % 20 == 0:
            log(f"  {n}/{len(ordered)} ({time.time()-t0:.0f}s)")
    write_tsv(pd.DataFrame(rows), "figure_manifest.tsv", folder=ISLANDS / "figures")

    print("\n===== QC 20_island_figures =====")
    print(f"figures drawn: {len(rows)}   skipped (empty or > {MAX_GENES} genes): {skipped}")
    if rows:
        man = pd.DataFrame(rows)
        for r, n in man["reason"].value_counts().items():
            print(f"  {r:22s} {n}")
        print(f"\noutput: {out_dir.relative_to(ISLANDS.parent.parent)}  (PNG + SVG each)")
    log(f"done in {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
