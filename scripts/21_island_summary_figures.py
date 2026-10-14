"""21 - the four summary figures that say whether the island set can be trusted.

F01  island tiers and how much of each window actually existed. The headline caveat of this dataset
     in one picture: acquired resistance sits on short, fragmented contigs, backbone does not.
F02  prevalence of the AMR-associated protein clusters, with the core threshold drawn on it. This is
     the figure that justifies the anchor tiers rather than asserting them.
F03  the most recurrent island types that actually CARRY acquired resistance, genomes against SNP
     clusters. Backbone types (gyrA and other chromosomal neighbourhoods) are excluded: they are the
     most repeated architectures in the dataset simply because those genes are in every genome, and
     they would otherwise fill the top of a figure about resistance islands. The gap between the two
     bars is clonal expansion - counting genomes would rank one outbreak as the commonest island.
F04  which resistance gene families travel together, counted in distinct SNP clusters.

Usage: python3 scripts/21_island_summary_figures.py [--group G]
"""
import time

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from amr_common import ISLANDS, PANGENOME, read_tsv

TIER_COLOURS = {"tier1_multi_acquired": "#c0392b", "tier2_single_acquired": "#e59866",
                "tier3_core_only": "#85929e"}


def save(fig, name, out):
    fig.tight_layout()
    fig.savefig(out / f"{name}.png", dpi=150, bbox_inches="tight")
    fig.savefig(out / f"{name}.svg", bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {name}.png/.svg")


def main():
    t0 = time.time()
    out = ISLANDS / "figures" / "summary"
    out.mkdir(parents=True, exist_ok=True)
    isl = read_tsv(ISLANDS / "islands.tsv")
    for c in ("n_genes", "frac_window", "window_complete", "n_cargo_clusters", "contig_n_genes"):
        isl[c] = pd.to_numeric(isl[c])

    # ---- F01 tiers and truncation ----------------------------------------------------------------
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4))
    counts = isl["island_tier"].value_counts()
    a1.bar(range(len(counts)), counts.values,
           color=[TIER_COLOURS.get(t, "#bdc3c7") for t in counts.index])
    a1.set_xticks(range(len(counts)))
    a1.set_xticklabels([t.replace("_", "\n") for t in counts.index], fontsize=8)
    for i, v in enumerate(counts.values):
        a1.text(i, v, f"{v:,}", ha="center", va="bottom", fontsize=8)
    a1.set_ylabel("islands")
    a1.set_title("Island tiers", loc="left", fontsize=10)

    for t in counts.index:
        sub = isl[isl["island_tier"] == t]
        a2.hist(sub["frac_window"], bins=25, alpha=0.6, label=f"{t} (n={len(sub):,})",
                color=TIER_COLOURS.get(t, "#bdc3c7"))
    a2.set_xlabel("fraction of the requested ±10-gene window that existed")
    a2.set_ylabel("islands")
    a2.legend(fontsize=7)
    a2.set_title("Acquired resistance sits on fragmented contigs", loc="left", fontsize=10)
    save(fig, "F01_tiers_and_truncation", out)

    # ---- F02 prevalence of AMR clusters -----------------------------------------------------------
    ann = read_tsv(PANGENOME / "cluster_annotation.tsv")
    ann["prevalence_eff"] = pd.to_numeric(ann["prevalence_eff"], errors="coerce").fillna(0)
    ann["n_amr_hits"] = pd.to_numeric(ann["n_amr_hits"])
    ann["n_point_hits"] = pd.to_numeric(ann["n_point_hits"])
    amr = ann[(ann["n_amr_hits"] > 0) | (ann["n_point_hits"] > 0)]
    fig, ax = plt.subplots(figsize=(9, 4.2))
    ax.hist(amr["prevalence_eff"], bins=40, color="#5d6d7e")
    ax.axvline(0.95, color="#c0392b", ls="--", lw=1.2)
    ax.text(0.95, ax.get_ylim()[1] * 0.9, " core threshold 0.95", color="#c0392b", fontsize=8)
    for sym, lbl in [("blaCMY", "blaCMY (chromosomal AmpC)"), ("sul1", "sul1 (acquired)"),
                     ("blaKPC-2", "blaKPC-2"), ("blaNDM-1", "blaNDM-1")]:
        hit = amr[amr["amr_symbols"].str.contains(sym, regex=False, na=False)]
        if len(hit):
            x = hit["prevalence_eff"].max()
            ax.annotate(lbl, xy=(x, 2), xytext=(x, ax.get_ylim()[1] * (0.35 + 0.12 * hash(sym) % 4)),
                        fontsize=7, ha="center",
                        arrowprops=dict(arrowstyle="-", lw=0.6, color="#566573"))
    gy = amr[amr["point_symbols"].str.contains("gyrA", regex=False, na=False)]
    if len(gy):
        ax.annotate("gyrA (point mutation,\nin every genome)",
                    xy=(gy["prevalence_eff"].max(), 2),
                    xytext=(0.78, ax.get_ylim()[1] * 0.62), fontsize=7, ha="center",
                    arrowprops=dict(arrowstyle="-", lw=0.6, color="#566573"))
    ax.set_xlabel("prevalence of the protein cluster (max of genus-wide and within-species)")
    ax.set_ylabel("AMR-associated clusters")
    ax.set_title(f"Why the anchor tiers exist ({len(amr):,} AMR-associated clusters)",
                 loc="left", fontsize=10)
    save(fig, "F02_amr_cluster_prevalence", out)

    # ---- F03 recurrence: genomes vs SNP clusters ---------------------------------------------------
    types = read_tsv(ISLANDS / "island_types.tsv")
    for c in ("n_genomes", "n_pds", "n_islands"):
        types[c] = pd.to_numeric(types[c])
    cargo = types[types["amr_gene_families"] != ""]
    n_backbone = len(types) - len(cargo)
    top = cargo.nlargest(20, "n_pds")
    fig, ax = plt.subplots(figsize=(10, 5.5))
    y = np.arange(len(top))
    ax.barh(y + 0.2, top["n_genomes"], height=0.4, color="#aeb6bf", label="genomes")
    ax.barh(y - 0.2, top["n_pds"], height=0.4, color="#c0392b", label="distinct SNP clusters")
    ax.set_yticks(y)
    ax.set_yticklabels([f"{r.island_type_id}  {(r.amr_gene_families or '-')[:34]}"
                        for r in top.itertuples()], fontsize=7)
    ax.invert_yaxis()
    ax.set_xlabel("count")
    ax.legend(fontsize=8)
    ax.set_title("Most recurrent island types carrying acquired resistance\n"
                 f"(the bar gap is clonal expansion; {n_backbone:,} backbone types excluded)",
                 loc="left", fontsize=10)
    save(fig, "F03_type_recurrence", out)

    # ---- F04 co-occurrence -------------------------------------------------------------------------
    co = read_tsv(ISLANDS / "amr_family_cooccurrence.tsv")
    if len(co):
        co["n_pds"] = pd.to_numeric(co["n_pds"])
        top_co = co.nlargest(25, "n_pds")
        fams = sorted({f for r in top_co.itertuples() for f in (r.family_a, r.family_b)})
        idx = {f: i for i, f in enumerate(fams)}
        mat = np.zeros((len(fams), len(fams)))
        for r in top_co.itertuples():
            mat[idx[r.family_a], idx[r.family_b]] = r.n_pds
            mat[idx[r.family_b], idx[r.family_a]] = r.n_pds
        fig, ax = plt.subplots(figsize=(8.5, 7.5))
        im = ax.imshow(mat, cmap="Reds")
        ax.set_xticks(range(len(fams)))
        ax.set_xticklabels(fams, rotation=90, fontsize=7)
        ax.set_yticks(range(len(fams)))
        ax.set_yticklabels(fams, fontsize=7)
        fig.colorbar(im, ax=ax, shrink=0.7, label="distinct SNP clusters")
        ax.set_title("Resistance gene families that travel together", loc="left", fontsize=10)
        save(fig, "F04_family_cooccurrence", out)

    print(f"\n===== QC 21_island_summary_figures =====")
    print(f"figures written to {out.relative_to(ISLANDS.parent.parent)}")
    print(f"done in {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
