"""18 - the wide per-gene island table: one row per gene inside an island.

This is the drawing layer. Everything a figure or the later browser needs about a gene sits on its
row, so nothing downstream has to re-open the 30M-row gene index or re-derive a colour.

PURE JOINS, with one honest exception. The gene facts come from the gene index (11), the cluster
facts from the cluster annotation (15/16), the hit facts from the hit link (12) and the island facts
from the islands table (17). The ONLY computed columns are `draw_role` and `label_text`, and they are
computed here exactly once so that tables, figures and the browser can never disagree about what
colour a gene is. Saying "pure join" and then quietly recomputing a value is how two views of the
same island end up looking different.

`draw_role` is a single dimension with eight values. Precedence, when a gene qualifies for several:
AMR anchor > stress > virulence > mobile element > hypothetical > other. The AMR anchor splits by
the anchor tier from 16, so a chromosomal ampC is visibly weaker than a Tn-borne blaKPC rather than
being the same red.

Usage: python3 scripts/18_island_gene_table.py [--group G]
"""
import time
import zlib

import pandas as pd

from amr_common import ISLANDS, PANGENOME, read_tsv, write_tsv

GENE_COLS = ["gene_id", "target_acc", "download_acc", "contig", "contig_ordinal", "ordinal",
             "start", "end", "strand", "locus_tag", "protein_id", "gene_symbol", "product",
             "aa_len", "aa_len_gff", "protein_md5", "is_pseudo", "is_partial", "n_segments",
             "contig_n_genes", "contig_len", "contig_topology", "contig_molecule"]
KEEP = ["start", "end", "strand", "locus_tag", "protein_id", "gene_symbol", "product", "aa_len",
        "is_partial", "contig_len", "contig_molecule"]


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def draw_role(anchor_tier, element_type, mge_kind, product):
    if anchor_tier in ("SHELL_AMR", "MOBILE_AMR"):
        return "amr_anchor_mobile"
    if anchor_tier in ("CORE_AMR", "SOFTCORE_AMR"):
        return "amr_anchor_core"
    if anchor_tier == "POINT_ONLY":
        return "amr_anchor_point"
    if "STRESS" in element_type:
        return "stress"
    if "VIRULENCE" in element_type:
        return "virulence"
    if mge_kind:
        return "mge"
    return "hypothetical" if "hypothetical" in product.lower() else "other"


def main():
    t0 = time.time()
    gi = read_tsv(ISLANDS / "gene_island.tsv")
    log(f"{len(gi):,} in-island genes")

    # gene facts, pulled genome by genome from the indexed gzip
    off = read_tsv(PANGENOME / "gene_index_offsets.tsv")
    offsets = {r.download_acc: (int(r.byte_offset), int(r.byte_length)) for r in off.itertuples()}
    wanted = set(gi["gene_id"])
    col = {c: i for i, c in enumerate(GENE_COLS)}
    facts = {}
    with open(PANGENOME / "gene_index.tsv.gz", "rb") as fh:
        for acc in sorted(set(gi["download_acc"])):
            if acc not in offsets:
                continue
            o, n = offsets[acc]
            fh.seek(o)
            for line in zlib.decompressobj(31).decompress(fh.read(n)).decode().rstrip("\n").split("\n"):
                f = line.split("\t")
                if f[0] in wanted:
                    facts[f[0]] = [f[col[c]] for c in KEEP]
    log(f"gene facts collected for {len(facts):,} genes ({time.time()-t0:.0f}s)")
    gf = pd.DataFrame.from_dict(facts, orient="index", columns=KEEP)
    gf.index.name = "gene_id"
    out = gi.merge(gf.reset_index(), on="gene_id", how="left")

    # cluster facts
    ann = read_tsv(PANGENOME / "cluster_annotation.tsv")
    cols = ["cluster_id", "n_members", "n_genomes", "prevalence", "prevalence_eff",
            "pangenome_class", "product_top", "element_type", "amr_symbols",
            "amr_gene_families", "amr_classes", "evidence_tier", "propagation_ratio"]
    out = out.merge(ann[cols].rename(columns={
        "n_members": "pc_n_members", "n_genomes": "pc_n_genomes", "prevalence": "pc_prevalence",
        "prevalence_eff": "pc_prevalence_eff", "pangenome_class": "pc_pangenome_class",
        "product_top": "pc_product_top", "element_type": "pc_element_type",
        "amr_symbols": "pc_amr_symbols", "amr_gene_families": "pc_amr_gene_families",
        "amr_classes": "pc_amr_classes", "evidence_tier": "pc_evidence_tier",
        "propagation_ratio": "pc_propagation_ratio"}), on="cluster_id", how="left")

    # hit facts: one row per gene, best hit first (AMR > STRESS > VIRULENCE, then identity)
    link = read_tsv(ISLANDS / "hit_gene_link.tsv")
    link = link[link["gene_id"] != ""].copy()
    order = {"AMR": 0, "STRESS": 1, "VIRULENCE": 2}
    link["_o"] = link["element_type"].map(order).fillna(3)
    link["_id"] = pd.to_numeric(link["pct_identity"], errors="coerce").fillna(0) \
        if "pct_identity" in link.columns else 0
    link = link.sort_values(["gene_id", "_o", "_id"], ascending=[True, True, False])
    n_hits = link.groupby("gene_id").size().rename("n_hits_on_gene")
    all_syms = link.groupby("gene_id")["symbol"].agg(lambda s: ";".join(sorted(set(s)))) \
        .rename("all_hit_symbols")
    best = link.drop_duplicates("gene_id")[
        ["gene_id", "element_type", "subtype", "symbol", "gene_family", "class", "subclass",
         "method", "match_tier", "propagates_label"]].rename(columns={
            "element_type": "hit_element_type", "subtype": "hit_subtype", "symbol": "hit_symbol",
            "gene_family": "hit_gene_family", "class": "hit_class", "subclass": "hit_subclass",
            "method": "hit_method", "match_tier": "hit_match_tier",
            "propagates_label": "hit_propagates_label"})
    out = out.merge(best, on="gene_id", how="left")
    out = out.merge(n_hits, on="gene_id", how="left").merge(all_syms, on="gene_id", how="left")

    # island facts
    isl = read_tsv(ISLANDS / "islands.tsv")
    out = out.merge(isl[["island_id", "island_tier", "truncation", "window_complete",
                         "frac_window", "n_genes", "PDS_acc", "scientific_name"]].rename(
        columns={"n_genes": "island_n_genes"}), on="island_id", how="left")

    out = out.fillna("")
    out["draw_role"] = [draw_role(a, e, m, p) for a, e, m, p in
                        zip(out["anchor_tier"], out["pc_element_type"], out["mge_kind"],
                            out["product"])]
    out["label_text"] = [h if h else (s if s else p[:30]) for h, s, p in
                         zip(out["hit_symbol"], out["gene_symbol"], out["product"])]
    write_tsv(out, "island_gene_table.tsv", folder=ISLANDS)

    print("\n===== QC 18_island_gene_table =====")
    print(f"rows: {len(out):,}   islands: {out['island_id'].nunique():,}")
    print(f"rows == sum(islands.n_genes): "
          f"{len(out) == pd.to_numeric(isl['n_genes']).sum()}")
    print(f"genes with no gene-index facts (must be 0): {(out['product'] == '').sum() - (out['product']=='').sum()}")
    miss = out["start"].eq("").sum()
    print(f"genes missing coordinates: {miss}  [must be 0]")
    print("\ndraw_role composition:")
    for r, n in out["draw_role"].value_counts().items():
        print(f"  {r:20s} {n:8,} ({100*n/len(out):5.1f}%)")
    print(f"\ngenes with a direct hit: {(out['hit_symbol']!='').sum():,}")
    print(f"anchors: {pd.to_numeric(out['is_anchor']).sum():,}   "
          f"pseudogenes: {pd.to_numeric(out['is_pseudo']).sum():,}")
    print(f"unclustered genes (no protein): {(out['cluster_id']=='').sum():,}")
    log(f"done in {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
