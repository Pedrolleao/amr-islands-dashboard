"""17 - build the resistance islands.

DEFINITION (the field standard, Doron 2018, and the same rule and merge code as the DefenseFinder
precedent): an island is the union of +/-K ORF windows around anchors on ONE contig, with
overlapping or touching windows merged. K = 10. Islands never span contigs, because gene order does
not exist across an assembly gap.

An anchor is any gene whose protein cluster is AMR-associated (16). That is cluster-level by design:
one annotated hit makes every member of the cluster an anchor, in every genome, which propagates the
annotation to copies the annotation missed - measured here, 3,203 gyrA hits become 6,055 anchors.

TRUNCATION IS LABELLED, NEVER DROPPED. Most of these assemblies are draft, median 84 contigs, so a
+/-10 ORF window frequently runs off a contig edge. `n_flank_lo`, `n_flank_hi`, `window_complete`,
`frac_window` and `truncation` record exactly how much of the requested window existed. A 6-gene
contig fragment and a real 21-gene neighbourhood must not be indistinguishable downstream.

TIER separates cargo from backbone, using the anchor tiers from 16:

    tier1_multi_acquired   >= 2 distinct acquired AMR gene FAMILIES among cargo anchors
    tier2_single_acquired  exactly 1
    tier3_core_only        no cargo anchor at all - gyrA, chromosomal ampC and friends

THE UNIT OF "A DIFFERENT GENE" IS A CLUSTER, NOT A LABEL. Counting the gene_family strings a cluster
carries looks equivalent and is not: one cluster routinely carries several MicroBIGG-E family labels
for what is one gene - aac(6')-Ib, aac(6')-Ib', aac(6')-Ib-cr and aac(6')-Ib11 are the same
aminoglycoside acetyltransferase, and bla/blaTEM the same beta-lactamase. Tiering on labels promoted
7,358 single-gene islands to "multi-acquired" here, a 1.42x inflation of the headline count. The
tier therefore counts DISTINCT CARGO CLUSTERS; the family labels are kept as description.

Nothing is filtered. tier3 islands are kept and are the declared backbone control: if the tier label
is doing its job, a tier3 island drawn next to a tier1 island looks obviously different.

Usage: python3 scripts/17_define_islands.py [--group G] [--k 10]
"""
import gzip
import re
import sys
import time
import zlib

import pandas as pd

from amr_common import GROUP, GROUPS, ISLANDS, PANGENOME, read_tsv, write_tsv

K = 10                     # ORFs each side of an anchor
CARGO_TIERS = {"SHELL_AMR", "MOBILE_AMR"}
CORE_TIERS = {"CORE_AMR", "SOFTCORE_AMR", "POINT_ONLY"}

# Mobile-element machinery, recognised from the product name. Evaluated in this order, so each gene
# gets exactly one kind. `integrase_recombinase` is deliberately loose (it catches chromosomal
# XerC/D too), which is why `has_integron_integrase` is kept as its own strict flag.
MGE_PATTERNS = [
    ("integron_integrase", re.compile(r"class [12] integron integrase|integron integrase|\bintI[123]\b", re.I)),
    ("transposase", re.compile(r"transposase|\btnpA\b|Tn3 family", re.I)),
    ("IS_element", re.compile(r"\bIS\d{1,4}[A-Za-z-]*|insertion sequence|IS\d+-like element", re.I)),
    ("integrase_recombinase", re.compile(r"site-specific integrase|tyrosine-type recombinase|"
                                         r"\bintegrase\b|recombinase family|resolvase|\btnpR\b", re.I)),
    ("conjugation", re.compile(r"conjugal transfer|conjugative|relaxase|\bMob[A-Z]\b|"
                               r"type IV secretion|\btra[A-Z]\b", re.I)),
    ("plasmid_replication", re.compile(r"plasmid replication|replication initiat|\bRepA\b|RepB family", re.I)),
]
GENE_COLS = ["gene_id", "target_acc", "download_acc", "contig", "contig_ordinal", "ordinal",
             "start", "end", "strand", "locus_tag", "protein_id", "gene_symbol", "product",
             "aa_len", "aa_len_gff", "protein_md5", "is_pseudo", "is_partial", "n_segments",
             "contig_n_genes", "contig_len", "contig_topology", "contig_molecule"]


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def mge_kind(product):
    for kind, rx in MGE_PATTERNS:
        if rx.search(product):
            return kind
    return ""


def merge(intervals):
    """Merge [lo, hi] windows that overlap or touch - the precedent's rule, unchanged."""
    out = []
    for lo, hi in sorted(intervals):
        if out and lo <= out[-1][1] + 1:
            out[-1][1] = max(out[-1][1], hi)
        else:
            out.append([lo, hi])
    return out


class BlockReader:
    """Random access to one genome's block of a gzip file written as one member per genome."""

    def __init__(self, path, offsets_path):
        self.fh = open(path, "rb")
        off = read_tsv(offsets_path)
        self.offsets = {r.download_acc: (int(r.byte_offset), int(r.byte_length))
                        for r in off.itertuples()}

    def block(self, acc):
        if acc not in self.offsets:
            return []
        o, n = self.offsets[acc]
        self.fh.seek(o)
        return zlib.decompressobj(31).decompress(self.fh.read(n)).decode().rstrip("\n").split("\n")


def main():
    argv = sys.argv[1:]
    k = int(argv[argv.index("--k") + 1]) if "--k" in argv else int(
        GROUPS[GROUP].get("island_k", K))
    prefix = GROUPS[GROUP].get("island_prefix", GROUPS[GROUP].get("cluster_prefix", "ISL") + "I")

    anchors = read_tsv(ISLANDS / "anchor_clusters.tsv")
    tier_of = dict(zip(anchors["cluster_id"], anchors["anchor_tier"]))
    fam_of = {r.cluster_id: set(filter(None, r.amr_gene_families.split(";")))
              for r in anchors.itertuples()}
    cls_of = {r.cluster_id: set(filter(None, r.amr_classes.split(";")))
              for r in anchors.itertuples()}
    log(f"{len(tier_of):,} anchor clusters, K={k}")

    ann = read_tsv(PANGENOME / "cluster_annotation.tsv")
    prev_of = dict(zip(ann["cluster_id"], ann["prevalence_eff"]))
    etype_of = dict(zip(ann["cluster_id"], ann["element_type"]))
    del ann

    link = read_tsv(ISLANDS / "hit_gene_link.tsv")
    hits_by_gene = {}
    for r in link.itertuples():
        if r.gene_id:
            hits_by_gene.setdefault(r.gene_id, []).append(
                (r.element_type, r.subtype, r.symbol, r.gene_family, r.match_tier))
    del link

    man = read_tsv(ISLANDS.parent / "Genomes" / "manifest" / "genome_manifest.tsv")
    meta = {r.download_acc: (r.target_acc, r.scientific_name, r.PDS_acc) for r in man.itertuples()}

    genes = BlockReader(PANGENOME / "gene_index.tsv.gz", PANGENOME / "gene_index_offsets.tsv")
    cmap = BlockReader(PANGENOME / "gene_cluster_map.tsv.gz",
                       PANGENOME / "gene_cluster_map_offsets.tsv")

    isl_rows, gene_rows = [], []
    n_isl = 0
    t0 = time.time()
    col = {c: i for i, c in enumerate(GENE_COLS)}

    for gi, acc in enumerate(genes.offsets, 1):
        rows = [r.split("\t") for r in genes.block(acc)]
        cl = dict(line.split("\t", 1) for line in cmap.block(acc) if line)
        target_acc, species, pds = meta.get(acc, ("", "", ""))

        by_contig = {}
        for r in rows:
            by_contig.setdefault(r[col["contig"]], []).append(r)

        for contig, lst in by_contig.items():
            lst.sort(key=lambda r: int(r[col["ordinal"]]))
            n_genes_contig = len(lst)
            anchor_idx = []
            for i, r in enumerate(lst):
                cid = cl.get(r[col["gene_id"]], "")
                if cid and cid in tier_of:
                    anchor_idx.append(i)
            if not anchor_idx:
                continue

            for lo, hi in merge([[i - k, i + k] for i in anchor_idx]):
                lo_c, hi_c = max(lo, 0), min(hi, n_genes_contig - 1)
                window = lst[lo_c:hi_c + 1]
                n_isl += 1
                island_id = f"{prefix}{n_isl:06d}"
                a_lo = min(i for i in anchor_idx if lo_c <= i <= hi_c)
                a_hi = max(i for i in anchor_idx if lo_c <= i <= hi_c)
                n_flank_lo, n_flank_hi = a_lo - lo_c, hi_c - a_hi
                requested = (a_hi + k) - (a_lo - k) + 1
                trunc = ("both" if n_flank_lo < k and n_flank_hi < k else
                         "lo" if n_flank_lo < k else "hi" if n_flank_hi < k else "none")

                cargo_fams, core_fams, classes = set(), set(), set()
                cargo_clusters = set()
                n_a = {"cargo": 0, "core": 0, "point": 0}
                n_mge = n_mge_pseudo = n_pseudo = n_hypo = 0
                mge_kinds, has_intI = set(), 0
                n_stress = n_virulence = 0
                stress_syms, virulence_syms, direct_syms = set(), set(), set()
                clusters_seen = set()

                for pos, r in enumerate(window):
                    gene_id = r[col["gene_id"]]
                    cid = cl.get(gene_id, "")
                    product = r[col["product"]]
                    kind = mge_kind(product)
                    is_pseudo = r[col["is_pseudo"]] == "1"
                    tier = tier_of.get(cid, "")
                    if cid:
                        clusters_seen.add(cid)
                    if tier:
                        if tier in CARGO_TIERS:
                            n_a["cargo"] += 1
                            cargo_fams |= fam_of.get(cid, set())
                            cargo_clusters.add(cid)
                        elif tier == "POINT_ONLY":
                            n_a["point"] += 1
                        else:
                            n_a["core"] += 1
                            core_fams |= fam_of.get(cid, set())
                        classes |= cls_of.get(cid, set())
                    et = etype_of.get(cid, "")
                    if "STRESS" in et:
                        n_stress += 1
                    if "VIRULENCE" in et:
                        n_virulence += 1
                    if kind:
                        n_mge += 1
                        mge_kinds.add(kind)
                        n_mge_pseudo += is_pseudo
                        has_intI |= int(kind == "integron_integrase")
                    n_pseudo += is_pseudo
                    n_hypo += int("hypothetical" in product.lower())
                    for et2, sub, sym, fam, mt in hits_by_gene.get(gene_id, ()):
                        if et2 == "AMR":
                            direct_syms.add(sym)
                        elif et2 == "STRESS":
                            stress_syms.add(sym)
                        else:
                            virulence_syms.add(sym)
                    gene_rows.append((island_id, pos, gene_id, target_acc, acc, contig,
                                      int(r[col["ordinal"]]), cid, tier, kind,
                                      1 if tier else 0, int(is_pseudo)))

                tier_isl = ("tier1_multi_acquired" if len(cargo_clusters) >= 2 else
                            "tier2_single_acquired" if len(cargo_clusters) == 1 else
                            "tier3_core_only")
                isl_rows.append({
                    "island_id": island_id, "target_acc": target_acc, "download_acc": acc,
                    "scientific_name": species, "PDS_acc": pds, "contig": contig,
                    "contig_n_genes": n_genes_contig,
                    "contig_molecule": window[0][col["contig_molecule"]],
                    "ord_lo": int(window[0][col["ordinal"]]),
                    "ord_hi": int(window[-1][col["ordinal"]]),
                    "n_genes": len(window),
                    "bp_lo": int(window[0][col["start"]]), "bp_hi": int(window[-1][col["end"]]),
                    "span_kb": round((int(window[-1][col["end"]]) - int(window[0][col["start"]]) + 1) / 1000, 2),
                    "n_flank_lo": n_flank_lo, "n_flank_hi": n_flank_hi,
                    "window_complete": int(n_flank_lo >= k and n_flank_hi >= k),
                    "frac_window": round(len(window) / requested, 3),
                    "truncation": trunc,
                    "island_tier": tier_isl,
                    "n_anchor": n_a["cargo"] + n_a["core"] + n_a["point"],
                    "n_anchor_cargo": n_a["cargo"], "n_anchor_core": n_a["core"],
                    "n_anchor_point": n_a["point"],
                    "anchor_profile": ("cargo_only" if n_a["cargo"] and not (n_a["core"] + n_a["point"])
                                       else "mixed" if n_a["cargo"] else
                                       "point_only" if n_a["point"] and not n_a["core"] else "core_only"),
                    "n_cargo_clusters": len(cargo_clusters),
                    "cargo_cluster_ids": ";".join(sorted(cargo_clusters)),
                    "amr_gene_families": ";".join(sorted(cargo_fams)),
                    "n_amr_families": len(cargo_fams),
                    "core_amr_families": ";".join(sorted(core_fams)),
                    "amr_classes": ";".join(sorted(classes)), "n_amr_classes": len(classes),
                    "direct_hit_symbols": ";".join(sorted(direct_syms)),
                    "n_stress_genes": n_stress, "stress_symbols": ";".join(sorted(stress_syms)),
                    "n_virulence_genes": n_virulence,
                    "virulence_symbols": ";".join(sorted(virulence_syms)),
                    "n_mge_genes": n_mge, "frac_mge": round(n_mge / len(window), 3),
                    "n_mge_pseudo": n_mge_pseudo, "mge_kinds": ";".join(sorted(mge_kinds)),
                    "has_integron_integrase": has_intI,
                    "n_pseudo": n_pseudo, "frac_hypothetical": round(n_hypo / len(window), 3),
                    "n_clusters": len(clusters_seen),
                })
        if gi % 1000 == 0:
            log(f"  {gi:,} genomes, {n_isl:,} islands ({time.time()-t0:.0f}s)")

    isl = pd.DataFrame(isl_rows)
    write_tsv(isl, "islands.tsv", folder=ISLANDS)
    gi_df = pd.DataFrame(gene_rows, columns=["island_id", "pos_in_island", "gene_id", "target_acc",
                                             "download_acc", "contig", "ordinal", "cluster_id",
                                             "anchor_tier", "mge_kind", "is_anchor", "is_pseudo"])
    write_tsv(gi_df, "gene_island.tsv", folder=ISLANDS)

    # ---- QC ------------------------------------------------------------------------------------
    print("\n===== QC 17_define_islands =====")
    print(f"islands: {len(isl):,}   in-island genes: {len(gi_df):,}   "
          f"genomes with >=1 island: {isl['download_acc'].nunique():,}")
    print(f"mean genes per island: {isl['n_genes'].mean():.1f}   max: {isl['n_genes'].max()}")
    bad = (isl["ord_hi"] - isl["ord_lo"] + 1 != isl["n_genes"]).sum()
    print(f"ordinal contiguity violations: {bad}  [must be 0]")
    print(f"sum(n_genes) == rows in gene_island: {isl['n_genes'].sum() == len(gi_df)}")
    print("\nisland tiers:")
    for t, n in isl["island_tier"].value_counts().items():
        print(f"  {t:22s} {n:8,} ({100*n/len(isl):5.1f}%)")
    print("\nanchor profile:")
    for t, n in isl["anchor_profile"].value_counts().items():
        print(f"  {t:12s} {n:8,}")
    print("\ntruncation:")
    for t, n in isl["truncation"].value_counts().items():
        print(f"  {t:6s} {n:8,} ({100*n/len(isl):5.1f}%)")
    print(f"window_complete: {isl['window_complete'].sum():,} "
          f"({100*isl['window_complete'].mean():.1f}%)   "
          f"median frac_window {isl['frac_window'].median():.3f}")
    t1 = isl[isl["island_tier"] == "tier1_multi_acquired"]
    print(f"\ntier1 islands: {len(t1):,}; median distinct cargo clusters "
          f"{pd.to_numeric(t1['n_cargo_clusters']).median():.0f}, "
          f"max {pd.to_numeric(isl['n_cargo_clusters']).max()}")
    print(f"  (counting family LABELS instead would give "
          f"{(pd.to_numeric(isl['n_amr_families'])>=2).sum():,} - the inflation this avoids)")
    print(f"islands with a mobile element: {(isl['n_mge_genes']>0).sum():,} "
          f"({100*(isl['n_mge_genes']>0).mean():.1f}%);  with an integron integrase: "
          f"{int(isl['has_integron_integrase'].sum()):,}")
    print(f"islands on a plasmid-labelled contig: {(isl['contig_molecule']=='plasmid').sum():,}")
    log(f"done in {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
