"""12 - link every MicroBIGG-E hit to the CDS it sits on.

The hit table gives a resistance gene's coordinates; the gene index (11) gives every CDS and its
order on the contig. This joins the two, which is what turns "isolate X carries blaKPC-2" into
"blaKPC-2 is gene 63 of contig NZ_..., with these genes either side".

WHY THIS IS NOT A ONE-LINE COORDINATE JOIN. Measured on a 6,501-hit sample:

    T1_exact            identical (contig, start, stop)            90.9%
    T2_stop_codon       same strand, both ends within 3 bp          6.9%
    T3_overlap50        same strand, reciprocal overlap >= 50%      2.1%
    T4_overlap_partial  same strand, any overlap                    0.1%
    T5_none             no overlapping CDS                          ~0

Hits that carry a protein accession match exactly 99.96% of the time. The 8,027 AMR hits with NO
protein accession - AMRFinderPlus' `*X` methods, i.e. the gene was found in translated DNA because
the annotation missed it or it is pseudogenised - match exactly only 3.5% of the time. Those are
precisely the interesting ones (disrupted transposons, truncated genes), so an exact-only join would
silently discard the hardest and most informative cases.

THE CONTAMINATION TRAP, AND WHY `propagates_label` EXISTS. Anchors are defined at CLUSTER level in
16: one hit anywhere in a cluster makes every member of that cluster an anchor. A cluster can hold
6,000 members. So ONE wrong overlap match - an X-method hit landing on a neighbouring CDS that is a
different gene - could turn an entire transposase or hypothetical-protein family into an "AMR"
anchor. Measured: within T3/T4 the matched CDS agrees with the hit's own protein accession only
17.5% of the time, versus 99.95% at T1. `propagates_label` is therefore 1 for T1/T2 and 0 for T3/T4:
the weak evidence is kept, joined and drawn, but 15 counts it separately so a cluster that is only
"AMR" through weak matches is visible as `evidence_tier = PROPAGATED_ONLY`. Nothing is dropped.

POINT MUTATIONS are carried through with `subtype = POINT` and are deliberately NOT mixed into the
acquired-AMR counts downstream: `gyrA_T83I` means an essential gene has a resistance substitution,
not that a resistance element has been acquired.

Usage: python3 scripts/12_link_hits_to_genes.py [--group G]
"""
import gzip
import sys
import time
import zlib
from bisect import bisect_left

import pandas as pd

from amr_common import ISLANDS, OUT, PANGENOME, read_tsv, write_tsv

NEAR_BP = 3               # T2: how far a CDS end may differ from the hit end (the stop codon)
MIN_OVERLAP = 0.50        # T3: reciprocal overlap fraction

GENE_COLS = ["gene_id", "target_acc", "download_acc", "contig", "contig_ordinal", "ordinal",
             "start", "end", "strand", "locus_tag", "protein_id", "gene_symbol", "product",
             "aa_len", "aa_len_gff", "protein_md5", "is_pseudo", "is_partial", "n_segments",
             "contig_n_genes", "contig_len", "contig_topology", "contig_molecule"]
OUT_COLS = ["hit_id", "target_acc", "download_acc", "contig", "hit_start", "hit_stop", "hit_strand",
            "element_type", "subtype", "symbol", "gene_family", "class", "subclass", "method",
            "hit_quality", "hit_protein", "gene_id", "ordinal", "contig_n_genes", "match_tier",
            "delta_start", "delta_end", "overlap_bp", "overlap_frac_hit", "overlap_frac_cds",
            "n_candidate_cds", "protein_agrees", "gene_protein_md5", "gene_is_pseudo",
            "propagates_label"]


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


class GeneIndexReader:
    """Random access into gene_index.tsv.gz, one independent gzip member per genome."""

    def __init__(self, folder):
        self.path = folder / "gene_index.tsv.gz"
        off = read_tsv(folder / "gene_index_offsets.tsv")
        self.offsets = {r.download_acc: (int(r.byte_offset), int(r.byte_length))
                        for r in off.itertuples()}
        self.fh = open(self.path, "rb")

    def genome(self, acc):
        """-> list of dicts for one genome, or [] if it was never indexed."""
        if acc not in self.offsets:
            return []
        off, length = self.offsets[acc]
        self.fh.seek(off)
        raw = zlib.decompressobj(31).decompress(self.fh.read(length))
        return [dict(zip(GENE_COLS, line.split("\t")))
                for line in raw.decode().rstrip("\n").split("\n")]


def choose(cands, hit_start, hit_stop, hit_protein):
    """Pick one CDS among overlapping candidates, deterministically.

    Protein-accession agreement first (it is decisive when present), then overlap length, then the
    smaller start offset, then the lower ordinal so the result never depends on input order.
    """
    def key(c):
        return (0 if (hit_protein and c["protein_id"] == hit_protein) else 1,
                -c["_ovl"], abs(c["_start"] - hit_start), c["_ord"])
    return sorted(cands, key=key)[0]


def main():
    hits = read_tsv(OUT / "microbigge_hits.tsv")
    # `class` is a Python keyword, so itertuples() would rename it positionally (_14) - which
    # silently breaks the day a column is inserted. Rename it once, explicitly.
    hits = hits.rename(columns={"class": "amr_class"})
    reader = GeneIndexReader(PANGENOME)
    log(f"{len(hits):,} hits, {len(reader.offsets):,} genomes in the gene index")

    man = read_tsv((PANGENOME / "genome_index.tsv"))
    acc_of = dict(zip(man["target_acc"], man["download_acc"]))

    rows, unmatched, synthetic = [], [], []
    tiers = {}
    n_genome_absent = 0
    t0 = time.time()
    for n, (target, grp) in enumerate(hits.groupby("target_acc", sort=False), 1):
        acc = acc_of.get(target, "")
        genes = reader.genome(acc) if acc else []
        if not genes:
            n_genome_absent += len(grp)
            for h in grp.itertuples():
                unmatched.append((f"{h.target_acc}:{h.contig}:{h.start}-{h.stop}:{h.strand}:{h.symbol}",
                                  target, acc, h.contig, h.start, h.stop, h.symbol,
                                  h.element_type, "genome_absent"))
            continue

        by_contig = {}
        for g in genes:
            g["_start"], g["_end"] = int(g["start"]), int(g["end"])
            g["_ord"] = int(g["ordinal"])
            by_contig.setdefault(g["contig"], []).append(g)
        exact = {}
        for contig, lst in by_contig.items():
            lst.sort(key=lambda g: g["_start"])
            for g in lst:
                exact[(contig, g["_start"], g["_end"])] = g
        starts = {c: [g["_start"] for g in lst] for c, lst in by_contig.items()}

        for h in grp.itertuples():
            hs, he = int(h.start), int(h.stop)
            hid = f"{h.target_acc}:{h.contig}:{h.start}-{h.stop}:{h.strand}:{h.symbol}"
            lst = by_contig.get(h.contig)
            if lst is None:
                unmatched.append((hid, target, acc, h.contig, h.start, h.stop, h.symbol,
                                  h.element_type, "contig_absent"))
                tiers["T5_none"] = tiers.get("T5_none", 0) + 1
                synthetic.append((acc, target, h.contig, hs, he, h.strand, h.symbol,
                                  h.element_name, h.element_type))
                continue

            g = exact.get((h.contig, hs, he))
            tier, cands = None, []
            if g is not None:
                tier, best = "T1_exact", g
                best["_ovl"] = he - hs + 1
            else:
                # candidates: CDS overlapping the hit span (scan a small window around it)
                i = bisect_left(starts[h.contig], hs) - 3
                for g2 in lst[max(i, 0): i + 12 if i > 0 else 12]:
                    ovl = min(he, g2["_end"]) - max(hs, g2["_start"]) + 1
                    if ovl > 0:
                        g2["_ovl"] = ovl
                        cands.append(g2)
                if cands:
                    same = [c for c in cands if c["strand"] == h.strand] or cands
                    near = [c for c in same if abs(c["_start"] - hs) <= NEAR_BP
                            and abs(c["_end"] - he) <= NEAR_BP]
                    if near:
                        tier, best = "T2_stop_codon", choose(near, hs, he, h.protein)
                    else:
                        hit_len, best = he - hs + 1, choose(same, hs, he, h.protein)
                        cds_len = best["_end"] - best["_start"] + 1
                        if best["_ovl"] / hit_len >= MIN_OVERLAP and best["_ovl"] / cds_len >= MIN_OVERLAP:
                            tier = "T3_overlap50"
                        else:
                            tier = "T4_overlap_partial"
                else:
                    tiers["T5_none"] = tiers.get("T5_none", 0) + 1
                    unmatched.append((hid, target, acc, h.contig, h.start, h.stop, h.symbol,
                                      h.element_type, "no_overlapping_cds"))
                    synthetic.append((acc, target, h.contig, hs, he, h.strand, h.symbol,
                                      h.element_name, h.element_type))
                    continue

            tiers[tier] = tiers.get(tier, 0) + 1
            hit_len = he - hs + 1
            cds_len = best["_end"] - best["_start"] + 1
            rows.append((hid, target, acc, h.contig, hs, he, h.strand, h.element_type, h.subtype,
                         h.symbol, h.gene_family, h.amr_class, h.subclass, h.method,
                         h.hit_quality, h.protein, best["gene_id"], best["_ord"],
                         best["contig_n_genes"], tier, best["_start"] - hs, best["_end"] - he,
                         best["_ovl"], round(best["_ovl"] / hit_len, 4),
                         round(best["_ovl"] / cds_len, 4), len(cands) or 1,
                         (1 if best["protein_id"] == h.protein else 0) if h.protein else "",
                         best["protein_md5"], best["is_pseudo"],
                         1 if tier in ("T1_exact", "T2_stop_codon") else 0))
        if n % 1000 == 0:
            log(f"  {n:,} genomes ({time.time()-t0:.0f}s)")

    ISLANDS.mkdir(parents=True, exist_ok=True)
    link = pd.DataFrame(rows, columns=OUT_COLS)
    write_tsv(link, "hit_gene_link.tsv", folder=ISLANDS)
    write_tsv(pd.DataFrame(unmatched, columns=["hit_id", "target_acc", "download_acc", "contig",
                                               "start", "stop", "symbol", "element_type", "reason"]),
              "hit_gene_unmatched.tsv", folder=ISLANDS)
    write_tsv(pd.DataFrame(synthetic, columns=["download_acc", "target_acc", "contig", "start",
                                               "stop", "strand", "symbol", "element_name",
                                               "element_type"]),
              "synthetic_features.tsv", folder=PANGENOME)

    # ---- QC ------------------------------------------------------------------------------------
    print("\n===== QC 12_link_hits_to_genes =====")
    total = len(hits)
    for tier in ("T1_exact", "T2_stop_codon", "T3_overlap50", "T4_overlap_partial", "T5_none"):
        n = tiers.get(tier, 0)
        print(f"  {tier:20s} {n:8,}  {100*n/total:6.2f}%")
    # NCBI ships a handful of assemblies with no genomic.gff at all; their hits cannot be placed
    # on a gene and are listed in hit_gene_unmatched.tsv. Anything beyond those IS a bug.
    absent_accs = sorted({r[1] for r in unmatched if r[8] == "genome_absent"})
    print(f"hits whose genome has no annotation: {n_genome_absent} "
          f"in {len(absent_accs)} genome(s) {', '.join(absent_accs[:4])}")

    if len(link):
        has_prot = link[link["hit_protein"] != ""]
        t1 = has_prot[has_prot["match_tier"] == "T1_exact"]
        print(f"hits with a protein accession: {len(has_prot):,}; exact-matched "
              f"{100*len(t1)/max(len(has_prot),1):.2f}%  [expect ~99.96%]")
        agree = pd.to_numeric(t1["protein_agrees"], errors="coerce")
        print(f"T1 protein accession agreement: {100*agree.mean():.2f}%  [expect ~99.95%]")
        weak = link[link["match_tier"].isin(["T3_overlap50", "T4_overlap_partial"])]
        wagree = pd.to_numeric(weak["protein_agrees"], errors="coerce").dropna()
        if len(wagree):
            print(f"T3/T4 protein agreement: {100*wagree.mean():.2f}%  "
                  f"[low is expected - this is why propagates_label exists]")
        noprot = link[link["hit_protein"] == ""]
        print(f"hits without a protein accession: {len(noprot):,}; of those "
              f"{100*(noprot['match_tier'] != 'T1_exact').mean():.1f}% needed a non-exact tier")
        print(f"matched CDS that are pseudogenes: "
              f"{(pd.to_numeric(link['gene_is_pseudo'])==1).sum():,}")
        print(f"propagates_label=1: {int(link['propagates_label'].sum()):,} "
              f"({100*link['propagates_label'].mean():.2f}%)")
        for et in ("AMR", "STRESS", "VIRULENCE"):
            sub = link[link["element_type"] == et]
            print(f"  {et:10s} linked {len(sub):,}")

    # independent cross-check against the md5 computed by a completely different code path (08)
    seq_index = ISLANDS.parent / "AMR_sequences" / "sequence_index.tsv"
    if seq_index.exists() and len(link):
        si = read_tsv(seq_index)
        si["hit_id"] = (si["target_acc"] + ":" + si["contig"] + ":" + si["start"] + "-" + si["stop"]
                        + ":" + si["strand"] + ":" + si["symbol"])
        m = link.merge(si[["hit_id", "protein_md5"]], on="hit_id", how="inner")
        m = m[(m["protein_md5"] != "") & (m["gene_protein_md5"] != "")]
        if len(m):
            same = (m["protein_md5"] == m["gene_protein_md5"]).mean()
            print(f"\nmd5 cross-check vs AMR_sequences/sequence_index.tsv: "
                  f"{100*same:.3f}% agree on {len(m):,} hits  [expect >= 99.9%]")
    log("done")


if __name__ == "__main__":
    main()
