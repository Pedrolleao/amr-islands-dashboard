"""15 - expand clusters back to every CDS, and describe what each cluster is.

13 clustered ~1.5M distinct sequences; this maps them back onto all ~30.8M CDS copies and builds the
per-cluster table that 16 turns into island anchors.

HOW THE COUNTS STAY EXACT WITHOUT A 30M-ROW DATAFRAME. `gene_index.tsv.gz` is written as one gzip
member per genome, in genome order, so this script walks it genome by genome. Distinct-genome counts
per cluster then need no global sets at all: within one genome, collect the clusters seen, then add
one to each of their genome counters. Peak memory is one genome's worth of rows plus the counter
dicts.

WHAT THE COLUMNS ARE FOR

  prevalence / pangenome_class
      A cluster in nearly every genome is a core chromosomal gene. This is the number that separates
      a mobile resistance cassette from the backbone, and 16 tiers anchors on it.

  n_amr_members vs n_point_members
      Kept apart on purpose. A point mutation (gyrA_T83I) means an essential gene carries a
      resistance substitution; it is not an acquired resistance element, and mixing the two would
      make every genome's gyrA look like cargo.

  evidence_tier / propagation_ratio
      Anchoring is by cluster, so ONE hit can label thousands of members. `propagation_ratio` is
      members per annotated hit - about 1.0 for a genuinely mobile gene, large for a backbone family
      that acquired its label from a handful of hits. `evidence_tier = PROPAGATED_ONLY` marks a
      cluster whose AMR status rests entirely on weak (T3/T4 overlap) matches, which 12 already
      flagged with `propagates_label = 0`. Both are reported, never filtered.

  product_top_frac / aa_len_cv
      Cluster hygiene. A clean cluster agrees on its product name and its length; a low agreement or
      a wide length spread means the cluster merged things it should not have.

Usage: python3 scripts/15_cluster_annotation.py [--group G]
"""
import gzip
import time
from collections import Counter, defaultdict

import pandas as pd

from amr_common import ISLANDS, PANGENOME, read_tsv, write_tsv

CORE = 0.99           # Roary-style pangenome cut-offs, so the words mean what readers expect
SOFT_CORE = 0.95
SHELL = 0.15
MAX_PRODUCTS = 40     # per cluster, keep at most this many distinct product names while counting


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def join_top(counter, n=12):
    return ";".join(s for s, _ in counter.most_common(n) if s)


def main():
    t0 = time.time()
    reg = read_tsv(PANGENOME / "cluster_registry.tsv")
    log(f"{len(reg):,} clusters")

    # sequence md5 -> cluster id
    cid_of_md5 = {}
    with open(PANGENOME / "clusters_raw.tsv") as fh:
        rep_to_cid = dict(zip(reg["rep_md5"], reg["cluster_id"]))
        for line in fh:
            rep, member = line.rstrip("\n").split("\t")
            cid_of_md5[member] = rep_to_cid[rep]
    log(f"{len(cid_of_md5):,} sequences mapped to a cluster")

    n_members = Counter()          # CDS copies
    n_genomes = Counter()          # distinct genomes
    max_in_genome = Counter()
    n_single_copy = Counter()      # genomes where the cluster has exactly one copy
    products = defaultdict(Counter)
    symbols = defaultdict(Counter)
    aa_sum, aa_min, aa_max, aa_sqsum = Counter(), {}, {}, Counter()
    near_end = Counter()
    on_plasmid = Counter()
    n_pseudo_total = 0
    n_rows = n_no_md5 = 0

    header = None
    cur_genome, local = None, Counter()
    local_products = []
    map_path = PANGENOME / "gene_cluster_map.tsv.gz"
    offsets = []

    def close_genome():
        for cid, cnt in local.items():
            n_genomes[cid] += 1
            if cnt > max_in_genome[cid]:
                max_in_genome[cid] = cnt
            if cnt == 1:
                n_single_copy[cid] += 1

    with gzip.open(PANGENOME / "gene_index.tsv.gz", "rt") as fh, open(map_path, "wb") as mout:
        header = fh.readline().rstrip("\n").split("\t")
        col = {c: i for i, c in enumerate(header)}
        buf = []
        for line in fh:
            f = line.rstrip("\n").split("\t")
            n_rows += 1
            acc = f[col["download_acc"]]
            if acc != cur_genome:
                if cur_genome is not None:
                    close_genome()
                    start = mout.tell()
                    comp = gzip.compress(("\n".join(buf) + "\n").encode(), 6)
                    mout.write(comp)
                    offsets.append((cur_genome, start, len(comp), len(buf)))
                cur_genome, local, buf = acc, Counter(), []
            md5 = f[col["protein_md5"]]
            cid = cid_of_md5.get(md5, "") if md5 else ""
            buf.append(f'{f[col["gene_id"]]}\t{cid}')
            if f[col["is_pseudo"]] == "1":
                n_pseudo_total += 1
            if not cid:
                n_no_md5 += 1
                continue
            local[cid] += 1
            n_members[cid] += 1
            prod = f[col["product"]]
            pc = products[cid]
            if prod in pc or len(pc) < MAX_PRODUCTS:
                pc[prod] += 1
            sym = f[col["gene_symbol"]]
            if sym:
                sc = symbols[cid]
                if sym in sc or len(sc) < MAX_PRODUCTS:
                    sc[sym] += 1
            aa = int(f[col["aa_len"]] or 0)
            aa_sum[cid] += aa
            aa_sqsum[cid] += aa * aa
            if cid not in aa_min or aa < aa_min[cid]:
                aa_min[cid] = aa
            if cid not in aa_max or aa > aa_max[cid]:
                aa_max[cid] = aa
            ordn, ncg = int(f[col["ordinal"]]), int(f[col["contig_n_genes"]] or 0)
            if ordn <= 3 or (ncg and ordn >= ncg - 2):
                near_end[cid] += 1
            if f[col["contig_molecule"]] == "plasmid":
                on_plasmid[cid] += 1
            if n_rows % 5_000_000 == 0:
                log(f"  {n_rows:,} CDS ({time.time()-t0:.0f}s)")
        if cur_genome is not None:
            close_genome()
            start = mout.tell()
            comp = gzip.compress(("\n".join(buf) + "\n").encode(), 6)
            mout.write(comp)
            offsets.append((cur_genome, start, len(comp), len(buf)))

    n_genomes_total = len(offsets)
    log(f"{n_rows:,} CDS over {n_genomes_total:,} genomes")
    write_tsv(pd.DataFrame(offsets, columns=["download_acc", "byte_offset", "byte_length", "n_rows"]),
              "gene_cluster_map_offsets.tsv", folder=PANGENOME)

    # ---- AMR / stress / virulence evidence per cluster ------------------------------------------
    link = read_tsv(ISLANDS / "hit_gene_link.tsv")
    link["cluster_id"] = link["gene_protein_md5"].map(cid_of_md5).fillna("")
    linked = link[link["cluster_id"] != ""].copy()
    linked["prop"] = pd.to_numeric(linked["propagates_label"])
    ev = {}
    for cid, grp in linked.groupby("cluster_id", sort=False):
        amr = grp[(grp["element_type"] == "AMR") & (grp["subtype"] != "POINT")]
        point = grp[(grp["element_type"] == "AMR") & (grp["subtype"] == "POINT")]
        strong = amr[amr["prop"] == 1]
        ev[cid] = {
            "n_amr_hits": len(amr), "n_amr_hits_strong": len(strong),
            "n_point_hits": len(point),
            "n_stress_hits": int((grp["element_type"] == "STRESS").sum()),
            "n_virulence_hits": int((grp["element_type"] == "VIRULENCE").sum()),
            "n_genomes_with_amr_hit": amr["target_acc"].nunique(),
            "amr_symbols": ";".join(sorted(set(amr["symbol"]))[:20]),
            "amr_gene_families": ";".join(sorted(set(amr["gene_family"]) - {""})),
            "amr_classes": ";".join(sorted({c for v in amr["class"] for c in v.split("/") if c})),
            "amr_subclasses": ";".join(sorted({c for v in amr["subclass"]
                                               for c in v.split("/") if c})),
            "point_symbols": ";".join(sorted(set(point["symbol"]))),
            "stress_symbols": ";".join(sorted(set(grp.loc[grp["element_type"] == "STRESS",
                                                          "symbol"]))[:20]),
            "virulence_symbols": ";".join(sorted(set(grp.loc[grp["element_type"] == "VIRULENCE",
                                                             "symbol"]))[:20]),
        }
    log(f"{len(ev):,} clusters carry at least one hit")

    # ---- assemble -------------------------------------------------------------------------------
    rows = []
    for r in reg.itertuples():
        cid = r.cluster_id
        nm = n_members[cid]
        ng = n_genomes[cid]
        prev = ng / n_genomes_total if n_genomes_total else 0
        pc = products[cid]
        sc = symbols[cid]
        top_p, top_pn = (pc.most_common(1)[0] if pc else ("", 0))
        top_s, top_sn = (sc.most_common(1)[0] if sc else ("", 0))
        mean = aa_sum[cid] / nm if nm else 0
        var = (aa_sqsum[cid] / nm - mean * mean) if nm else 0
        e = ev.get(cid, {})
        n_amr = e.get("n_amr_hits", 0)
        n_amr_strong = e.get("n_amr_hits_strong", 0)
        n_point = e.get("n_point_hits", 0)
        etypes = [t for t, k in (("AMR", n_amr or n_point), ("STRESS", e.get("n_stress_hits", 0)),
                                 ("VIRULENCE", e.get("n_virulence_hits", 0))) if k]
        rows.append({
            "cluster_id": cid, "rep_md5": r.rep_md5, "n_seqs": r.n_seqs, "n_members": nm,
            "n_genomes": ng, "prevalence": round(prev, 6),
            "pangenome_class": ("core" if prev >= CORE else "soft_core" if prev >= SOFT_CORE
                                else "shell" if prev >= SHELL else "cloud"),
            "mean_copies_per_genome": round(nm / ng, 3) if ng else 0,
            "max_copies_in_one_genome": max_in_genome[cid],
            "single_copy_fraction": round(n_single_copy[cid] / ng, 3) if ng else 0,
            "product_top": top_p, "product_top_frac": round(top_pn / nm, 3) if nm else 0,
            "n_distinct_products": len(pc), "products_seen": join_top(pc, 5),
            "gene_symbol_top": top_s,
            "gene_symbol_top_frac": round(top_sn / nm, 3) if nm else 0,
            "aa_len_min": aa_min.get(cid, 0), "aa_len_mean": round(mean, 1),
            "aa_len_max": aa_max.get(cid, 0),
            "aa_len_cv": round((var ** 0.5) / mean, 3) if mean else 0,
            "frac_near_contig_end": round(near_end[cid] / nm, 3) if nm else 0,
            "frac_on_plasmid_contig": round(on_plasmid[cid] / nm, 3) if nm else 0,
            "element_type": ";".join(etypes),
            "n_amr_hits": n_amr, "n_amr_hits_strong": n_amr_strong, "n_point_hits": n_point,
            "n_stress_hits": e.get("n_stress_hits", 0),
            "n_virulence_hits": e.get("n_virulence_hits", 0),
            "n_genomes_with_amr_hit": e.get("n_genomes_with_amr_hit", 0),
            "amr_member_fraction": round(n_amr / nm, 4) if nm else 0,
            "propagation_ratio": round(nm / n_amr, 1) if n_amr else "",
            "evidence_tier": ("NONE" if not n_amr else
                              "DIRECT" if n_amr_strong else "PROPAGATED_ONLY"),
            "amr_symbols": e.get("amr_symbols", ""),
            "amr_gene_families": e.get("amr_gene_families", ""),
            "amr_classes": e.get("amr_classes", ""), "amr_subclasses": e.get("amr_subclasses", ""),
            "point_symbols": e.get("point_symbols", ""),
            "stress_symbols": e.get("stress_symbols", ""),
            "virulence_symbols": e.get("virulence_symbols", ""),
        })
    ann = pd.DataFrame(rows)
    write_tsv(ann, "cluster_annotation.tsv", folder=PANGENOME)

    # ---- QC ---------------------------------------------------------------------------------------
    print("\n===== QC 15_cluster_annotation =====")
    print(f"CDS rows read: {n_rows:,}   without a protein (pseudogenes etc): {n_no_md5:,} "
          f"({100*n_no_md5/max(n_rows,1):.2f}%)")
    print(f"sum(n_members) + unclustered == CDS: "
          f"{ann['n_members'].sum() + n_no_md5 == n_rows}")
    print(f"genomes: {n_genomes_total:,}")
    print("\npangenome composition:")
    for k, v in ann["pangenome_class"].value_counts().items():
        print(f"  {k:10s} {v:8,} clusters, {ann.loc[ann.pangenome_class==k,'n_members'].sum():12,} CDS")
    amr = ann[ann["n_amr_hits"] > 0]
    print(f"\nclusters with acquired-AMR evidence: {len(amr):,} "
          f"({amr['n_members'].sum():,} CDS would anchor)")
    print(f"  of those core/soft_core: {(amr['pangenome_class'].isin(['core','soft_core'])).sum():,}"
          f"  -> these are the backbone anchors 16 tiers apart")
    print(f"  evidence_tier PROPAGATED_ONLY: {(amr['evidence_tier']=='PROPAGATED_ONLY').sum():,}")
    print(f"clusters with point-mutation evidence only: "
          f"{((ann['n_point_hits']>0) & (ann['n_amr_hits']==0)).sum():,}")
    print(f"clusters with stress evidence: {(ann['n_stress_hits']>0).sum():,}; "
          f"virulence: {(ann['n_virulence_hits']>0).sum():,}")
    print(f"\nproduct_top_frac < 0.5 (possibly chimeric): "
          f"{(ann['product_top_frac']<0.5).sum():,} clusters")
    print(f"aa_len_cv > 0.2 (wide length spread despite --cov-mode 0): "
          f"{(ann['aa_len_cv']>0.2).sum():,} clusters")
    print("\ntop 10 clusters by CDS count:")
    for r in ann.nlargest(10, "n_members").itertuples():
        print(f"  {r.cluster_id} {r.n_members:8,} CDS  prev={r.prevalence:.3f} "
              f"{r.pangenome_class:9s} {r.product_top[:50]}")
    log(f"done in {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
