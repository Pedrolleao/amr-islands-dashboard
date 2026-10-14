"""14 - mint stable protein-cluster ids (CFR-000001) from the raw MMseqs2 output.

MMseqs names a cluster after whichever member it picked as representative, which is an md5 here and
can change between runs. Every table, figure and future browser payload refers to clusters by id, so
the id must be stable and must not depend on run order.

TWO LESSONS TAKEN FROM THE PRECEDENT PROJECT, BOTH MEASURED THERE

1. THE SORT KEY MUST BE TOTAL. Its v1 sorted clusters by size alone; ties fell to dict insertion
   order, and a re-run silently moved 3,094 of 3,128 ids onto different clusters. Here the key is
   (-n_members, -n_seqs, rep_md5) and rep_md5 is unique per cluster, so the order is total.

2. AN ID IS KEYED ON CONTENT, NOT ON POSITION. `content_md5` is the md5 of the cluster's sorted
   member list. A cluster whose membership is unchanged keeps its id across a re-clustering;
   existing ids are never renumbered and retired ids are never reused. `cluster_registry_changes.tsv`
   records what was reused, minted and retired so a re-run is auditable rather than mysterious.

Note the name: these are PROTEIN clusters. `AMR_annotation/cluster_amr_summary.tsv` in this same
project means Pathogen Detection SNP clusters (PDS_acc) - a different thing entirely.

Usage: python3 scripts/14_cluster_registry.py [--group G]
"""
import hashlib
import time

import pandas as pd

from amr_common import CLUSTER_PREFIX, PANGENOME, read_tsv, write_tsv

REGISTRY = "cluster_registry.tsv"


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def main():
    raw = PANGENOME / "clusters_raw.tsv"
    if not raw.exists():
        raise SystemExit(f"{raw} missing - run 13_protein_clusters.py first")

    nr = read_tsv(PANGENOME / "protein_nr.tsv")
    copies_of = dict(zip(nr["protein_md5"], nr["n_copies"].astype(int)))
    genomes_of = dict(zip(nr["protein_md5"], nr["n_genomes"].astype(int)))

    members = {}
    with open(raw) as fh:
        for line in fh:
            rep, member = line.rstrip("\n").split("\t")
            members.setdefault(rep, []).append(member)
    log(f"{len(members):,} clusters, {sum(len(v) for v in members.values()):,} sequences")

    rows = []
    for rep, mem in members.items():
        mem_sorted = sorted(mem)
        rows.append({
            "rep_md5": rep,
            "content_md5": hashlib.md5("\n".join(mem_sorted).encode()).hexdigest(),
            "n_seqs": len(mem),
            "n_members": sum(copies_of.get(m, 0) for m in mem),
            # genome counts per sequence cannot simply be summed (the same genome may hold two
            # different sequences of the cluster), so this is an upper bound; 15 computes the exact
            # value from the gene index.
            "n_genomes_upper": sum(genomes_of.get(m, 0) for m in mem),
        })
    # total order: size, then distinct sequences, then the unique representative md5
    rows.sort(key=lambda r: (-r["n_members"], -r["n_seqs"], r["rep_md5"]))

    prev_path = PANGENOME / REGISTRY
    previous = read_tsv(prev_path) if prev_path.exists() else pd.DataFrame(
        columns=["cluster_id", "content_md5"])
    known = dict(zip(previous.get("content_md5", []), previous.get("cluster_id", [])))
    used_numbers = {int(cid.rsplit("-", 1)[1]) for cid in known.values()} if len(known) else set()
    next_number = max(used_numbers) + 1 if used_numbers else 1
    today = time.strftime("%Y-%m-%d")
    first_seen_of = dict(zip(previous.get("content_md5", []), previous.get("first_seen", [])))

    n_reused = n_minted = 0
    for r in rows:
        cid = known.get(r["content_md5"])
        if cid:
            n_reused += 1
            r["first_seen"] = first_seen_of.get(r["content_md5"], today)
        else:
            cid = f"{CLUSTER_PREFIX}-{next_number:06d}"
            next_number += 1
            n_minted += 1
            r["first_seen"] = today
        r["cluster_id"] = cid
        r["last_seen"] = today

    reg = pd.DataFrame(rows)[["cluster_id", "content_md5", "rep_md5", "n_seqs", "n_members",
                              "n_genomes_upper", "first_seen", "last_seen"]]
    retired = sorted(set(known) - set(reg["content_md5"]))
    write_tsv(reg, REGISTRY, folder=PANGENOME)
    if len(known):
        write_tsv(pd.DataFrame({"retired_content_md5": retired,
                                "retired_cluster_id": [known[c] for c in retired]}),
                  "cluster_registry_changes.tsv", folder=PANGENOME)

    print("\n===== QC 14_cluster_registry =====")
    print(f"clusters: {len(reg):,}   reused ids: {n_reused:,}   newly minted: {n_minted:,}   "
          f"retired: {len(retired):,}")
    print(f"cluster_id unique: {reg['cluster_id'].is_unique}   "
          f"content_md5 unique: {reg['content_md5'].is_unique}")
    print(f"sum(n_seqs) == input sequences: "
          f"{reg['n_seqs'].sum() == sum(len(v) for v in members.values())}")
    print(f"largest cluster: {reg['n_members'].max():,} CDS copies "
          f"({reg['n_seqs'].max():,} distinct sequences)")
    print(f"singleton clusters: {(reg['n_seqs'] == 1).sum():,}")
    log("done")


if __name__ == "__main__":
    main()
