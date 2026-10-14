"""28 - what each predicted structure looks like, by structure comparison against the PDB.

Sequence search finds a homolog only when the sequence still resembles one. Structure survives far
longer than sequence does, so Foldseek finds relatives that Pfam and the NCBI annotation miss - which
is exactly the situation of the clusters that came out of step 24 with no name and low confidence.

It only applies to clusters that HAVE a structure, so it complements step 27 rather than replacing
it: Pfam covers all 200,000 clusters, this covers the few hundred that were folded.

Read the two together. A Pfam hit says which family the sequence belongs to; a Foldseek hit says
which solved structure the fold resembles, with a TM-score for how well. TM-score above 0.5 is the
usual line for "the same fold"; below 0.3 is noise. A low-pLDDT model gives unreliable hits, so the
model's own confidence is carried into the output and must be read alongside the score.

Usage: tools/viz_env/bin/python is NOT needed; plain python3.
  python3 scripts/28_foldseek_hits.py [--group G] [--set amr|stress|virulence] [--threads N]
                                      [--evalue E] [--max-hits N]
  --db pdb,afdb_sp   which structure databases to search (default: both)
  Output in Domains[_<set>]/:
    foldseek_hits.tsv   every hit kept, best first per cluster, with the database it came from
    foldseek_best.tsv   one row per cluster: its best hit, or an explicit "no hit" row
"""
import subprocess
import sys
import time

import pandas as pd

from amr_common import (DOMAINS, FOLDSEEK, FOLDSEEK_DBS, SET, STRUCTURES, read_tsv, set_dir,
                        write_tsv)

DOMAINS = set_dir(DOMAINS)
STRUCTURES = set_dir(STRUCTURES)

COLS = ["query", "target", "theader", "fident", "alnlen", "qlen", "tlen", "qstart", "qend",
        "tstart", "tend", "evalue", "bits", "alntmscore", "qtmscore", "ttmscore", "lddt", "qcov",
        "tcov", "prob"]
SAME_FOLD = 0.5          # TM-score at or above this is normally read as the same fold


def main():
    t0 = time.time()
    threads = int(sys.argv[sys.argv.index("--threads") + 1]) if "--threads" in sys.argv else 32
    evalue = sys.argv[sys.argv.index("--evalue") + 1] if "--evalue" in sys.argv else "0.001"
    max_hits = sys.argv[sys.argv.index("--max-hits") + 1] if "--max-hits" in sys.argv else "10"
    if SET == "all":
        raise SystemExit("--set all has no structures; steps 24-25 run on amr, stress, virulence")
    which = (sys.argv[sys.argv.index("--db") + 1].split(",") if "--db" in sys.argv
             else list(FOLDSEEK_DBS))
    unknown = [w for w in which if w not in FOLDSEEK_DBS]
    if unknown:
        raise SystemExit(f"unknown --db {unknown}; known: {', '.join(FOLDSEEK_DBS)}")
    if not FOLDSEEK.exists():
        raise SystemExit(f"the foldseek binary is not at {FOLDSEEK}")
    for w in which:
        if not FOLDSEEK_DBS[w].with_suffix(".index").exists():
            raise SystemExit(f"the {w} database is not at {FOLDSEEK_DBS[w]}")
    pdb_dir = STRUCTURES / "pdb"
    structures = sorted(pdb_dir.glob("*.pdb"))
    if not structures:
        raise SystemExit(f"no structures in {pdb_dir} - run step 24 first")
    DOMAINS.mkdir(parents=True, exist_ok=True)

    tmp = DOMAINS / "_foldseek_tmp"
    print(f"set {SET!r}: {len(structures)} structures vs {', '.join(which)}, {threads} threads",
          flush=True)
    parts = []
    for w in which:
        raw = DOMAINS / f"_foldseek_raw_{w}.tsv"
        cmd = [str(FOLDSEEK), "easy-search", str(pdb_dir), str(FOLDSEEK_DBS[w]), str(raw), str(tmp),
               "--format-output", ",".join(COLS), "--format-mode", "4",
               "-e", evalue, "--max-seqs", "1000", "--threads", str(threads)]
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL)
        part = pd.read_csv(raw, sep="\t")
        part["database"] = w
        parts.append(part)
        print(f"  {w}: {len(part):,} hits ({time.time() - t0:.0f} s)", flush=True)
    h = pd.concat(parts, ignore_index=True)
    h["name"] = h["query"].str.replace(r"\.pdb.*$", "", regex=True)

    # Map the file back to its row by NAME, not by splitting the file name: a protein folded in
    # modules by step 29 has one file per module and a cluster_id like CFR-014064.m2, which a split
    # on "__" would collapse onto the parent and silently lose every module's hits.
    summ = read_tsv(STRUCTURES / "structure_summary.tsv")
    by_name = dict(zip(summ.name, summ.cluster_id))
    h["cluster_id"] = h.name.map(by_name)
    lost = h.cluster_id.isna()
    if lost.any():
        h.loc[lost, "cluster_id"] = h.loc[lost, "name"].str.split("__").str[0]
        print(f"  {h.loc[lost, 'name'].nunique()} structure(s) not in the summary, "
              f"id taken from the file name", flush=True)
    h = h.sort_values(["cluster_id", "alntmscore"], ascending=[True, False])
    h = h.groupby("cluster_id", sort=False).head(int(max_hits))

    # the model's own confidence belongs next to the hit: a bad model gives a meaningless match
    conf = dict(zip(summ.cluster_id, pd.to_numeric(summ.mean_plddt, errors="coerce")))
    fam = dict(zip(summ.cluster_id, summ.family.fillna("")))
    prod = dict(zip(summ.cluster_id, summ.product_top.fillna("")))
    h["mean_plddt"] = h.cluster_id.map(conf)
    h["gene"] = h.cluster_id.map(fam)
    h["same_fold"] = h.alntmscore >= SAME_FOLD
    write_tsv(h[["cluster_id", "name", "gene", "mean_plddt", "database", "target", "theader",
                 "alntmscore", "qtmscore", "ttmscore", "lddt", "fident", "evalue", "bits", "qcov",
                 "tcov", "prob", "same_fold"]], "foldseek_hits.tsv", DOMAINS)

    # one row per structure, including the ones that matched nothing
    best = h.groupby("cluster_id", sort=False).head(1).set_index("cluster_id")
    rows = []
    for c in summ.itertuples():
        # the AMR summary predates the status column: everything in it was predicted
        if str(getattr(c, "status", "predicted")) != "predicted":
            continue
        b = best.loc[c.cluster_id] if c.cluster_id in best.index else None
        rows.append({
            "cluster_id": c.cluster_id, "name": c.name, "gene": fam.get(c.cluster_id, ""),
            "product_top": prod.get(c.cluster_id, ""), "mean_plddt": conf.get(c.cluster_id),
            "best_database": "" if b is None else b.database,
            "best_target": "" if b is None else b.target,
            "best_target_desc": "" if b is None else b.theader,
            "tm_score": "" if b is None else round(float(b.alntmscore), 3),
            "seq_identity": "" if b is None else round(float(b.fident), 3),
            "evalue": "" if b is None else b.evalue,
            "verdict": ("no structural hit" if b is None else
                        "same fold" if b.alntmscore >= SAME_FOLD else "weak, below TM 0.5"),
        })
    best_df = pd.DataFrame(rows)
    write_tsv(best_df, "foldseek_best.tsv", DOMAINS)

    print("\nchecks")
    n = len(best_df)
    print(f"  structures searched: {n}; with a hit: {(best_df.best_target != '').sum()}; "
          f"same fold (TM >= {SAME_FOLD}): {(best_df.verdict == 'same fold').sum()}")
    lowconf = best_df[pd.to_numeric(best_df.mean_plddt, errors='coerce') < 70]
    print(f"  of the {len(lowconf)} models below pLDDT 70, "
          f"{(lowconf.verdict == 'same fold').sum()} still match a known fold")
    print(f"\ndone in {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
