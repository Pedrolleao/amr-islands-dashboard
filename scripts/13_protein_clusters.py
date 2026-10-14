"""13 - cluster the non-redundant proteome with MMseqs2 (90% identity, 80% coverage).

Input is `proteins_nr.faa` from 11: one record per DISTINCT protein sequence, named by its md5.
For Citrobacter that is ~1.5M sequences standing in for ~30.8M CDS copies. Clustering the deduplicated
set is not an approximation - identical sequences must land in the same cluster by definition, and 15
expands membership back to every copy.

PARAMETER CHOICES, AND WHY

  easy-linclust, not easy-cluster
      At 90% identity two sequences share very many exact k-mers, which is the regime linclust was
      built for; the cascaded all-vs-all search buys sensitivity that only matters at low identity.
      The precedent project's 4-day run was 38M sequences at 30% identity with -s 7.5 - three
      compounding cost multipliers, none of which apply here.

  --min-seq-id 0.90 -c 0.80 --cov-mode 0
      `--cov-mode 0` requires 80% coverage of BOTH sequences. This matters more than usual here,
      because cluster membership PROPAGATES AN AMR LABEL (see 16): a one-sided coverage mode would
      let a 150 aa fragment join a 600 aa family, chaining pseudogenised transposase fragments into
      the intact family and letting a truncated resistance gene drag unrelated short proteins into
      an "AMR" cluster. Bidirectional coverage enforces length homogeneity.

  --cluster-reassign belongs to easy-cluster ONLY
      MMseqs' own help explains why: "cascaded clustering can cluster sequences that do not fulfil
      the clustering criteria". That is a property of the multi-round cascade, not of linclust,
      which aligns each member against its representative directly. So linclust needs no
      reassignment pass and the flag does not exist for it (mmseqs 17 rejects it outright).
      `--method cluster` runs the cascaded algorithm instead, with reassignment enabled, so the two
      can be compared on the same input - see the QC block.

  tmp on /tmp, never on /mnt/e
      MMseqs' temporary directory is thousands of files with heavy random I/O; on the NTFS mount
      that is an order of magnitude slower. The script refuses to run if /tmp is short on space.

Calibration: --subsample 0.1 runs the identical command on a 10% sample first. Linclust is roughly
linear, so that 2-minute run predicts the real one better than any estimate.

Usage: python3 scripts/13_protein_clusters.py [--group G] [--subsample F] [--method linclust|cluster]
                                              [--tmp DIR] [--threads N]
"""
import os
import random
import shutil
import subprocess
import sys
import time
from pathlib import Path

from amr_common import PANGENOME

MIN_SEQ_ID = "0.90"
COVERAGE = "0.80"
COV_MODE = "0"            # bidirectional: 80% of query AND target
THREADS = "32"
SPLIT_MEMORY = "30G"
KMER_PER_SEQ = "80"       # above the default 21; nearly free at high identity, closes the gap
MIN_FREE_GB = 60          # refuse to start if the tmp filesystem has less than this


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def free_gb(path):
    st = os.statvfs(path)
    return st.f_bavail * st.f_frsize / 1e9


def count_fasta(path):
    n = 0
    with open(path) as fh:
        for line in fh:
            n += line.startswith(">")
    return n


def subsample(src, dst, fraction, seed=20260921):
    rng = random.Random(seed)
    kept = 0
    with open(src) as fh, open(dst, "w") as out:
        keep = False
        for line in fh:
            if line.startswith(">"):
                keep = rng.random() < fraction
                kept += keep
            if keep:
                out.write(line)
    return kept


def main():
    argv = sys.argv[1:]
    threads = argv[argv.index("--threads") + 1] if "--threads" in argv else THREADS
    tmp_root = Path(argv[argv.index("--tmp") + 1] if "--tmp" in argv else "/tmp/mmseqs_amr")
    frac = float(argv[argv.index("--subsample") + 1]) if "--subsample" in argv else None
    method = argv[argv.index("--method") + 1] if "--method" in argv else "linclust"
    if method not in ("linclust", "cluster"):
        raise SystemExit("--method must be linclust or cluster")

    mmseqs = shutil.which("mmseqs")
    if not mmseqs:
        raise SystemExit("mmseqs not found on PATH")

    src = PANGENOME / "proteins_nr.faa"
    if not src.exists():
        raise SystemExit(f"{src} missing - run 11_gene_index.py first")

    tmp_root.mkdir(parents=True, exist_ok=True)
    free = free_gb(tmp_root)
    log(f"mmseqs {subprocess.run([mmseqs, 'version'], capture_output=True, text=True).stdout.strip()}"
        f"   tmp={tmp_root} ({free:.0f} GB free)")
    if free < MIN_FREE_GB:
        raise SystemExit(f"only {free:.0f} GB free on {tmp_root}; need >= {MIN_FREE_GB} GB. "
                         f"Pass --tmp DIR to use another filesystem (but not one on /mnt/e).")

    if frac:
        sample = tmp_root / "sample.faa"
        n = subsample(src, sample, frac)
        log(f"calibration run on {n:,} sequences ({frac:.0%} sample)")
        src, out_prefix = sample, tmp_root / f"sample_clu_{method}"
    else:
        out_prefix = tmp_root / f"nr90_{method}"

    n_in = count_fasta(src)
    work = tmp_root / "work"
    if work.exists():
        shutil.rmtree(work)

    cmd = [mmseqs, f"easy-{method}", str(src), str(out_prefix), str(work),
           "--min-seq-id", MIN_SEQ_ID, "-c", COVERAGE, "--cov-mode", COV_MODE,
           "--cluster-mode", "0", "--threads", threads,
           "--split-memory-limit", SPLIT_MEMORY, "--remove-tmp-files", "1", "-v", "3"]
    if method == "linclust":
        cmd += ["--kmer-per-seq", KMER_PER_SEQ]
    else:
        cmd += ["--cluster-reassign", "1", "-s", "4.0"]
    log("$ " + " ".join(cmd))
    t0 = time.time()
    log_path = PANGENOME / (f"mmseqs_sample_{method}.log" if frac else f"mmseqs_{method}.log")
    with open(log_path, "w") as lf:
        lf.write(" ".join(cmd) + "\n\n")
        lf.flush()
        rc = subprocess.run(cmd, stdout=lf, stderr=subprocess.STDOUT).returncode
    elapsed = time.time() - t0
    if rc:
        raise SystemExit(f"mmseqs failed (rc={rc}); see {log_path}")

    tsv = Path(f"{out_prefix}_cluster.tsv")
    reps = {}
    n_members = 0
    with open(tsv) as fh:
        for line in fh:
            rep, member = line.rstrip("\n").split("\t")
            reps[rep] = reps.get(rep, 0) + 1
            n_members += 1
    sizes = sorted(reps.values(), reverse=True)
    singletons = sum(1 for s in sizes if s == 1)

    print(f"\n===== QC 13_protein_clusters [{method}]"
          f"{' (calibration)' if frac else ''} =====")
    print(f"input sequences: {n_in:,}")
    print(f"members in the output: {n_members:,}   [must equal the input: "
          f"{n_members == n_in}]")
    print(f"clusters: {len(reps):,}   singletons: {singletons:,} "
          f"({100*singletons/max(len(reps),1):.1f}% of clusters)")
    print(f"largest clusters: {', '.join(f'{s:,}' for s in sizes[:5])}")
    print(f"runtime: {elapsed/60:.1f} min on {threads} threads")

    if frac:
        print(f"\nextrapolated to the full set: ~{elapsed*(1/frac)/60:.0f} min "
              f"(linclust is ~linear; add headroom)")
        print("calibration only - no output copied. Re-run without --subsample.")
        return

    shutil.copy(tsv, PANGENOME / "clusters_raw.tsv")
    shutil.copy(f"{out_prefix}_rep_seq.fasta", PANGENOME / "cluster_reps.faa")
    print(f"\nwrote {(PANGENOME/'clusters_raw.tsv').relative_to(PANGENOME.parent.parent)} "
          f"and cluster_reps.faa")
    log("done")


if __name__ == "__main__":
    main()
