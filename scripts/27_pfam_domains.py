"""27 - Pfam domain architecture of every cluster representative.

Answers what each protein is BUILT FROM, which the NCBI product name and the AMRFinderPlus symbol do
not: the ordered list of Pfam domains and where each one starts and ends. That is what lets a
multi-domain protein be read - and, for clbB, cut into foldable modules.

Pfam 38.0 ships with the InterProScan installed on this machine, already pressed, so nothing is
downloaded and the version is pinned by that path.

Two Pfam conventions are applied, and both matter:
  * the family's own curated **gathering threshold** (--cut_ga) decides what counts as a hit, not a
    generic E-value - a permissive E-value would fill the table with junk;
  * **clan competition** - families of one clan are homologous and their hits overlap, so where two
    hits of the same clan overlap by more than 50% of the shorter one, only the better-scoring
    survives. Without this the architecture double-counts the same region.

hmmsearch, not hmmscan: same models, same answers, an order of magnitude faster, because the
sequences stay resident while the profiles stream past.

Usage: python3 scripts/27_pfam_domains.py [--group G] [--set amr|stress|virulence|all] [--threads N]
  Output in Domains[_<set>]/:
    pfam_domains.tsv      one row per retained domain hit
    cluster_domains.tsv   one row per cluster: the ordered architecture and how much it covers
    _hmmsearch.domtbl     the raw hmmsearch table, kept so the filtering can be re-derived
"""
import subprocess
import sys
import time
from collections import defaultdict

import pandas as pd

from amr_common import DOMAINS, PANGENOME, PFAM_DAT, PFAM_HMM, SET, read_tsv, set_dir, write_tsv
from cluster_sets import load_cluster_set

DOMAINS = set_dir(DOMAINS)
OVERLAP_FRAC = 0.50          # same-clan hits overlapping more than this: keep the better one


def load_pfam_meta():
    """accession -> (name, description, clan, type) from pfam_a.dat."""
    meta, cur = {}, {}
    for line in open(PFAM_DAT):
        if line.startswith("#=GF ID"):
            cur = {"name": line.split(None, 2)[2].strip(), "clan": "", "desc": "", "type": ""}
        elif line.startswith("#=GF AC"):
            cur["acc"] = line.split(None, 2)[2].strip()
        elif line.startswith("#=GF DE"):
            cur["desc"] = line.split(None, 2)[2].strip()
        elif line.startswith("#=GF TP"):
            cur["type"] = line.split(None, 2)[2].strip()
        elif line.startswith("#=GF CL"):
            cur["clan"] = line.split(None, 2)[2].strip()
        elif line.startswith("//") and cur.get("acc"):
            meta[cur["acc"].split(".")[0]] = cur
            cur = {}
    return meta


def write_query_fasta(d, path):
    """The representative of every cluster in the set, named by cluster id."""
    want = dict(zip(d.rep_md5, d.cluster_id))
    lengths, cur = {}, None
    with open(path, "w") as out:
        for line in open(PANGENOME / "proteins_nr.faa"):
            if line.startswith(">"):
                md5 = line[1:].split()[0]
                cur = want.get(md5)
                if cur:
                    out.write(f">{cur}\n")
                    lengths[cur] = 0
            elif cur:
                seq = line.strip().replace("*", "")
                out.write(seq + "\n")
                lengths[cur] += len(seq)
    return lengths


def parse_domtbl(path):
    rows = []
    for line in open(path):
        if line.startswith("#"):
            continue
        f = line.split()
        if len(f) < 23:
            continue
        rows.append({"cluster_id": f[0], "aa_len": int(f[2]), "pfam_name": f[3],
                     "pfam_acc": f[4].split(".")[0], "model_len": int(f[5]),
                     "full_evalue": float(f[6]), "dom_score": float(f[13]),
                     "dom_ievalue": float(f[12]),
                     "hmm_from": int(f[15]), "hmm_to": int(f[16]),
                     "start": int(f[19]), "end": int(f[20])})   # envelope coordinates
    return pd.DataFrame(rows)


def clan_filter(hits, meta):
    """Pfam clan competition: within a clan, overlapping hits compete and the best score wins."""
    keep = []
    for _, g in hits.groupby("cluster_id", sort=False):
        chosen = []
        for h in g.sort_values("dom_score", ascending=False).itertuples():
            clan = meta.get(h.pfam_acc, {}).get("clan", "")
            drop = False
            for c in chosen:
                if clan and clan == meta.get(c.pfam_acc, {}).get("clan", ""):
                    lo, hi = max(h.start, c.start), min(h.end, c.end)
                    ov = max(0, hi - lo + 1)
                    if ov > OVERLAP_FRAC * min(h.end - h.start + 1, c.end - c.start + 1):
                        drop = True
                        break
            if not drop:
                chosen.append(h)
        keep.extend(chosen)
    return pd.DataFrame(keep).drop(columns=["Index"], errors="ignore")


def main():
    t0 = time.time()
    threads = int(sys.argv[sys.argv.index("--threads") + 1]) if "--threads" in sys.argv else 32
    if not PFAM_HMM.exists():
        raise SystemExit(f"Pfam not found at {PFAM_HMM}")
    DOMAINS.mkdir(parents=True, exist_ok=True)

    d = load_cluster_set()
    fa = DOMAINS / "_representatives.faa"
    lengths = write_query_fasta(d, fa)
    n = len(lengths)
    print(f"set {SET!r}: {len(d):,} clusters, {n:,} representative sequences written", flush=True)
    if n != len(d):
        print(f"  WARNING: {len(d) - n} representative(s) not found in proteins_nr.faa", flush=True)

    domtbl = DOMAINS / "_hmmsearch.domtbl"
    cmd = ["hmmsearch", "--cut_ga", "--cpu", str(threads), "--domtblout", str(domtbl),
           "-o", "/dev/null", str(PFAM_HMM), str(fa)]
    print(f"hmmsearch: Pfam 38.0 ({PFAM_HMM.parent.name}) vs {n:,} proteins on {threads} threads",
          flush=True)
    subprocess.run(cmd, check=True)
    print(f"  done ({time.time() - t0:.0f} s)", flush=True)

    meta = load_pfam_meta()
    hits = parse_domtbl(domtbl)
    print(f"raw domain hits above the gathering threshold: {len(hits):,}", flush=True)
    hits = clan_filter(hits, meta)
    print(f"after clan competition: {len(hits):,}", flush=True)

    hits["pfam_desc"] = hits.pfam_acc.map(lambda a: meta.get(a, {}).get("desc", ""))
    hits["clan"] = hits.pfam_acc.map(lambda a: meta.get(a, {}).get("clan", ""))
    hits["pfam_type"] = hits.pfam_acc.map(lambda a: meta.get(a, {}).get("type", ""))
    hits["length"] = hits.end - hits.start + 1
    hits["model_coverage"] = ((hits.hmm_to - hits.hmm_from + 1) / hits.model_len).round(3)
    hits = hits.sort_values(["cluster_id", "start"])
    write_tsv(hits[["cluster_id", "aa_len", "pfam_acc", "pfam_name", "pfam_desc", "pfam_type",
                    "clan", "start", "end", "length", "model_coverage", "dom_score", "dom_ievalue"]],
              "pfam_domains.tsv", DOMAINS)

    # one row per cluster: the architecture, in order along the sequence
    arch = []
    by = {k: g for k, g in hits.groupby("cluster_id")}
    for c in d.itertuples():
        g = by.get(c.cluster_id)
        covered = int(g.length.sum()) if g is not None else 0
        L = lengths.get(c.cluster_id, 0)
        arch.append({
            "cluster_id": c.cluster_id, "gene": c.amr_gene_families, "element": c.amr_classes,
            "product_top": c.product_top, "n_members": c.n_members,
            "aa_len": lengths.get(c.cluster_id, ""),
            "n_domains": 0 if g is None else len(g),
            "n_distinct_families": 0 if g is None else g.pfam_acc.nunique(),
            "architecture": "" if g is None else " - ".join(g.sort_values("start").pfam_name),
            "pfam_accs": "" if g is None else ";".join(g.sort_values("start").pfam_acc),
            "residues_covered": covered,
            "fraction_covered": round(covered / L, 3) if L else "",
        })
    arch = pd.DataFrame(arch)
    write_tsv(arch, "cluster_domains.tsv", DOMAINS)

    print("\nchecks")
    print(f"  clusters with at least one domain: {(arch.n_domains > 0).sum():,} of {len(arch):,} "
          f"({100 * (arch.n_domains > 0).mean():.1f}%)")
    cov = pd.to_numeric(arch.fraction_covered, errors="coerce")
    print(f"  median fraction of the protein covered by domains: {cov.median():.2f}")
    print(f"  multi-domain clusters (2+ distinct families): "
          f"{(arch.n_distinct_families >= 2).sum():,}")
    top = hits.pfam_name.value_counts().head(8)
    print("  commonest families: " + ", ".join(f"{k} ({v})" for k, v in top.items()))
    print(f"\ndone in {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
