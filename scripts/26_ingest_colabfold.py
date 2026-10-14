"""26 - take in the structures that were folded elsewhere with ColabFold.

Step 24 defers any protein it cannot fold on this machine: too long for the card (`--max-len`) or
over the per-protein time cap (`--max-seconds`). Each one is written to
`Structures*/colabfold/<name>.fasta` and listed in `deferred_for_colabfold.tsv`.

Fold those FASTA files with ColabFold, drop whatever PDB it produces into
`Structures*/colabfold/returned/`, and run this script. It matches each file to its cluster by the
cluster id in the file name, in either spelling (`CFR-115601` or `CFR115601`), searching
subfolders too, so ColabFold's own layout and long names need no tidying, checks the residue count against the representative sequence, and writes the same
two files every other structure has:

  pdb/<name>.pdb            B-factor = pLDDT, taken from the ColabFold model as it comes
  pdb_variation/<name>.pdb  B-factor = the step-23 per-position variation, so the viewer can colour
                            by either

It then turns the cluster's `deferred_*` row in `structure_summary.tsv` into a `predicted` row, with
`source=colabfold`. Re-run step 25 afterwards to build the viewer pages for the new structures.

Usage: python3 scripts/26_ingest_colabfold.py [--group G] [--set amr|stress|virulence]
"""
import re
import sys

import numpy as np
import pandas as pd

from amr_common import (CLUSTER_PREFIX, PANGENOME, STRUCTURES, VARIATION, read_tsv, set_dir,
                        write_tsv)

STRUCTURES = set_dir(STRUCTURES)
VARIATION = set_dir(VARIATION)
# ColabFold and the download often drop the hyphen (CFR-115601 -> CFR115601), and write each
# model into its own folder, so match either spelling and search subfolders.
CLUSTER_ID = re.compile(rf"({CLUSTER_PREFIX}-?\d+)", re.I)


def count_ca(path):
    return sum(1 for line in open(path) if line.startswith("ATOM") and line[12:16].strip() == "CA")


def read_plddt(path):
    """One B-factor per residue, in residue order. ColabFold writes pLDDT there."""
    out = {}
    for line in open(path):
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            out[int(line[22:26])] = float(line[60:66])
    return np.array([out[k] for k in sorted(out)])


def rewrite_bfactors(pdb_text, values):
    out = []
    for line in pdb_text.splitlines():
        if line.startswith(("ATOM", "HETATM")):
            resi = int(line[22:26])
            v = values[resi - 1] if 0 < resi <= len(values) else 0.0
            line = f"{line[:60]}{min(v, 999.99):6.2f}{line[66:]}"
        out.append(line)
    return "\n".join(out) + "\n"


def load_sequences(md5s):
    seqs, keep, cur = {}, set(md5s), None
    with open(PANGENOME / "proteins_nr.faa") as fh:
        for line in fh:
            if line.startswith(">"):
                h = line[1:].split()[0]
                cur = h if h in keep else None
                if cur:
                    seqs[cur] = []
            elif cur:
                seqs[cur].append(line.strip())
    return {k: "".join(v).rstrip("*").upper() for k, v in seqs.items()}


def main():
    returned = STRUCTURES / "colabfold" / "returned"
    if not returned.exists():
        raise SystemExit(f"no folder {returned} - nothing to ingest")
    files = sorted(p for p in returned.rglob("*.pdb"))
    if not files:
        raise SystemExit(f"no .pdb files in {returned}")

    summ = read_tsv(STRUCTURES / "structure_summary.tsv")
    deferred = summ[summ.status.astype(str).str.startswith("deferred")]
    print(f"{len(files)} PDB file(s) returned; {len(deferred)} cluster(s) are waiting for one")
    if deferred.empty:
        raise SystemExit("nothing in the summary is deferred - is --set right?")

    want = dict(zip(deferred.cluster_id, deferred.name))
    seqs = load_sequences(set(deferred.representative_md5))
    md5 = dict(zip(deferred.cluster_id, deferred.representative_md5))

    var = None
    if (VARIATION / "position_variation.tsv").exists():
        v = pd.read_csv(VARIATION / "position_variation.tsv", sep="\t",
                        usecols=["cluster_id", "position", "var_frac_lineages"])
        var = {k: g.sort_values("position").var_frac_lineages.to_numpy() * 100
               for k, g in v.groupby("cluster_id")}

    # one file per cluster: if ColabFold returned several ranks, keep the best-scoring one
    best = {}
    for f in files:
        m = CLUSTER_ID.search(f.name)
        if not m:
            print(f"  {f.name}: no CFR- id in the file name, skipped")
            continue
        cid = m.group(1).upper()
        if "-" not in cid:                      # CFR115601 -> CFR-115601
            cid = f"{CLUSTER_PREFIX}-{cid[len(CLUSTER_PREFIX):]}"
        if cid not in want:
            print(f"  {f.name}: {cid} is not deferred in this set, skipped")
            continue
        mean = float(read_plddt(f).mean())
        if cid not in best or mean > best[cid][1]:
            best[cid] = (f, mean)

    rows, bad = [], []
    for cid, (f, mean) in sorted(best.items()):
        name, seq = want[cid], seqs[md5[cid]]
        n_res = count_ca(f)
        if n_res != len(seq):
            bad.append(f"{name}: {n_res} residues in {f.name}, expected {len(seq)}")
            continue
        text = f.read_text()
        plddt = read_plddt(f)
        (STRUCTURES / "pdb" / f"{name}.pdb").write_text(rewrite_bfactors(text, plddt))
        v = var.get(cid) if var else None
        if v is None or len(v) != len(seq):
            v = np.zeros(len(seq))
        (STRUCTURES / "pdb_variation" / f"{name}.pdb").write_text(rewrite_bfactors(text, v))
        row = summ[summ.cluster_id == cid].iloc[0].to_dict()
        row.update({"status": "predicted", "source": "colabfold", "defer_reason": "",
                    "mean_plddt": round(mean, 2), "min_plddt": round(float(plddt.min()), 2),
                    "frac_plddt_ge70": round(float((plddt >= 70).mean()), 4),
                    "max_variation_pct": round(float(np.max(v)), 2), "device": "colabfold",
                    "colabfold_file": f.name})
        rows.append(row)
        print(f"  {name}: {n_res} residues, mean pLDDT {mean:.1f}  <- {f.name}")

    if bad:
        print("\nNOT ingested - the model does not match the representative sequence:")
        for b in bad:
            print(f"  {b}")

    if rows:
        new = pd.DataFrame(rows)
        summ = pd.concat([summ[~summ.cluster_id.isin(new.cluster_id)], new]).sort_values("cluster_id")
        write_tsv(summ, "structure_summary.tsv", STRUCTURES)

    still = summ[summ.status.astype(str).str.startswith("deferred")]
    print(f"\ningested {len(rows)}; still waiting for ColabFold: {len(still)}"
          f"{' - ' + ', '.join(still.name) if len(still) else ''}")
    if rows:
        print(f"now re-run step 25 to build their viewer pages:\n"
              f"  tools/viz_env/bin/python scripts/25_structure_views.py --set <set>")


if __name__ == "__main__":
    main()
