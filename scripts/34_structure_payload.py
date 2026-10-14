"""34 - the structure payload: 808 backbone traces, small enough to inline.

The dashboard draws its own structures on a canvas rather than linking out to the 1.2 GB of py3Dmol
pages from step 25, so the pages stay self-contained and portable. All it needs is the CA trace.

THE PACKING, per chain:

  seq    one letter per residue
  coords base64-varint zigzag of the deltas between consecutive CA atoms, x/y/z interleaved, at
         0.1 A. Consecutive CA atoms are ~3.8 A apart, so a delta is one or two characters where an
         absolute coordinate would be five or six. This halves the coordinate cost.
  plddt  one printable character per residue, the model's own confidence
  var    the same, the step-23 per-position variation, and ABSENT where a cluster has none - 369 of
         the 808 chains are all zero, and a string of zeros is 86,000 wasted characters

Two structures need care. `clbB` has no whole-chain model at all: it was folded in three overlapping
modules, so its record carries the module ranges and three chain references instead of one. And the
resistance-set summary predates the `status` column, so its absence means "predicted", not missing.

Usage: python3 scripts/34_structure_payload.py [--group G]
"""
import importlib.util
import os
import time

import pandas as pd

from amr_common import STRUCTURES, read_tsv

_spec = importlib.util.spec_from_file_location(
    "dash", os.path.join(os.path.dirname(os.path.abspath(__file__)), "30_dashboard_common.py"))
dash = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(dash)

SETS = [("amr", ""), ("stress", "_stress"), ("virulence", "_virulence")]
AA3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
       "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
       "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}


def suffixed(base, suffix):
    return base if suffix == "" else base.with_name(base.name + suffix)


def read_ca(path):
    """CA atoms only: (x, y, z), one-letter residue, B-factor."""
    xyz, seq, b = [], [], []
    with open(path) as fh:
        for line in fh:
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                xyz.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
                seq.append(AA3.get(line[17:20].strip(), "X"))
                b.append(float(line[60:66]))
    return xyz, "".join(seq), b


def pack_coords(xyz):
    """Delta-encoded at 0.1 A, x/y/z interleaved, through the shared varint codec."""
    out, prev = [], (0, 0, 0)
    for c in xyz:
        v = [int(round(c[i] * 10)) for i in range(3)]
        out.append("".join(dash.vint(dash.zig(v[i] - prev[i])) for i in range(3)))
        prev = tuple(v)
    return "".join(out)


def main():
    t0 = time.time()
    chains, recs, missing, zero_var = [], {}, [], 0
    d = dash.Dicts()
    total_res = 0

    for name, suffix in SETS:
        folder = suffixed(STRUCTURES, suffix)
        summ = read_tsv(folder / "structure_summary.tsv")
        if "status" not in summ.columns:
            summ["status"] = "predicted"          # the amr summary predates the column
        for col in ("source", "module", "module_start", "module_end", "module_domains"):
            if col not in summ.columns:
                summ[col] = ""

        for r in summ.itertuples():
            pdb = folder / "pdb" / f"{r.name}.pdb"
            if not pdb.exists():
                # the parent row of a protein folded by module has no chain of its own, by design
                if str(r.status) != "folded_by_domain":
                    missing.append(r.name)
                continue
            xyz, seq, plddt = read_ca(pdb)
            if not xyz:
                missing.append(r.name)
                continue
            vfile = folder / "pdb_variation" / f"{r.name}.pdb"
            var = read_ca(vfile)[2] if vfile.exists() else []
            has_var = any(v > 0 for v in var)
            if not has_var:
                zero_var += 1

            idx = len(chains)
            chains.append({
                "c": dash.cluster_num(r.cluster_id), "m": dash.cluster_module(r.cluster_id),
                "n": len(xyz), "name": r.name, "seq": seq,
                "xyz": pack_coords(xyz), "pl": dash.qual_string(plddt),
                **({"var": dash.qual_string(var)} if has_var else {}),
                "set": {"amr": 0, "stress": 1, "virulence": 2}[name],
                "src": d.enc("source", r.source or "esmfold"),
                "frag": 1 if str(r.is_fragment) == "True" else 0,
            })
            total_res += len(xyz)

            cn = dash.cluster_num(r.cluster_id)
            rec = recs.setdefault(cn, {"chains": [], "modules": []})
            rec["chains"].append(idx)
            if dash.cluster_module(r.cluster_id):
                rec["modules"].append([int(float(r.module_start)), int(float(r.module_end)), idx,
                                       d.enc("moddom", r.module_domains)])

        # the parent row of a module-folded protein: record its true length and its module ranges
        for r in summ[summ.status.astype(str) == "folded_by_domain"].itertuples():
            cn = dash.cluster_num(r.cluster_id)
            rec = recs.setdefault(cn, {"chains": [], "modules": []})
            rec["whole"] = int(float(r.aa_len))
            rec["by_domain"] = 1

    for rec in recs.values():
        rec["modules"].sort()
        if not rec["modules"]:
            rec.pop("modules")

    payload = {"built": time.strftime("%Y-%m-%d"), "dicts": d.as_dict(),
               "chains": chains, "recs": recs}
    mb = dash.write_payload(payload, "structures.json")

    print("\nchecks")
    print(f"  chains: {len(chains)}; residues: {total_res:,}")
    print(f"  clusters represented: {len(recs)} "
          f"(expected 806); folded by module: {sum(1 for r in recs.values() if r.get('by_domain'))}")
    print(f"  chains with variation: {len(chains) - zero_var}; without: {zero_var}")
    by_src = pd.Series([c["src"] for c in chains]).value_counts().to_dict()
    names = d.as_dict().get("source", [])
    print("  sources: " + ", ".join(f"{names[k]} {v}" for k, v in by_src.items()))
    if missing:
        print(f"  MISSING PDB for {len(missing)}: {missing[:5]}")
    else:
        print("  every non-module-parent row has a PDB")
    print(f"  size {mb:.2f} MB ({1e6 * mb / max(1, total_res):.2f} bytes per residue)")
    print(f"\ndone in {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
