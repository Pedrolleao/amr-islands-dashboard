"""36 - decode the payloads back and check them against the tables they came from.

The dashboard's data is packed hard enough that a codec bug would not look like a crash. It would
look like a gene drawn in the wrong place, or a structure that is subtly the wrong shape, and nobody
would notice. So the packing is checked by decoding it again in Python and comparing against the
source TSVs and PDBs, rather than by trusting the encoder.

Usage: python3 scripts/36_verify_payloads.py [--group G] [--islands N] [--structures N]
"""
import importlib.util
import json
import os
import random
import sys

from amr_common import DASHBOARD, ISLANDS, PANGENOME, STRUCTURES, read_tsv

_spec = importlib.util.spec_from_file_location(
    "dash", os.path.join(os.path.dirname(os.path.abspath(__file__)), "30_dashboard_common.py"))
dash = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(dash)

SETS = ["", "_stress", "_virulence"]


def decode_island(packed, bp, flag_label, flag_pseudo, n_roles):
    """The Python mirror of what the page's JavaScript does, so both are checked against the same
    source of truth rather than against each other."""
    vals = dash.varints(packed)
    i = 0
    pos = vals[i] * bp
    i += 1
    genes, prev_end = [], vals[0]
    while i < len(vals):
        gap = dash.unzig(vals[i]); i += 1
        width = vals[i]; i += 1
        cs = vals[i]; i += 1
        flag = vals[i]; i += 1
        start = (prev_end + gap) * bp
        end = start + width * bp
        label_idx = None
        if flag & flag_label:
            label_idx = vals[i] - 1; i += 1
        genes.append({"start": start, "end": end,
                      "cluster": (cs >> 1) - 1, "strand": 1 if cs & 1 else -1,
                      "role": flag & (n_roles - 1 if n_roles & (n_roles - 1) == 0 else 15),
                      "pseudo": bool(flag & flag_pseudo), "label_idx": label_idx})
        prev_end = prev_end + gap + width
    return genes


def decode_coords(packed):
    vals = dash.varints(packed)
    out, cur = [], [0, 0, 0]
    for j in range(0, len(vals), 3):
        for k in range(3):
            cur[k] += dash.unzig(vals[j + k])
        out.append((cur[0] / 10, cur[1] / 10, cur[2] / 10))
    return out


def check_islands(n_sample):
    gm = json.loads((DASHBOARD / "genemaps.json").read_text())
    bp, roles = gm["bp"], gm["roles"]
    fl, fp = gm["flag_label"], gm["flag_pseudo"]
    keys = random.sample(sorted(gm["g"]), min(n_sample, len(gm["g"])))
    want = {dash.ISLAND_PREFIX + k.zfill(6) for k in keys}

    src = read_tsv(ISLANDS / "island_gene_table.tsv",
                   usecols=["island_id", "pos_in_island", "cluster_id", "start", "end", "strand",
                            "draw_role", "is_pseudo"])
    src = src[src.island_id.isin(want)].copy()
    src["pos_in_island"] = src.pos_in_island.astype(int)      # numeric, or "10" sorts before "2"
    src = src.sort_values(["island_id", "pos_in_island"])
    # the payload must list genes in contig order, which a string sort silently would not

    bad_n = bad_coord = bad_role = bad_cluster = bad_pseudo = bad_order = 0
    worst = 0
    for key in keys:
        iid = dash.ISLAND_PREFIX + key.zfill(6)
        got = decode_island(gm["g"][key], bp, fl, fp, len(roles))
        if any(got[k]["start"] > got[k + 1]["start"] for k in range(len(got) - 1)):
            bad_order += 1                       # genes must come out in contig order
        rows = src[src.island_id == iid]
        if len(got) != len(rows):
            bad_n += 1
            continue
        for g, r in zip(got, rows.itertuples()):
            worst = max(worst, abs(g["start"] - int(r.start)), abs(g["end"] - int(r.end)))
            if abs(g["start"] - int(r.start)) > bp / 2 or abs(g["end"] - int(r.end)) > bp / 2:
                bad_coord += 1
            if roles[g["role"]] != str(r.draw_role):
                bad_role += 1
            cid = "" if g["cluster"] < 0 else f"{dash.CLUSTER_PREFIX}-{gm['cluster'][g['cluster']]:06d}"
            if cid != ("" if str(r.cluster_id) in ("", "nan") else str(r.cluster_id)):
                bad_cluster += 1
            if g["pseudo"] != (str(r.is_pseudo) in ("1", "True", "1.0")):
                bad_pseudo += 1
    # the per-cluster AMR class arrays: every referenced cluster, against the annotation table
    ann = read_tsv(PANGENOME / "cluster_annotation.tsv",
                   usecols=["cluster_id", "amr_classes", "amr_subclasses"])
    ann = ann.set_index("cluster_id")
    bad_class = 0
    if len(gm["cclass"]) != len(gm["cluster"]) or len(gm["csub"]) != len(gm["cluster"]):
        bad_class += 1
    else:
        for c, ci, si in zip(gm["cluster"], gm["cclass"], gm["csub"]):
            cid = f"{dash.CLUSTER_PREFIX}-{c:06d}"
            want_c = str(ann.amr_classes.get(cid, "")) if cid in ann.index else ""
            want_s = str(ann.amr_subclasses.get(cid, "")) if cid in ann.index else ""
            got_c = gm["dicts"]["amrclass"][ci] if ci >= 0 else ""
            got_s = gm["dicts"]["amrsubclass"][si] if si >= 0 else ""
            if got_c != want_c or got_s != want_s:
                bad_class += 1
    print(f"  islands sampled: {len(keys)}")
    print(f"    cluster AMR class/subclass mismatches over all {len(gm['cluster']):,} referenced: {bad_class}")
    print(f"    gene count mismatches: {bad_n}")
    print(f"    islands whose genes are out of contig order: {bad_order}")
    print(f"    coordinate mismatches beyond {bp / 2:.0f} bp: {bad_coord} (worst {worst} bp)")
    print(f"    role mismatches: {bad_role}; cluster mismatches: {bad_cluster}; "
          f"pseudogene mismatches: {bad_pseudo}")
    return bad_n + bad_coord + bad_role + bad_cluster + bad_pseudo + bad_order + bad_class


def check_structures(n_sample):
    st = json.loads((DASHBOARD / "structures.json").read_text())
    chains = random.sample(st["chains"], min(n_sample, len(st["chains"])))
    bad_len = bad_xyz = bad_pl = 0
    worst_xyz = worst_pl = 0.0
    for c in chains:
        path = None
        for suffix in SETS:
            p = (STRUCTURES if suffix == "" else STRUCTURES.with_name(STRUCTURES.name + suffix))
            cand = p / "pdb" / f"{c['name']}.pdb"
            if cand.exists():
                path = cand
                break
        if path is None:
            bad_len += 1
            continue
        xyz, seq, plddt = read_ca_local(path)
        got = decode_coords(c["xyz"])
        if len(got) != len(xyz) or c["n"] != len(xyz) or len(c["seq"]) != len(xyz):
            bad_len += 1
            continue
        for a, b in zip(got, xyz):
            e = max(abs(a[0] - b[0]), abs(a[1] - b[1]), abs(a[2] - b[2]))
            worst_xyz = max(worst_xyz, e)
            if e > 0.06:
                bad_xyz += 1
        for a, b in zip(dash.qual_values(c["pl"]), plddt):
            worst_pl = max(worst_pl, abs(a - b))
            if abs(a - b) > 0.6:
                bad_pl += 1
    print(f"  structures sampled: {len(chains)}")
    print(f"    length mismatches: {bad_len}")
    print(f"    coordinates off by more than 0.06 A: {bad_xyz} (worst {worst_xyz:.3f} A)")
    print(f"    confidence off by more than 0.6: {bad_pl} (worst {worst_pl:.3f})")
    return bad_len + bad_xyz + bad_pl


def read_ca_local(path):
    xyz, seq, b = [], [], []
    with open(path) as fh:
        for line in fh:
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                xyz.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
                seq.append(line[17:20].strip())
                b.append(float(line[60:66]))
    return xyz, seq, b


def main():
    random.seed(0)
    n_isl = int(sys.argv[sys.argv.index("--islands") + 1]) if "--islands" in sys.argv else 100
    n_str = int(sys.argv[sys.argv.index("--structures") + 1]) if "--structures" in sys.argv else 20
    print("decoding the payloads and comparing against the source tables")
    bad = check_islands(n_isl) + check_structures(n_str)
    sizes = {p.name: p.stat().st_size / 1e6 for p in sorted(DASHBOARD.glob("*.json"))}
    print("\n  payload sizes: " + ", ".join(f"{k} {v:.2f} MB" for k, v in sizes.items()))
    print(f"  total {sum(sizes.values()):.2f} MB")
    print("\nPASS - every sampled record decodes back to its source" if not bad
          else f"\nFAIL - {bad} mismatches")
    raise SystemExit(1 if bad else 0)


if __name__ == "__main__":
    main()
