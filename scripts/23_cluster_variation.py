"""23 - per-position variation of every resistance-associated protein cluster against one reference.

For every anchor cluster (Islands/anchor_clusters.tsv) with at least 2 distinct protein sequences:

1. Reference = the most common FULL-LENGTH sequence, counted in lineages. Full length means at least
   90% of the copy-weighted median length of the cluster, so a common truncated form cannot win.
2. Copies shorter than 90% of the reference are FRAGMENTS (mostly genes cut by a contig edge). They
   are aligned too, but only feed the coverage track: a cut-off gene must never read as a deletion.
3. MAFFT aligns the distinct full-length sequences (never all copies: identical proteins are one
   sequence with a weight); fragments are added with --addfragments --keeplength.
4. Every column is projected onto reference coordinates and, per position, the share of copies
   carrying a substitution, a deletion or an insertion after it is computed twice:
     * per LINEAGE (default): a lineage is a Pathogen Detection SNP cluster (PDS_acc), or the isolate
       itself when it is in no SNP cluster. Each lineage weighs 1, split equally over its copies, so
       an outbreak of 500 genomes counts once - the same rule the island types use.
     * per protein COPY: the raw count, kept to show how much clonal expansion inflates it.
   Leading / trailing gaps are 'not covered', not deletions.
5. The reference is named after its closest catalog protein (MMseqs2 against AMRProt.fa), and for
   clusters with catalog point mutations (gyrA_T83I ...) the catalog reference protein is aligned to
   the cluster reference so the known sites are placed and labelled in catalog numbering.

Outputs in <group>/Variation/:
  cluster_variation_summary.tsv   one row per cluster
  position_variation.tsv          one row per (cluster, reference position)
  alignments/<cluster>.fasta      the alignment; header = md5|copies|lineages|full|fragment
  figures/<cluster>__<family>.png/.svg

Usage: python3 scripts/23_cluster_variation.py [--group G] [--refresh] [--only CFR-000010,...]
  --refresh   re-extract the protein copies from the 1.8 GB gene index (cached otherwise)
  --set       amr (default), stress or virulence; output goes to Variation[_<set>]/
"""
import math
import re
import shutil
import subprocess
import sys
import time
from collections import Counter, defaultdict

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from Bio import Align
from Bio.Align import substitution_matrices

from amr_common import (ISLANDS, OUT, PANGENOME, ROOT, SET, TAXGROUP, VARIATION, read_tsv, ref_dir,
                        set_dir, write_tsv)
from cluster_sets import load_cluster_set

VARIATION = set_dir(VARIATION)

FULL_LENGTH = 0.90          # of the reference length; below = fragment
REF_CANDIDATE = 0.90        # of the copy-weighted median length; below cannot be reference
ROW_RESIDUES = 100          # residues per row of the figure
N_PEAK_LABELS = 3

# validated with the dataviz skill's validate_palette.js (light surface #fcfcfb): all checks pass
C_SUB = "#2a78d6"           # substitution
C_DEL = "#eb6834"           # deletion
INK = "#1f1f1e"
INK_2 = "#52514e"
MUTED = "#a3a29c"
SURFACE = "#fcfcfb"
GRID = "#e6e5e0"
COV_FULL = "#b5b4ae"
COV_CMAP = matplotlib.colors.LinearSegmentedColormap.from_list("cov", [SURFACE, COV_FULL])


# ---------------------------------------------------------------------------- data

def load_clusters():
    cl = load_cluster_set()
    print(f"cluster set '{SET}': {len(cl)} clusters; with >=2 distinct sequences: {(cl.n_seqs >= 2).sum()}")
    cl = cl[cl.n_seqs >= 2].copy()
    raw = pd.read_csv(PANGENOME / "clusters_raw.tsv", sep="\t", header=None, names=["rep_md5", "md5"],
                      dtype=str)
    raw = raw[raw.rep_md5.isin(set(cl.rep_md5))]
    members = raw.merge(cl[["cluster_id", "rep_md5"]], on="rep_md5")[["cluster_id", "md5"]]
    return cl, members


def load_copies(members, refresh):
    """Every protein copy (gene) whose sequence belongs to one of the clusters, with its lineage."""
    VARIATION.mkdir(parents=True, exist_ok=True)
    cache = VARIATION / "_protein_copies.tsv"
    if refresh or not cache.exists():
        t0 = time.time()
        keys = VARIATION / "_md5_keys.txt"
        keys.write_text("\n".join(members.md5) + "\n")
        # awk streams the 1.8 GB index far faster than pandas; column 16 = protein_md5
        cmd = (f"zcat '{PANGENOME / 'gene_index.tsv.gz'}' | "
               f"awk -F'\\t' 'NR==FNR{{k[$1];next}} FNR==1{{print $1\"\\t\"$2\"\\t\"$16;next}} "
               f"($16 in k){{print $1\"\\t\"$2\"\\t\"$16}}' '{keys}' - > '{cache}'")
        subprocess.run(cmd, shell=True, check=True, executable="/bin/bash")
        keys.unlink()
        print(f"extracted protein copies from the gene index in {time.time() - t0:.0f} s")
    copies = read_tsv(cache).rename(columns={"protein_md5": "md5"})
    prof = read_tsv(OUT / "isolate_resistance_profile.tsv", usecols=["target_acc", "PDS_acc"])
    copies = copies.merge(prof, on="target_acc", how="left")
    copies["lineage"] = np.where(copies.PDS_acc.fillna("") != "", copies.PDS_acc, copies.target_acc)
    copies = copies.merge(members, on="md5")
    return copies


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
    return {k: "".join(v).rstrip("*") for k, v in seqs.items()}


def read_fasta(path):
    out, name = {}, None
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            name = line[1:].split()[0]
            out[name] = []
        elif name:
            out[name].append(line.strip())
    return {k: "".join(v) for k, v in out.items()}


def load_catalog():
    """AMRProt sequences by accession + point-mutation rows of this organism."""
    rd = ref_dir()
    prot, hdr = {}, {}
    for name, seq in read_fasta(rd / "AMRProt.fa").items():
        f = name.split("|")
        prot[f[0]] = seq
        hdr[name] = f
    cat = read_tsv(rd / "ReferenceGeneCatalog.txt")
    point = cat[(cat.subtype == "POINT") & (cat.whitelisted_taxa == TAXGROUP)]
    point_ref = dict(zip(point.allele, point.refseq_protein_accession))
    return prot, point_ref


# ---------------------------------------------------------------------------- per cluster

def weights(cc):
    """Lineage weight and copy count per distinct sequence of one cluster."""
    per_lin = cc.groupby("lineage").size()
    w = cc.lineage.map(1.0 / per_lin)
    return (cc.assign(w=w).groupby("md5").agg(copies=("gene_id", "size"), lineages=("w", "sum")))


def pick_reference(wt, seqs, wt_filter=None):
    lens = pd.Series({m: len(seqs[m]) for m in wt.index})
    order = lens.sort_values()
    cum = wt.copies.reindex(order.index).cumsum()
    med = order[cum >= cum.iloc[-1] / 2].iloc[0]
    cand = wt[lens >= REF_CANDIDATE * med]
    # clusters with catalog point mutations: the reference must be wild type at every site, otherwise
    # the resistant residue becomes the baseline and wild type reads as 'variation'
    wt_ok = wt_filter(list(cand.index)) if wt_filter else list(cand.index)
    if wt_ok:
        cand = cand.loc[wt_ok]
    ref = cand.sort_values(["lineages", "copies"], ascending=False).index[0]
    return ref, lens


def mafft(in_fa, out_fa, extra=()):
    cmd = ["mafft", "--quiet", "--anysymbol", "--thread", "8", *extra]
    with open(out_fa, "w") as fh:
        subprocess.run(cmd + [str(in_fa)] if not extra or extra[0] != "--addfragments" else cmd,
                       stdout=fh, check=True)


def align(cid, full, frags, seqs, tmp):
    """Alignment dict md5 -> gapped string. The reference is first in `full`."""
    fa = tmp / f"{cid}_full.fa"
    fa.write_text("".join(f">{m}\n{seqs[m]}\n" for m in full))
    aln_fa = tmp / f"{cid}_full.aln"
    if len(full) > 1:
        mafft(fa, aln_fa, ("--auto",))
    else:
        shutil.copy(fa, aln_fa)
    if frags:
        fr = tmp / f"{cid}_frag.fa"
        fr.write_text("".join(f">{m}\n{seqs[m]}\n" for m in frags))
        out = tmp / f"{cid}_all.aln"
        with open(out, "w") as fh:
            subprocess.run(["mafft", "--quiet", "--anysymbol", "--thread", "8", "--addfragments", str(fr),
                            "--keeplength", str(aln_fa)], stdout=fh, check=True)
        aln_fa = out
    aln = {k: v.upper() for k, v in read_fasta(aln_fa).items()}
    missing = set(full) | set(frags)
    missing -= set(aln)
    if missing:
        raise SystemExit(f"{cid}: {len(missing)} sequences missing from the alignment")
    return aln


def profile(ref, aln, wt, full_set):
    """Per reference position: coverage (all copies) and substitution / deletion / insertion shares."""
    ref_row = aln[ref]
    ref_cols = [i for i, c in enumerate(ref_row) if c != "-"]
    L = len(ref_cols)
    col_to_pos = np.full(len(ref_row), -1)
    for p, c in enumerate(ref_cols):
        col_to_pos[c] = p
    # the reference position an insertion column follows (-1 = before the first residue)
    after = np.maximum.accumulate(np.where(col_to_pos >= 0, col_to_pos, -1))

    z = lambda: np.zeros(L)
    cov_l, cov_c = z(), z()                       # all copies incl. fragments
    den_l, den_c = z(), z()                       # full-length copies covering the position
    sub_l, sub_c, del_l, del_c, ins_l, ins_c = z(), z(), z(), z(), z(), z()
    alt = [Counter() for _ in range(L)]
    ref_residues = [ref_row[c] for c in ref_cols]
    identity = {}

    for m, row in aln.items():
        wl, wc = wt.at[m, "lineages"], wt.at[m, "copies"]
        res_idx = [i for i, ch in enumerate(row) if ch != "-"]
        first, last = res_idx[0], res_idx[-1]
        covered = np.array([first <= c <= last for c in ref_cols])
        cov_l += wl * covered
        cov_c += wc * covered
        if m not in full_set:
            continue
        den_l += wl * covered
        den_c += wc * covered
        same = n = 0
        for p, c in enumerate(ref_cols):
            if not covered[p]:
                continue
            ch, r = row[c], ref_residues[p]
            n += 1
            if ch == "-":
                del_l[p] += wl; del_c[p] += wc
                alt[p]["del"] += wl
            elif ch != r:
                sub_l[p] += wl; sub_c[p] += wc
                alt[p][ch] += wl
            else:
                same += 1
        ins_pos = {after[c] for c in range(first, last + 1) if col_to_pos[c] < 0 and row[c] != "-"}
        for p in ins_pos:
            if p >= 0:
                ins_l[p] += wl; ins_c[p] += wc
        identity[m] = same / n if n else np.nan

    safe = lambda a, b: np.divide(a, b, out=np.zeros_like(a), where=b > 0)
    tot_l, tot_c = wt.lineages.sum(), wt.copies.sum()
    df = pd.DataFrame({
        "position": np.arange(1, L + 1), "ref_residue": ref_residues,
        "coverage_lineages": cov_l / tot_l, "coverage_copies": cov_c / tot_c,
        "n_lineages_scored": den_l.round(3), "n_copies_scored": den_c.astype(int),
        "sub_frac_lineages": safe(sub_l, den_l), "del_frac_lineages": safe(del_l, den_l),
        "ins_after_frac_lineages": safe(ins_l, den_l),
        "sub_frac_copies": safe(sub_c, den_c), "del_frac_copies": safe(del_c, den_c),
        "ins_after_frac_copies": safe(ins_c, den_c),
    })
    df["var_frac_lineages"] = df.sub_frac_lineages + df.del_frac_lineages
    df["var_frac_copies"] = df.sub_frac_copies + df.del_frac_copies
    df["alternatives"] = [";".join(f"{k}:{v / den_l[p]:.4f}" for k, v in alt[p].most_common(5))
                          if den_l[p] else "" for p in range(L)]
    df["n_alternatives"] = [len(a) for a in alt]
    return df, identity


def map_catalog(ref_seq, cat_seq):
    """reference position (1-based) -> catalog position (1-based), by global alignment."""
    al = Align.PairwiseAligner()
    al.substitution_matrix = substitution_matrices.load("BLOSUM62")
    al.open_gap_score, al.extend_gap_score = -11, -1
    al.mode = "global"
    a = al.align(ref_seq, cat_seq)[0]
    m = {}
    for (r0, r1), (c0, c1) in zip(*a.aligned):
        for k in range(r1 - r0):
            m[r0 + k + 1] = c0 + k + 1
    return m


def closest_catalog(refs, tmp):
    """MMseqs2 top hit of each reference against AMRProt.fa."""
    q = tmp / "refs.fa"
    q.write_text("".join(f">{cid}\n{s}\n" for cid, s in refs.items()))
    out = tmp / "refs_vs_amrprot.tsv"
    subprocess.run(["mmseqs", "easy-search", str(q), str(ref_dir() / "AMRProt.fa"), str(out), str(tmp / "mm"),
                    "-s", "7.5", "--threads", "8", "-v", "1",
                    "--format-output", "query,target,pident,qcov,tcov,bits"], check=True)
    hits = pd.read_csv(out, sep="\t", header=None, names=["cluster_id", "target", "pident", "qcov", "tcov", "bits"])
    hits = hits.sort_values(["cluster_id", "bits"], ascending=[True, False]).drop_duplicates("cluster_id")
    f = hits.target.str.split("|")
    hits["catalog_accession"] = f.str[0]
    hits["catalog_symbol"] = f.str[3]
    hits["catalog_product"] = f.str[-1].str.replace("_", " ")
    hits["catalog_identity"] = hits.pident.round(1)
    hits["catalog_coverage"] = hits[["qcov", "tcov"]].min(axis=1).round(3)
    return hits.set_index("cluster_id")[["catalog_accession", "catalog_symbol", "catalog_product",
                                         "catalog_identity", "catalog_coverage"]]


SITE_RE = re.compile(r"^(.+)_([A-Z])(\d+)([A-Z]+|DEL|STOP)$")


def place_sites(point_symbols, ref_seq, pos, prot, point_ref):
    """Catalog point-mutation sites -> reference coordinates, with a wild-type check."""
    known, notes = {}, []
    for sym in filter(None, (point_symbols or "").split(";")):
        mt = SITE_RE.match(sym)
        acc = point_ref.get(sym)
        if not mt or acc not in prot:
            notes.append(f"{sym}:unplaced")
            continue
        cat_pos = int(mt.group(3))
        p = {v: k for k, v in map_catalog(ref_seq, prot[acc]).items()}.get(cat_pos)
        if p is None:
            notes.append(f"{sym}:unplaced")
            continue
        known[p] = cat_pos
        res = ref_seq[p - 1]
        notes.append(f"{sym}:ref_pos={p}:ref_residue={res}{'' if res == mt.group(2) else '(NOT wild type ' + mt.group(2) + ')'}")
    return known, notes


def wild_type_md5s(point_symbols, candidates, seqs, prot, point_ref):
    """Sequences carrying the catalog wild-type residue at every point-mutation site."""
    sites = []
    for sym in filter(None, (point_symbols or "").split(";")):
        mt = SITE_RE.match(sym)
        if mt and point_ref.get(sym) in prot:
            sites.append((prot[point_ref[sym]], int(mt.group(3)), mt.group(2)))
    if not sites:
        return list(candidates)
    ok = []
    for m in candidates:
        good = True
        for cat_seq, cat_pos, wt in sites:
            p = {v: k for k, v in map_catalog(seqs[m], cat_seq).items()}.get(cat_pos)
            if p is None or seqs[m][p - 1] != wt:
                good = False
                break
        if good:
            ok.append(m)
    return ok


# ---------------------------------------------------------------------------- figure

def flat_note(s):
    """Say why a figure is empty, so a flat plot is never mistaken for missing data."""
    if s["n_variable_positions"] > 0:
        return ""
    why = []
    if s["n_distinct_fragments"]:
        why.append(f"{s['n_distinct_fragments']} fragment sequence(s)")
    if s["n_insertion_sites"]:
        why.append(f"insertions at {s['n_insertion_sites']} site(s) (table only)")
    if s["n_distinct_full_length"] > 1 and not s["n_insertion_sites"]:
        why.append(f"length outside the reference (full-length range {s['full_length_len_range']} aa)")
    return ("\nNo substitution or deletion in any full-length copy: the other sequences differ only by "
            + "; ".join(why) + ".")


def nice_ceiling(x):
    for v in (0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 0.75, 1.0):
        if x <= v:
            return v
    return 1.0


def draw(row, pos, known, path_stem):
    L = len(pos)
    n_rows = math.ceil(L / ROW_RESIDUES)
    ymax = nice_ceiling(max(pos.var_frac_lineages.max(), pos.var_frac_copies.max(), 0.001) * 1.2)
    variable = (pos.var_frac_lineages > 0).to_numpy()

    peaks = pos[pos.var_frac_lineages > 0].nlargest(N_PEAK_LABELS, "var_frac_lineages").position.tolist()
    labelled = sorted(set(peaks) | set(known))

    W = 15.5
    h_head, h_let, h_bar, h_cov, h_gap, h_foot = 0.62 + 0.15 * (row.subtitle.count("\n") + 1), 0.2, 1.0, 0.07, 0.42, 0.45
    H = h_head + n_rows * (h_let + h_bar + h_cov + h_gap) + h_foot
    fig = plt.figure(figsize=(W, H), facecolor=SURFACE)
    left, right = 1.05, 0.25
    unit = (W - left - right) / ROW_RESIDUES

    fig.text(left / W, 1 - 0.28 / H, row.title, fontsize=11, color=INK, weight="bold", va="top")
    fig.text(left / W, 1 - 0.55 / H, row.subtitle, fontsize=7.5, color=INK_2, va="top", linespacing=1.5)

    y_top = H - h_head
    for r in range(n_rows):
        lo, hi = r * ROW_RESIDUES, min(L, (r + 1) * ROW_RESIDUES)
        seg = pos.iloc[lo:hi]
        w_in = (hi - lo) * unit
        x = seg.position.to_numpy()

        # reference letters
        y_let = y_top - h_let
        axl = fig.add_axes([left / W, y_let / H, w_in / W, h_let / H])
        axl.set_xlim(lo + 0.5, hi + 0.5); axl.set_ylim(0, 1); axl.axis("off")
        for p, res in zip(x, seg.ref_residue):
            v = variable[p - 1]
            axl.text(p, 0.35, res, ha="center", va="center", fontsize=6.3, family="monospace",
                     color=INK if v else MUTED, weight="bold" if v else "normal")

        # variation bars
        y_bar = y_let - h_bar
        ax = fig.add_axes([left / W, y_bar / H, w_in / W, h_bar / H], facecolor=SURFACE)
        ax.set_xlim(lo + 0.5, hi + 0.5); ax.set_ylim(0, ymax)
        ax.bar(x, seg.sub_frac_lineages, width=0.8, color=C_SUB, linewidth=0)
        ax.bar(x, seg.del_frac_lineages, bottom=seg.sub_frac_lineages, width=0.8, color=C_DEL, linewidth=0)
        ax.step(x, seg.var_frac_copies, where="mid", color=INK_2, linewidth=0.8)
        ax.set_yticks([0, ymax / 2, ymax])
        ax.set_yticklabels([f"{t:.0%}" if ymax >= 0.05 else f"{t:.1%}" for t in (0, ymax / 2, ymax)],
                           fontsize=6, color=INK_2)
        ax.yaxis.grid(True, color=GRID, linewidth=0.5); ax.set_axisbelow(True)
        ax.set_xticks([])
        for s in ("top", "right", "bottom"):
            ax.spines[s].set_visible(False)
        ax.spines["left"].set_color(GRID)
        ax.tick_params(axis="y", length=0)
        for p in labelled:
            if not lo < p <= hi:
                continue
            pr = pos.iloc[p - 1]
            top_alt = pr.alternatives.split(";")[0].split(":")[0] if pr.alternatives else ""
            num = known.get(p, p)
            txt = f"{pr.ref_residue}{num}{'' if top_alt == 'del' else top_alt}"
            if top_alt == "del":
                txt += " del"
            if pr.var_frac_lineages > 0:
                txt += f" {pr.var_frac_lineages:.0%}" if pr.var_frac_lineages >= 0.01 else " <1%"
            y = min(pr.var_frac_lineages, ymax) + ymax * 0.04
            if p in known:
                txt = "\u25BC " + txt + " (catalog site)"
            ax.text(p, y, txt, ha="center", va="bottom", fontsize=6.5, color=INK,
                    weight="bold" if p in known else "normal", clip_on=False)

        # coverage strip (all copies incl. fragments), grey sequential
        y_cov = y_bar - h_cov
        axc = fig.add_axes([left / W, y_cov / H, w_in / W, h_cov / H])
        axc.imshow(seg.coverage_lineages.to_numpy()[None, :], aspect="auto", cmap=COV_CMAP, vmin=0, vmax=1,
                   extent=(lo + 0.5, hi + 0.5, 0, 1))
        axc.set_yticks([])
        ticks = [t for t in range(10 * math.ceil((lo + 1) / 10), hi + 1, 10)]
        axc.set_xticks(ticks); axc.set_xticklabels(ticks, fontsize=6, color=INK_2)
        axc.tick_params(axis="x", length=2, color=MUTED, pad=1.5)
        for s in axc.spines.values():
            s.set_visible(False)
        if r == 0:
            ax.set_ylabel("differs from\nreference", fontsize=6.5, color=INK_2, rotation=0, ha="right",
                          va="center", labelpad=8)
            axc.set_ylabel("coverage", fontsize=6, color=INK_2, rotation=0, ha="right", va="center", labelpad=4)
        y_top = y_cov - h_gap

    handles = [Patch(color=C_SUB, label="substitution, share of lineages"),
               Patch(color=C_DEL, label="deletion, share of lineages"),
               Line2D([], [], color=INK_2, lw=0.8, label="any difference, share of protein copies"),
               Patch(color=COV_FULL, label="coverage strip: share of lineages whose copy reaches the "
                                            "position, fragments included (grey = all, pale = few)")]
    fig.legend(handles=handles, loc="lower left", bbox_to_anchor=(left / W, 0.08 / H), ncol=4, fontsize=6.8,
               frameon=False, handlelength=1.4, labelcolor=INK_2)
    fig.savefig(f"{path_stem}.png", dpi=150, facecolor=SURFACE)
    fig.savefig(f"{path_stem}.svg", facecolor=SURFACE)
    plt.close(fig)


# ---------------------------------------------------------------------------- main

def safe_name(s):
    s = re.sub(r"[^A-Za-z0-9._-]+", "_", s).strip("_")
    return s[:60] or "unlabelled"


def main():
    t0 = time.time()
    refresh = "--refresh" in sys.argv
    only = set(sys.argv[sys.argv.index("--only") + 1].split(",")) if "--only" in sys.argv else None

    cl, members = load_clusters()
    if only:
        cl = cl[cl.cluster_id.isin(only)]
    copies = load_copies(members, refresh)
    copies = copies[copies.cluster_id.isin(set(cl.cluster_id))]
    seqs = load_sequences(set(copies.md5))
    prot, point_ref = load_catalog()
    print(f"clusters to process: {len(cl)}; protein copies: {len(copies):,}; distinct sequences: {len(seqs):,}")

    tmp = VARIATION / "_tmp"
    if tmp.exists():
        shutil.rmtree(tmp)
    tmp.mkdir(parents=True)
    (VARIATION / "alignments").mkdir(exist_ok=True)
    (VARIATION / "figures").mkdir(exist_ok=True)

    by_cluster = dict(tuple(copies.groupby("cluster_id")))
    work = {}
    for c in cl.itertuples():
        cc = by_cluster[c.cluster_id]
        wt = weights(cc)
        ref, lens = pick_reference(wt, seqs, lambda cand: wild_type_md5s(c.point_symbols, cand, seqs, prot,
                                                                          point_ref))
        full = [m for m in wt.index if lens[m] >= FULL_LENGTH * lens[ref]]
        full = [ref] + [m for m in full if m != ref]
        frags = [m for m in wt.index if m not in set(full)]
        work[c.cluster_id] = (cc, wt, ref, lens, full, frags)

    closest = closest_catalog({cid: seqs[w[2]] for cid, w in work.items()}, tmp)

    summaries, positions = [], []
    for i, c in enumerate(cl.itertuples(), 1):
        cc, wt, ref, lens, full, frags = work[c.cluster_id]
        aln = align(c.cluster_id, full, frags, seqs, tmp)
        if aln[ref].replace("-", "") != seqs[ref].upper():
            raise SystemExit(f"{c.cluster_id}: reference changed during alignment")
        pos, identity = profile(ref, aln, wt, set(full))

        known, site_notes = place_sites(c.point_symbols, seqs[ref], pos, prot, point_ref)
        pos["catalog_site"] = pos.position.map(lambda p: known.get(p, ""))

        with open(VARIATION / "alignments" / f"{c.cluster_id}.fasta", "w") as fh:
            for m in [ref] + [m for m in aln if m != ref]:
                kind = "full" if m in set(full) else "fragment"
                tag = "REFERENCE|" if m == ref else ""
                fh.write(f">{tag}{m}|copies={wt.at[m, 'copies']}|lineages={wt.at[m, 'lineages']:.2f}|{kind}\n"
                         f"{aln[m]}\n")

        fw = wt.loc[full]
        ident = pd.Series(identity).reindex(fw.index)
        frag_copies = int(wt.loc[frags, "copies"].sum()) if frags else 0
        point_genes = sorted({x.rsplit("_", 1)[0] for x in (c.point_symbols or "").split(";") if x})
        fam = c.amr_gene_families or ";".join(point_genes) or c.product_top or "unlabelled"
        cat = closest.loc[c.cluster_id] if c.cluster_id in closest.index else None
        s = {
            "cluster_id": c.cluster_id, "anchor_tier": c.anchor_tier, "amr_gene_families": c.amr_gene_families,
            "amr_classes": c.amr_classes, "product_top": c.product_top,
            "n_distinct_seqs": len(wt), "n_distinct_full_length": len(full), "n_distinct_fragments": len(frags),
            "n_copies": int(wt.copies.sum()), "n_fragment_copies": frag_copies,
            "fragment_copy_share": round(frag_copies / wt.copies.sum(), 4),
            "n_lineages": int(round(wt.lineages.sum())),
            "reference_md5": ref, "reference_len": lens[ref],
            "reference_copies": int(wt.at[ref, "copies"]),
            "reference_lineage_share": round(wt.at[ref, "lineages"] / wt.lineages.sum(), 4),
            "catalog_symbol": cat.catalog_symbol if cat is not None else "",
            "catalog_accession": cat.catalog_accession if cat is not None else "",
            "catalog_identity": cat.catalog_identity if cat is not None else np.nan,
            "catalog_coverage": cat.catalog_coverage if cat is not None else np.nan,
            "mean_identity_to_reference_lineages": round(float(np.average(ident, weights=fw.lineages)), 5),
            "min_identity_to_reference": round(float(ident.min()), 4),
            "n_variable_positions": int((pos.var_frac_lineages > 0).sum()),
            "n_positions_ge_1pct": int((pos.var_frac_lineages >= 0.01).sum()),
            "n_positions_ge_5pct": int((pos.var_frac_lineages >= 0.05).sum()),
            "n_insertion_sites": int((pos.ins_after_frac_lineages > 0).sum()),
            "full_length_len_range": f"{min(lens[m] for m in full)}-{max(lens[m] for m in full)}",
            "n_positions_ge_1pct_copies": int((pos.var_frac_copies >= 0.01).sum()),
            "top_position": int(pos.var_frac_lineages.idxmax() + 1) if pos.var_frac_lineages.max() > 0 else "",
            "top_position_frac_lineages": round(pos.var_frac_lineages.max(), 4),
            "catalog_sites": ";".join(site_notes),
            "figure": f"figures/{c.cluster_id}__{safe_name(fam.split(';')[0])}",
        }
        summaries.append(s)
        pos.insert(0, "cluster_id", c.cluster_id)
        positions.append(pos)

        cat_txt = (f"closest catalog protein {s['catalog_symbol']} ({s['catalog_accession']}), "
                   f"{s['catalog_identity']}% identity" if cat is not None else "no catalog hit")
        row = pd.Series({
            "title": f"{c.cluster_id}  \u00b7  {fam.replace(';', ', ')}  \u00b7  {c.amr_classes or 'point mutation'}"
                     f"  \u00b7  {c.anchor_tier}",
            "subtitle": (f"Reference: the most common full-length sequence, {lens[ref]} aa, carried by "
                         f"{s['reference_lineage_share']:.0%} of lineages; {cat_txt}.\n"
                         f"{s['n_distinct_seqs']} distinct sequences, {s['n_copies']:,} protein copies, "
                         f"{s['n_lineages']:,} lineages (SNP clusters, or single isolates outside one). "
                         f"{frag_copies:,} {'copy' if frag_copies == 1 else 'copies'} shorter than {FULL_LENGTH:.0%} of the reference are fragments: "
                         f"they count in the coverage strip only.\n"
                         f"Bars count lineages, so one outbreak counts once. Letters in bold vary in at least "
                         f"one lineage." + flat_note(s)),
        })
        draw(row, pos, known, VARIATION / s["figure"])
        if i % 25 == 0 or i == len(cl):
            print(f"  {i}/{len(cl)} clusters  ({time.time() - t0:.0f} s)", flush=True)

    shutil.rmtree(tmp)
    summ = pd.DataFrame(summaries)
    write_tsv(summ, "cluster_variation_summary.tsv", VARIATION)
    allpos = pd.concat(positions)
    for col in [c for c in allpos.columns if "frac" in c or "coverage" in c]:
        allpos[col] = allpos[col].round(5)
    write_tsv(allpos, "position_variation.tsv", VARIATION)

    # ---- checks
    print("\nchecks")
    ann_members = cl.set_index("cluster_id").n_members
    got = copies.groupby("cluster_id").size().reindex(ann_members.index)
    print(f"  copies per cluster = cluster n_members: {(got == ann_members).sum()}/{len(ann_members)}")
    print(f"  every cluster processed: {len(summ) == len(cl)}")
    g = summ[summ.catalog_sites.str.contains("gyrA", na=False)]
    for r in g.itertuples():
        p = allpos[allpos.cluster_id == r.cluster_id]
        top = p.nlargest(3, "var_frac_lineages")
        print(f"  {r.cluster_id} gyrA sites: {r.catalog_sites}")
        print("    top positions (ref pos, catalog site, share of lineages): " +
              ", ".join(f"{t.position}/{t.catalog_site or '-'}/{t.var_frac_lineages:.3f}" for t in top.itertuples()))
    print(f"\ndone in {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
