"""24 - ESMFold structure of the representative protein of every resistance-associated cluster.

One structure per anchor cluster (Islands/anchor_clusters.tsv), 411 for Citrobacter freundii:
  * the 253 clusters of step 23 use that step's reference (most common full-length sequence, wild type
    at catalog point-mutation sites);
  * the other 158 contain a single distinct sequence, which is that cluster's representative.

ESMFold runs locally (no MSA, nothing leaves the machine). The model is loaded with the language model
in half precision - in full precision it fills a 16 GB card and every long protein spills into system
memory. Proteins above ~900 residues still need a smaller chunk size, and the run falls back
automatically: GPU chunk 32 -> GPU chunk 8 -> CPU.

Two PDB files are written per cluster. They hold the same coordinates and differ in the B-factor column,
which is what molecular viewers colour by:
  pdb/           B-factor = pLDDT (0-100), the model's own confidence
  pdb_variation/ B-factor = share of lineages differing from the reference at that position (0-100),
                 from step 23 - colour by B-factor to see where the variable sites sit in the fold.
                 0 everywhere for the 158 single-sequence clusters: nothing varies in them.

Usage: python3 scripts/24_structure_prediction.py [--group G] [--only CFR-000010,...] [--cpu] [--redo]
  Resumable: clusters whose PDB already exists are skipped unless --redo.
  --set amr|stress|virulence   which cluster group; output goes to Structures[_<set>]/
  --max-len N       fold up to N residues (default 1500).
  --max-seconds N   drop any protein that takes longer than this (default 3600, one hour), and with
                    it every protein still queued - the queue is sorted by length, so the rest are
                    longer and would only be slower.
  --defer CFR-x,... drop these clusters without trying them, with --defer-reason to say why.
  --no-offload      keep the ESM-2 language model on the card for the whole fold (the old behaviour).

  By default the language model is evicted to the host as soon as it has produced the embeddings.
  It runs once per protein and the trunk, which does the recycling, never touches it, so the 6 GB it
  occupies is dead weight through the expensive part. Evicting it is bit-identical - the embeddings
  are still computed on the GPU, the weights are merely released afterwards - and measured on this
  16 GB card it cut silA (1,048 aa) from 328 s to 131 s at an unchanged pLDDT of 88.30, and turned
  sigA (1,287 aa) from over an hour with no result into 562 s.

  Nothing dropped is ever silently missing. Every one of them lands in deferred_for_colabfold.tsv,
  in deferred_for_colabfold.faa, and as its own FASTA in colabfold/, and carries
  status=deferred_too_long or deferred_too_slow with a reason in structure_summary.tsv. Fold them
  with ColabFold, put the PDB in colabfold/returned/, then run scripts/26_ingest_colabfold.py.
"""
import gc
import os
import re
import sys
import threading
import time

import numpy as np
import pandas as pd
import torch
from transformers import AutoTokenizer, EsmForProteinFolding

from amr_common import (DATA, ISLANDS, PANGENOME, SET, STRUCTURES, VARIATION, read_tsv, set_dir,
                        write_tsv)
from cluster_sets import load_cluster_set

STRUCTURES = set_dir(STRUCTURES)
VARIATION = set_dir(VARIATION)

LONG_AA = 900            # above this, start at the smaller chunk size
MAX_AA = 1500            # above this, ESMFold does not fit this 16 GB card: defer to local ColabFold
MAX_SECONDS = 3600       # a single protein may not take longer than this: defer it instead
FRAGMENT_OF_FAMILY = 0.80


def load_targets():
    """One representative sequence per anchor cluster."""
    d = load_cluster_set()
    var_file = VARIATION / "cluster_variation_summary.tsv"
    if var_file.exists():
        v = read_tsv(var_file, usecols=["cluster_id", "reference_md5", "catalog_symbol", "catalog_identity"])
        d = d.merge(v, on="cluster_id", how="left")
    else:
        d["reference_md5"] = d["catalog_symbol"] = d["catalog_identity"] = ""
    # step-23 reference where there is one, the single sequence otherwise
    d["md5"] = np.where(d.reference_md5.fillna("") != "", d.reference_md5, d.rep_md5)
    return d


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


def flag_fragments(d, seqs):
    """A representative much shorter than the longest cluster of the same gene family is a fragment."""
    d = d.copy()
    d["aa_len"] = d.md5.map(lambda m: len(seqs[m]))
    d["family"] = d.amr_gene_families.fillna("").str.split(";").str[0]
    point_gene = d.point_symbols.fillna("").str.split(";").str[0].str.rsplit("_", n=1).str[0]
    d["family"] = np.where(d.family != "", d.family, point_gene)
    longest = d[d.family != ""].groupby("family").aa_len.max()
    d["family_max_len"] = d.family.map(longest)
    d["is_fragment"] = (d.aa_len < FRAGMENT_OF_FAMILY * d.family_max_len).fillna(False)
    d["only_cluster_of_family"] = d.family.map(d[d.family != ""].family.value_counts()).eq(1).fillna(False)
    return d


def safe_name(s):
    return (re.sub(r"[^A-Za-z0-9._-]+", "_", s).strip("_") or "unlabelled")[:60]


def rewrite_bfactors(pdb, values):
    """Replace the B-factor column; `values` is one number per residue (1-based index)."""
    out = []
    for line in pdb.splitlines():
        if line.startswith(("ATOM", "HETATM")):
            resi = int(line[22:26])
            v = values[resi - 1] if 0 < resi <= len(values) else 0.0
            line = f"{line[:60]}{min(v, 999.99):6.2f}{line[66:]}"
        out.append(line)
    return "\n".join(out) + "\n"


def install_lm_offload(model):
    """Free the ESM-2 weights from the card once they have produced the embeddings.

    EsmForProteinFolding calls the language model exactly once per protein, then hands its output to
    the trunk, which does all the recycling and holds the L x L pair representation. The 6 GB of
    language-model weights sit unused through that, which on a 16 GB card is the difference between
    computing and spilling into host memory over PCIe.

    This does not change the result: the embeddings are still computed on the GPU with the same
    weights in the same precision, and only then are the weights moved off. Measured here, silA
    (1,048 aa) went from 328 s to 131 s with an identical mean pLDDT of 88.30.
    """
    original = model.compute_language_model_representations

    def offloading(esmaa):
        out = original(esmaa)
        model.esm.to("cpu")
        gc.collect()
        torch.cuda.empty_cache()
        return out

    model.compute_language_model_representations = offloading
    model._lm_offloaded = True


def fold(model, tok, seq, chunk, device):
    model.trunk.set_chunk_size(chunk)
    with torch.no_grad():
        ids = tok([seq], return_tensors="pt", add_special_tokens=False)["input_ids"].to(device)
        out = model(ids)
    plddt = out["plddt"][0, :, 1].float().cpu().numpy()          # CA
    if plddt.max() <= 1.0:
        plddt = plddt * 100
    return model.output_to_pdb(out)[0], plddt, float(out["ptm"])


def read_bfactors(path):
    """One B-factor per residue, in residue order (CA atoms)."""
    out = {}
    for line in open(path):
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            out[int(line[22:26])] = float(line[60:66])
    return np.array([out[k] for k in sorted(out)])


def recover_row(c, structures):
    """Rebuild the summary row of a cluster whose PDB exists but whose row was lost.

    A run interrupted after writing structures but before writing the table would otherwise leave
    those clusters out of structure_summary.tsv for good, since a resumed run skips them. Everything
    except ptm is in the two PDB files; ptm is not, and stays empty.
    """
    plddt = read_bfactors(structures / "pdb" / f"{c.name}.pdb")
    vfile = structures / "pdb_variation" / f"{c.name}.pdb"
    v = read_bfactors(vfile) if vfile.exists() else np.zeros(len(plddt))
    return {
        "status": "predicted",
        "cluster_id": c.cluster_id, "name": c.name, "family": c.family,
        "anchor_tier": c.anchor_tier, "amr_classes": c.amr_classes, "product_top": c.product_top,
        "n_members": c.n_members, "n_distinct_seqs": c.n_seqs, "prevalence": c.prevalence,
        "representative_md5": c.md5, "aa_len": c.aa_len,
        "from_step23_reference": isinstance(getattr(c, "reference_md5", ""), str) and bool(getattr(c, "reference_md5", "")),
        "catalog_symbol": getattr(c, "catalog_symbol", ""), "is_fragment": bool(c.is_fragment),
        "only_cluster_of_family": bool(c.only_cluster_of_family),
        "mean_plddt": round(float(plddt.mean()), 2), "min_plddt": round(float(plddt.min()), 2),
        "frac_plddt_ge70": round(float((plddt >= 70).mean()), 4), "ptm": "",
        "max_variation_pct": round(float(v.max()), 2), "device": "recovered_from_pdb",
    }


def defer_record(c, reason):
    """One row of the deferral tables, from a cluster tuple."""
    return {"cluster_id": c.cluster_id, "name": c.name, "family": c.family, "aa_len": int(c.aa_len),
            "md5": c.md5, "product_top": c.product_top, "n_members": c.n_members,
            "n_distinct_seqs": c.n_seqs, "anchor_tier": c.anchor_tier, "amr_classes": c.amr_classes,
            "is_fragment": bool(c.is_fragment), "reason": reason}


def write_deferred(records, seqs, structures):
    """The hand-off to ColabFold: one table, one combined FASTA, and one FASTA per protein.

    Anything ESMFold cannot do here lands in these files rather than going missing. Drop the PDB
    ColabFold returns into colabfold/returned/ and run scripts/26_ingest_colabfold.py.
    """
    if not records:
        return
    recs = sorted(records, key=lambda r: r["aa_len"])
    cols = ["cluster_id", "name", "family", "aa_len", "md5", "product_top", "n_members", "reason"]
    write_tsv(pd.DataFrame(recs)[cols], "deferred_for_colabfold.tsv", structures)
    (structures / "deferred_for_colabfold.faa").write_text(
        "".join(f">{r['name']}\n{seqs[r['md5']]}\n" for r in recs))
    one = structures / "colabfold"
    (one / "returned").mkdir(parents=True, exist_ok=True)
    for r in recs:
        (one / f"{r['name']}.fasta").write_text(f">{r['name']}\n{seqs[r['md5']]}\n")
    print(f"deferred to ColabFold: {len(recs)} -> {structures.name}/colabfold/ "
          f"(one FASTA each, plus deferred_for_colabfold.faa)", flush=True)


def deferred_summary_rows(records):
    return [{"status": "deferred_too_long" if "max-len" in r["reason"] else "deferred_too_slow",
             "cluster_id": r["cluster_id"], "name": r["name"], "family": r["family"],
             "anchor_tier": r["anchor_tier"], "amr_classes": r["amr_classes"],
             "product_top": r["product_top"], "n_members": r["n_members"],
             "n_distinct_seqs": r["n_distinct_seqs"], "representative_md5": r["md5"],
             "aa_len": r["aa_len"], "is_fragment": r["is_fragment"], "defer_reason": r["reason"]}
            for r in records]


def write_summary(rows, deferred_records, structures):
    """Merge this run's rows into structure_summary.tsv, keeping clusters the run did not touch."""
    summ = pd.DataFrame(rows + deferred_summary_rows(deferred_records))
    old = structures / "structure_summary.tsv"
    if old.exists():
        prev = read_tsv(old)
        touched = set(summ.cluster_id) if len(summ) else set()
        summ = pd.concat([prev[~prev.cluster_id.isin(touched)], summ]) if len(summ) else prev
    summ = summ.sort_values("cluster_id")
    write_tsv(summ, "structure_summary.tsv", structures)
    return summ


def fmt_cap(seconds):
    return f"{seconds:.0f} s" if seconds < 60 else f"{seconds / 60:.0f} min"


def main():
    t0 = time.time()
    only = set(sys.argv[sys.argv.index("--only") + 1].split(",")) if "--only" in sys.argv else None
    force_cpu, redo = "--cpu" in sys.argv, "--redo" in sys.argv
    max_aa = int(sys.argv[sys.argv.index("--max-len") + 1]) if "--max-len" in sys.argv else MAX_AA
    max_seconds = (float(sys.argv[sys.argv.index("--max-seconds") + 1])
                   if "--max-seconds" in sys.argv else MAX_SECONDS)
    defer_ids = (set(sys.argv[sys.argv.index("--defer") + 1].split(","))
                 if "--defer" in sys.argv else set())
    defer_reason = (sys.argv[sys.argv.index("--defer-reason") + 1]
                    if "--defer-reason" in sys.argv else "named in --defer; fold with ColabFold")
    offload_lm = "--no-offload" not in sys.argv

    d = load_targets()
    seqs = load_sequences(set(d.md5))
    missing = set(d.md5) - set(seqs)
    if missing:
        raise SystemExit(f"{len(missing)} representative sequences not found in proteins_nr.faa")
    d = flag_fragments(d, seqs)
    print(f"clusters: {len(d)}  (fragment representatives: {int(d.is_fragment.sum())}, "
          f"longest {d.aa_len.max()} aa, total {d.aa_len.sum():,} residues)")

    if only:
        d = d[d.cluster_id.isin(only)]
    (STRUCTURES / "pdb").mkdir(parents=True, exist_ok=True)
    (STRUCTURES / "pdb_variation").mkdir(parents=True, exist_ok=True)

    var = None
    if (VARIATION / "position_variation.tsv").exists():
        var = pd.read_csv(VARIATION / "position_variation.tsv", sep="\t",
                          usecols=["cluster_id", "position", "var_frac_lineages"])
        var = {k: g.sort_values("position").var_frac_lineages.to_numpy() * 100
               for k, g in var.groupby("cluster_id")}

    d = d.sort_values("aa_len")                      # short first: fast feedback, memory grows slowly
    d["name"] = d.cluster_id + "__" + d.family.where(d.family != "", d.product_top.fillna("")).map(safe_name)

    # Deferrals decided before folding: too long for this card, or named on the command line.
    d["defer_reason"] = ""
    d.loc[d.aa_len > max_aa, "defer_reason"] = f"longer than --max-len {max_aa}; fold with ColabFold"
    if defer_ids:
        d.loc[d.cluster_id.isin(defer_ids), "defer_reason"] = defer_reason
    deferred_recs = [defer_record(c, c.defer_reason) for c in d[d.defer_reason != ""].itertuples()]
    d = d[d.defer_reason == ""]
    write_deferred(deferred_recs, seqs, STRUCTURES)

    has_pdb = d.name.map(lambda n: (STRUCTURES / "pdb" / f"{n}.pdb").exists())
    todo = list((d if redo else d[~has_pdb]).itertuples())
    print(f"to predict: {len(todo)} ({len(d) - len(todo)} already done); "
          f"time cap {fmt_cap(max_seconds)} per protein", flush=True)

    # A cluster whose PDB is on disk but which has no row in the summary - a run killed between the
    # two - would be skipped for good. Rebuild its row from the files instead of re-folding it.
    rows = []
    if not redo:
        old_summ = STRUCTURES / "structure_summary.tsv"
        known = set(read_tsv(old_summ).cluster_id) if old_summ.exists() else set()
        orphans = [c for c in d[has_pdb].itertuples() if c.cluster_id not in known]
        if orphans:
            rows = [recover_row(c, STRUCTURES) for c in orphans]
            print(f"rebuilt {len(rows)} summary row(s) from structures already on disk "
                  f"(ptm is not stored in a PDB and stays empty for them)", flush=True)

    device = "cpu" if force_cpu else ("cuda" if torch.cuda.is_available() else "cpu")
    tok = AutoTokenizer.from_pretrained("facebook/esmfold_v1")
    model = EsmForProteinFolding.from_pretrained("facebook/esmfold_v1", low_cpu_mem_usage=True).eval()
    if device == "cuda":
        model.esm = model.esm.half()                 # keeps the model at ~8.4 GB instead of filling the card
        torch.backends.cuda.matmul.allow_tf32 = True
    model = model.to(device)
    if offload_lm and device == "cuda":
        install_lm_offload(model)
    print(f"model on {device} ({time.time() - t0:.0f} s"
          f"{'; language model evicted after embedding' if offload_lm and device == 'cuda' else ''})",
          flush=True)

    failed = []
    for i, c in enumerate(todo, 1):
        seq = seqs[c.md5]
        if getattr(model, "_lm_offloaded", False):
            model.esm = model.esm.to(device)      # evicted again by the hook once it has embedded
        plan = [(32, device)] if c.aa_len < LONG_AA else [(16, device), (4, device)]
        if device == "cuda":
            plan.append((8, "cpu"))
        pdb = None

        # A protein that blows the time cap is dropped, and so is every protein still to come: the
        # list is sorted by length, so they are all longer and would only be slower. A CUDA call
        # cannot be interrupted from Python, so the watchdog saves the tables and ends the process.
        def on_timeout(current=c, rest=todo[i:]):
            recs = deferred_recs + [defer_record(x, f"exceeded the {fmt_cap(max_seconds)} "
                                                    f"time cap; fold with ColabFold")
                                    for x in [current] + list(rest)]
            print(f"\n  {current.name} ({current.aa_len} aa) passed the {fmt_cap(max_seconds)} "
                  f"cap - dropping it and the {len(rest)} longer protein(s) still queued", flush=True)
            write_deferred(recs, seqs, STRUCTURES)
            write_summary(rows, recs, STRUCTURES)
            print(f"  wrote {len(rows)} predicted and {len(recs)} deferred rows to "
                  f"{STRUCTURES.name}/structure_summary.tsv", flush=True)
            sys.stdout.flush()
            os._exit(3)

        watchdog = threading.Timer(max_seconds, on_timeout)
        watchdog.daemon = True
        watchdog.start()
        try:
            for chunk, dev in plan:
                try:
                    if dev != device:
                        model_cpu = model.float().to("cpu") if dev == "cpu" else model
                        pdb, plddt, ptm = fold(model_cpu, tok, seq, chunk, dev)
                        model = model.to(device)
                        if device == "cuda":
                            model.esm = model.esm.half()
                    else:
                        pdb, plddt, ptm = fold(model, tok, seq, chunk, dev)
                    used_dev, used_chunk = dev, chunk
                    break
                except torch.cuda.OutOfMemoryError:
                    torch.cuda.empty_cache(); gc.collect()
                    print(f"  {c.name}: out of memory at chunk {chunk} on {dev}, retrying", flush=True)
        finally:
            watchdog.cancel()

        if pdb is None:
            failed.append(c.cluster_id)
            print(f"  {c.name}: FAILED", flush=True)
            continue

        (STRUCTURES / "pdb" / f"{c.name}.pdb").write_text(rewrite_bfactors(pdb, plddt))
        v = var.get(c.cluster_id) if var else None
        if v is None or len(v) != len(seq):
            v = np.zeros(len(seq))
        (STRUCTURES / "pdb_variation" / f"{c.name}.pdb").write_text(rewrite_bfactors(pdb, v))

        rows.append({
            "status": "predicted",
            "cluster_id": c.cluster_id, "name": c.name, "family": c.family,
            "anchor_tier": c.anchor_tier, "amr_classes": c.amr_classes, "product_top": c.product_top,
            "n_members": c.n_members, "n_distinct_seqs": c.n_seqs, "prevalence": c.prevalence,
            "representative_md5": c.md5, "aa_len": c.aa_len,
            "from_step23_reference": isinstance(getattr(c, "reference_md5", ""), str) and bool(getattr(c, "reference_md5", "")),
            "catalog_symbol": getattr(c, "catalog_symbol", ""), "is_fragment": bool(c.is_fragment),
            "only_cluster_of_family": bool(c.only_cluster_of_family),
            "mean_plddt": round(float(plddt.mean()), 2), "min_plddt": round(float(plddt.min()), 2),
            "frac_plddt_ge70": round(float((plddt >= 70).mean()), 4), "ptm": round(ptm, 4),
            "max_variation_pct": round(float(v.max()), 2),
            "device": used_dev, "chunk": used_chunk, "seconds": round(time.time() - t0, 1),
        })
        if device == "cuda":
            torch.cuda.empty_cache()
        if i % 20 == 0 or i == len(todo):
            print(f"  {i}/{len(todo)}  last: {c.name} {c.aa_len} aa "
                  f"pLDDT {rows[-1]['mean_plddt']:.0f}  ({time.time() - t0:.0f} s)", flush=True)

    deferred = pd.DataFrame(deferred_recs)
    summ = write_summary(rows, deferred_recs, STRUCTURES)

    print("\nchecks")
    n_pdb = len(list((STRUCTURES / "pdb").glob("*.pdb")))
    print(f"  structures written: {n_pdb} of {len(d)} foldable clusters; failed: {len(failed)} "
          f"{failed if failed else ''}; deferred to ColabFold: {len(deferred)}")
    s = summ.copy()
    for col in ["mean_plddt", "aa_len", "ptm"]:
        s[col] = pd.to_numeric(s[col])
    print(f"  mean pLDDT: median {s.mean_plddt.median():.1f}; below 70: {(s.mean_plddt < 70).sum()}; "
          f"below 50: {(s.mean_plddt < 50).sum()}")
    print(f"  fragments: {int(s.is_fragment.astype(str).eq('True').sum())}; "
          f"their median pLDDT {s[s.is_fragment.astype(str).eq('True')].mean_plddt.median():.1f} "
          f"vs {s[~s.is_fragment.astype(str).eq('True')].mean_plddt.median():.1f} for the rest")
    # deferred clusters have no PDB, so the length check runs over the predicted ones only
    have = [r for r in s.itertuples() if (STRUCTURES / "pdb" / f"{r.name}.pdb").exists()]
    print(f"  residue count equals sequence length: "
          f"{sum(1 for r in have if count_ca(STRUCTURES / 'pdb' / f'{r.name}.pdb') == int(r.aa_len))}"
          f"/{len(have)}")
    print(f"\ndone in {time.time() - t0:.0f} s")


def count_ca(path):
    return sum(1 for line in open(path) if line.startswith("ATOM") and line[12:16].strip() == "CA")


if __name__ == "__main__":
    main()
