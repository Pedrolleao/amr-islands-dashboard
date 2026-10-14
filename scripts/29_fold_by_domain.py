"""29 - fold a protein that is too long to fold whole, one module at a time.

For a modular protein - an NRPS/PKS above all - a single-chain prediction is both the hardest thing to
compute and the least trustworthy to interpret, because the arrangement of the modules is genuinely
flexible. The confident part of any whole-chain model would be the individual domains anyway.

So: take the Pfam architecture from step 27, group consecutive domains into modules no longer than
--max-len, cut in the linkers between them, and fold each module with step 24's protocol. The pieces are
written as ordinary structures with their own rows, so steps 25 and 28 pick them up unchanged, and the
parent cluster keeps a row of its own marked `folded_by_domain` so it can never read as a whole-chain
model and can never read as missing.

Boundaries: a cut falls at the START of the first domain of the next module, so the linker stays with the
module before it. Each module is then extended by --overlap residues into the next, so the interface is
present in both pieces and no residue is lost at a join.

Usage: python3 scripts/29_fold_by_domain.py --cluster CFR-014064 [--group G] [--set virulence]
                                            [--max-len 1200] [--overlap 20] [--max-seconds 3600]
"""
import gc
import importlib.util
import os
import sys
import threading
import time

import numpy as np
import pandas as pd

from amr_common import DOMAINS, PANGENOME, STRUCTURES, VARIATION, read_tsv, set_dir, write_tsv
from cluster_sets import load_cluster_set

DOMAINS, STRUCTURES, VARIATION = set_dir(DOMAINS), set_dir(STRUCTURES), set_dir(VARIATION)


def step24():
    """Reuse step 24's folding helpers rather than copying them (its name starts with a digit)."""
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "24_structure_prediction.py")
    spec = importlib.util.spec_from_file_location("step24", path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def modules_from_domains(dom, length, max_len):
    """Group consecutive Pfam domains into modules no longer than max_len; cut before a domain."""
    cuts, start = [1], 1
    for r in dom.sort_values("start").itertuples():
        if r.end - start + 1 > max_len:
            cuts.append(int(r.start))
            start = int(r.start)
    cuts.append(length + 1)
    out = []
    for i in range(len(cuts) - 1):
        a, b = cuts[i], cuts[i + 1] - 1
        names = dom[(dom.start >= a) & (dom.start <= b)].sort_values("start").pfam_name
        out.append({"module": i + 1, "start": a, "end": b, "domains": " - ".join(names)})
    return out


def main():
    t0 = time.time()
    cid = sys.argv[sys.argv.index("--cluster") + 1] if "--cluster" in sys.argv else None
    if not cid:
        raise SystemExit("--cluster is required, e.g. --cluster CFR-014064")
    max_len = int(sys.argv[sys.argv.index("--max-len") + 1]) if "--max-len" in sys.argv else 1200
    overlap = int(sys.argv[sys.argv.index("--overlap") + 1]) if "--overlap" in sys.argv else 20
    cap = float(sys.argv[sys.argv.index("--max-seconds") + 1]) if "--max-seconds" in sys.argv else 3600

    s24 = step24()
    d = load_cluster_set()
    row = d[d.cluster_id == cid]
    if row.empty:
        raise SystemExit(f"{cid} is not in this cluster set")
    row = row.iloc[0]
    seqs = s24.load_sequences({row.rep_md5})
    seq = seqs[row.rep_md5]
    gene = (row.amr_gene_families or row.product_top or "").split(";")[0]
    base = f"{cid}__{s24.safe_name(gene)}"
    print(f"{base}: {len(seq)} aa", flush=True)

    dom = read_tsv(DOMAINS / "pfam_domains.tsv")
    dom = dom[dom.cluster_id == cid]
    if dom.empty:
        raise SystemExit(f"no Pfam domains for {cid} - run step 27 for this set first")
    dom["start"] = dom.start.astype(int); dom["end"] = dom.end.astype(int)
    mods = modules_from_domains(dom, len(seq), max_len)
    print(f"{len(dom)} domains -> {len(mods)} module(s), max {max_len} aa, {overlap} aa overlap:")
    for m in mods:
        hi = min(len(seq), m["end"] + overlap)
        m["fold_start"], m["fold_end"] = m["start"], hi
        print(f"  m{m['module']}: {m['start']}-{hi} ({hi - m['start'] + 1} aa)  {m['domains']}",
              flush=True)

    var = None
    vf = VARIATION / "position_variation.tsv"
    if vf.exists():
        v = pd.read_csv(vf, sep="\t", usecols=["cluster_id", "position", "var_frac_lineages"])
        v = v[v.cluster_id == cid].sort_values("position")
        if len(v):
            var = v.var_frac_lineages.to_numpy() * 100
            print(f"per-position variation available ({len(var)} positions)", flush=True)

    tok, model, device = s24.AutoTokenizer.from_pretrained("facebook/esmfold_v1"), None, "cuda"
    model = s24.EsmForProteinFolding.from_pretrained("facebook/esmfold_v1",
                                                     low_cpu_mem_usage=True).eval()
    model.esm = model.esm.half()
    s24.torch.backends.cuda.matmul.allow_tf32 = True
    model = model.to(device)
    s24.install_lm_offload(model)
    print(f"model on {device}, language model evicted after embedding "
          f"({time.time() - t0:.0f} s)", flush=True)

    (STRUCTURES / "pdb").mkdir(parents=True, exist_ok=True)
    (STRUCTURES / "pdb_variation").mkdir(parents=True, exist_ok=True)
    rows = []
    for m in mods:
        a, b = m["fold_start"], m["fold_end"]
        sub = seq[a - 1:b]
        name = f"{base}_m{m['module']}_{a}-{b}"
        model.esm = model.esm.to(device)
        s24.torch.cuda.empty_cache(); s24.torch.cuda.reset_peak_memory_stats()
        timer = threading.Timer(cap, lambda n=name: (print(f"  {n}: past the cap, stopping",
                                                           flush=True), os._exit(3)))
        timer.daemon = True; timer.start()
        t1 = time.time()
        pdb = plddt = None
        for chunk in (16, 4):
            try:
                pdb, plddt, ptm = s24.fold(model, tok, sub, chunk, device)
                break
            except s24.torch.cuda.OutOfMemoryError:
                s24.torch.cuda.empty_cache(); gc.collect()
                print(f"  {name}: out of memory at chunk {chunk}, retrying", flush=True)
        timer.cancel()
        if pdb is None:
            print(f"  {name}: FAILED", flush=True)
            continue
        (STRUCTURES / "pdb" / f"{name}.pdb").write_text(s24.rewrite_bfactors(pdb, plddt))
        v = var[a - 1:b] if var is not None and len(var) >= b else np.zeros(len(sub))
        (STRUCTURES / "pdb_variation" / f"{name}.pdb").write_text(s24.rewrite_bfactors(pdb, v))
        print(f"  {name}: {time.time() - t1:.0f} s, mean pLDDT {plddt.mean():.2f}, "
              f"peak {s24.torch.cuda.max_memory_allocated() / 1e9:.2f} GB", flush=True)
        rows.append({
            "status": "predicted", "cluster_id": f"{cid}.m{m['module']}", "name": name,
            "family": gene, "anchor_tier": row.anchor_tier, "amr_classes": row.amr_classes,
            "product_top": row.product_top, "n_members": row.n_members,
            "n_distinct_seqs": row.n_seqs, "prevalence": row.prevalence,
            "representative_md5": row.rep_md5, "aa_len": len(sub), "is_fragment": False,
            "mean_plddt": round(float(plddt.mean()), 2), "min_plddt": round(float(plddt.min()), 2),
            "frac_plddt_ge70": round(float((plddt >= 70).mean()), 4), "ptm": round(ptm, 4),
            "max_variation_pct": round(float(np.max(v)), 2), "device": "cuda",
            "source": "esmfold_by_domain", "parent_cluster_id": cid,
            "module": m["module"], "module_start": a, "module_end": b,
            "module_domains": m["domains"],
        })

    if not rows:
        raise SystemExit("no module folded")
    w = np.array([r["aa_len"] for r in rows], float)
    parent = {
        "status": "folded_by_domain", "cluster_id": cid, "name": base, "family": gene,
        "anchor_tier": row.anchor_tier, "amr_classes": row.amr_classes,
        "product_top": row.product_top, "n_members": row.n_members,
        "n_distinct_seqs": row.n_seqs, "representative_md5": row.rep_md5, "aa_len": len(seq),
        "is_fragment": False, "source": "esmfold_by_domain",
        "mean_plddt": round(float(np.average([r["mean_plddt"] for r in rows], weights=w)), 2),
        "defer_reason": "",
        "module_domains": "; ".join(f"m{r['module']} {r['module_start']}-{r['module_end']}"
                                    for r in rows),
    }
    summ = pd.DataFrame(rows + [parent])
    old = STRUCTURES / "structure_summary.tsv"
    if old.exists():
        prev = read_tsv(old)
        summ = pd.concat([prev[~prev.cluster_id.isin(set(summ.cluster_id))], summ])
    write_tsv(summ.sort_values("cluster_id"), "structure_summary.tsv", STRUCTURES)

    print("\nchecks")
    covered = set()
    for r in rows:
        covered |= set(range(r["module_start"], r["module_end"] + 1))
    print(f"  modules folded: {len(rows)}; residues covered {len(covered)}/{len(seq)} "
          f"({'complete' if len(covered) == len(seq) else 'GAP'})")
    print(f"  length-weighted mean pLDDT {parent['mean_plddt']}; "
          f"per module " + ", ".join(f"m{r['module']} {r['mean_plddt']:.1f}" for r in rows))
    print(f"  parent row status: {parent['status']} (no whole-chain PDB exists, by design)")
    print(f"\ndone in {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
