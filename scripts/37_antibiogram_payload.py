"""37 - the antibiogram payload: measured susceptibility beside the genotype of the same isolate.

NCBI holds an antibiogram for a small share of the group (75 of 7,287 isolates for C. freundii), so
this layer is tiny and is kept on its own page rather than folded into the islands or clusters
pages. For every tested isolate it carries:

  tests   one row per antibiotic: the S/I/R call, the measured value with its sign and units, the
          method and the standard
  genes   the isolate's resistance genotype from AMR_annotation/isolate_gene_long.tsv (AMR
          elements, point mutations included), with the analysed cluster of the family where one
          exists, so a row can link into the clusters page
  islands the isolate's tier-1 and tier-2 islands, so a row can link into the islands page, and
          their packed gene maps (a slice of genemaps.json, re-indexed to a local cluster list and
          decoded back as a check) so the page can draw them

and, once, the map from each antibiotic to the AMRFinder subclasses whose genes act on it. That map
is hand-written below and is the only judgement in this step: it decides which genotype rows the
page flags as "expected" for an antibiotic. It errs broad for the beta-lactams (any beta-lactamase
counts for a penicillin) and narrow elsewhere.

Counts are per isolate. The page shows the number of SNP lineages beside every count, because one
outbreak clone tested many times would otherwise dominate an association.

Usage: python3 scripts/37_antibiogram_payload.py [--group G]
"""
import importlib.util
import os
import time

import pandas as pd

import json

from amr_common import DASHBOARD, DATA, ISLANDS, METADATA_TSV, OUT, read_tsv

ANTIBIOGRAMS = DATA / "Antibiograms"           # written by 09_fetch_antibiograms.py

_spec = importlib.util.spec_from_file_location(
    "dash", os.path.join(os.path.dirname(os.path.abspath(__file__)), "30_dashboard_common.py"))
dash = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(dash)

KEEP_TIERS = ("tier1_multi_acquired", "tier2_single_acquired")
TIER3 = "tier3_core_only"
PHENO = {"susceptible": 0, "intermediate": 1, "resistant": 2}      # "not defined" is dropped

# antibiotic -> AMRFinder subclass tokens (as in isolate_gene_long.subclass, split on ";" and "/")
# whose genes act on it. Used only to flag which genotype rows are EXPECTED to matter.
BETA_LACTAMASES = ["BETA-LACTAM", "CEPHALOSPORIN", "CARBAPENEM"]
ANTIBIOTIC_TARGETS = {
    "ampicillin": BETA_LACTAMASES, "amoxicillin-clavulanic acid": BETA_LACTAMASES,
    "ampicillin-sulbactam": BETA_LACTAMASES, "piperacillin-tazobactam": BETA_LACTAMASES,
    "ticarcillin-clavulanic acid": BETA_LACTAMASES,
    "cefazolin": ["CEPHALOSPORIN", "CARBAPENEM"], "cefoxitin": ["CEPHALOSPORIN", "CARBAPENEM"],
    "cefotaxime": ["CEPHALOSPORIN", "CARBAPENEM"], "ceftriaxone": ["CEPHALOSPORIN", "CARBAPENEM"],
    "ceftazidime": ["CEPHALOSPORIN", "CARBAPENEM"], "cefepime": ["CEPHALOSPORIN", "CARBAPENEM"],
    "cefuroxime": ["CEPHALOSPORIN", "CARBAPENEM"], "aztreonam": ["CEPHALOSPORIN", "CARBAPENEM"],
    "meropenem": ["CARBAPENEM"], "ertapenem": ["CARBAPENEM"], "imipenem": ["CARBAPENEM"],
    "doripenem": ["CARBAPENEM"],
    "amikacin": ["AMIKACIN", "AMINOGLYCOSIDE"], "gentamicin": ["GENTAMICIN", "AMINOGLYCOSIDE"],
    "tobramycin": ["TOBRAMYCIN", "AMINOGLYCOSIDE"], "kanamycin": ["KANAMYCIN", "AMINOGLYCOSIDE"],
    "streptomycin": ["STREPTOMYCIN"],
    "ciprofloxacin": ["QUINOLONE"], "levofloxacin": ["QUINOLONE"], "nalidixic acid": ["QUINOLONE"],
    "trimethoprim-sulfamethoxazole": ["TRIMETHOPRIM", "SULFONAMIDE"], "trimethoprim": ["TRIMETHOPRIM"],
    "tetracycline": ["TETRACYCLINE"], "doxycycline": ["TETRACYCLINE"], "minocycline": ["TETRACYCLINE"],
    "tigecycline": ["TIGECYCLINE"],
    "nitrofurantoin": ["NITROFURANTOIN", "NITROFURAN"], "colistin": ["COLISTIN"],
    "chloramphenicol": ["CHLORAMPHENICOL", "PHENICOL"], "fosfomycin": ["FOSFOMYCIN"],
    "erythromycin": ["ERYTHROMYCIN", "MACROLIDE"], "azithromycin": ["AZITHROMYCIN", "MACROLIDE"],
    "linezolid": ["LINEZOLID", "OXAZOLIDINONE"], "vancomycin": ["VANCOMYCIN", "GLYCOPEPTIDE"],
    "rifampin": ["RIFAMYCIN"], "clindamycin": ["CLINDAMYCIN", "LINCOSAMIDE"],
    "penicillin": BETA_LACTAMASES, "amoxicillin": BETA_LACTAMASES, "piperacillin": BETA_LACTAMASES,
    "ticarcillin": BETA_LACTAMASES, "carbenicillin": BETA_LACTAMASES, "oxacillin": BETA_LACTAMASES,
    "cephalothin": ["CEPHALOSPORIN", "CARBAPENEM"], "cefotetan": ["CEPHALOSPORIN", "CARBAPENEM"],
    "cefpodoxime": ["CEPHALOSPORIN", "CARBAPENEM"], "ceftaroline": ["CEPHALOSPORIN", "CARBAPENEM"],
    "ceftazidime-avibactam": ["CEPHALOSPORIN", "CARBAPENEM"], "ceftolozane-tazobactam": ["CEPHALOSPORIN", "CARBAPENEM"],
    "moxifloxacin": ["QUINOLONE"], "sulfamethoxazole": ["SULFONAMIDE"], "sulfisoxazole": ["SULFONAMIDE"],
    "plazomicin": ["AMINOGLYCOSIDE"], "daptomycin": [],        # daptomycin: no resistance gene in this genus
}


def tok(s):
    return [t for t in str(s or "").replace("/", ";").split(";") if t and t != "nan"]


def num(v, scale=1, default=-1):
    try:
        f = float(v)
    except (TypeError, ValueError):
        return default
    return default if f != f else int(round(f * scale))


def slice_genemaps(island_nums):
    """The packed gene maps of `island_nums`, with cluster indices remapped to a local list.

    The stream is decoded, its cluster references rewritten, and re-packed with the same codec;
    then decoded again and compared gene by gene with the original, so a slip here cannot draw a
    gene under the wrong name."""
    src = DASHBOARD / "genemaps.json"
    empty = {"bp": 10, "roles": [], "flag_label": 16, "flag_pseudo": 32, "cluster": [], "clabel": [],
             "cclass": [], "dicts": {"label": [], "amrclass": []}, "g": {}}
    if not src.exists() or not island_nums:
        return empty, 0
    gm = json.loads(src.read_text(encoding="utf-8"))
    fl = gm["flag_label"]
    local, out, n_genes = {}, {}, 0

    def decode(packed):
        v = dash.varints(packed)
        genes, i = [], 1
        while i < len(v):
            gap, w, cs, flag = v[i], v[i + 1], v[i + 2], v[i + 3]
            i += 4
            lab = None
            if flag & fl:
                lab = v[i]
                i += 1
            genes.append((gap, w, cs, flag, lab))
        return v[0], genes

    for n in sorted(island_nums):
        packed = gm["g"].get(str(n))
        if packed is None:
            continue
        first, genes = decode(packed)
        parts = [dash.vint(first)]
        for gap, w, cs, flag, lab in genes:
            ci = (cs >> 1) - 1
            if ci >= 0:
                if ci not in local:
                    local[ci] = len(local)
                cs = (local[ci] + 1) * 2 + (cs & 1)
            parts += [dash.vint(gap), dash.vint(w), dash.vint(cs), dash.vint(flag)]
            if lab is not None:
                parts.append(dash.vint(lab))
        out[n] = "".join(parts)
        n_genes += len(genes)
        # the check: same positions, flags, labels and (through the local list) the same clusters
        f2, g2 = decode(out[n])
        assert f2 == first and len(g2) == len(genes)
        for (gap, w, cs, flag, lab), (gap2, w2, cs2, flag2, lab2) in zip(genes, g2):
            ci, ci2 = (cs >> 1) - 1, (cs2 >> 1) - 1
            assert (gap, w, flag, lab, cs & 1) == (gap2, w2, flag2, lab2, cs2 & 1)
            assert (ci < 0 and ci2 < 0) or gm["cluster"][ci] == gm["cluster"][list(local)[ci2]]
    order = sorted(local, key=local.get)
    return {"bp": gm["bp"], "roles": gm["roles"], "flag_label": fl, "flag_pseudo": gm["flag_pseudo"],
            "cluster": [gm["cluster"][ci] for ci in order], "clabel": [gm["clabel"][ci] for ci in order],
            "cclass": [gm["cclass"][ci] for ci in order],
            "dicts": {"label": gm["dicts"]["label"], "amrclass": gm["dicts"]["amrclass"]},
            "g": out}, n_genes


def main():
    t0 = time.time()
    full = ANTIBIOGRAMS / "antibiogram_full.tsv"
    d = dash.Dicts()
    group_isolates = sum(1 for _ in open(METADATA_TSV, encoding="utf-8-sig")) - 1
    if not full.exists():
        dash.write_payload({"built": time.strftime("%Y-%m-%d"), "absent": 1, "group_isolates": group_isolates,
                            "dicts": {}, "targets": {}, "isolates": []}, "antibiograms.json")
        print(f"no antibiogram table at {full}; wrote an empty payload (run 09_fetch_antibiograms.py)")
        return

    ab = read_tsv(full)
    ab = ab[ab.resistance_phenotype.isin(PHENO)]
    tested = sorted(ab.target_acc.unique())
    print(f"{len(ab):,} defined tests on {len(tested)} isolates ({group_isolates:,} in the group)", flush=True)

    # islands of the tested isolates
    isl = read_tsv(ISLANDS / "islands.tsv", usecols=["island_id", "target_acc", "island_tier"])
    isl = isl[isl.target_acc.isin(tested)]
    islands_of = {t: sorted(dash.island_num(i) for i in grp.island_id) for t, grp in isl[isl.island_tier.isin(KEEP_TIERS)].groupby("target_acc")}
    islands3_of = {t: sorted(dash.island_num(i) for i in grp.island_id) for t, grp in isl[isl.island_tier == TIER3].groupby("target_acc")}

    # family -> the analysed cluster with the most members that lists it
    anchors = read_tsv(ISLANDS / "anchor_clusters.tsv", usecols=["cluster_id", "n_members", "amr_gene_families"])
    fam_cluster = {}
    for r in anchors.itertuples():
        for fam in tok(r.amr_gene_families):
            n = num(r.n_members, default=0)
            if fam not in fam_cluster or n > fam_cluster[fam][1]:
                fam_cluster[fam] = (dash.cluster_num(r.cluster_id), n)

    # the genotype of the tested isolates
    genes = read_tsv(OUT / "isolate_gene_long.tsv",
                     usecols=["target_acc", "element_type", "symbol", "call_qualifier", "gene_family",
                              "class", "subclass", "subclass_all"])
    genes = genes[(genes.element_type == "AMR") & genes.target_acc.isin(tested)]
    genes = genes.drop_duplicates(["target_acc", "symbol"]).rename(columns={"class": "cls"})

    meta_cols = ["target_acc", "biosample_acc", "asm_acc", "scientific_name", "PDS_acc"]
    meta = ab[meta_cols].drop_duplicates("target_acc").set_index("target_acc")
    has_genome = {}
    isolates_tsv = ANTIBIOGRAMS / "antibiogram_isolates.tsv"
    if isolates_tsv.exists():
        it = read_tsv(isolates_tsv, usecols=["target_acc", "has_genome"])
        has_genome = dict(zip(it.target_acc, it.has_genome))

    out = []
    for t in tested:
        m = meta.loc[t]
        tests = []
        for r in ab[ab.target_acc == t].sort_values("antibiotic").itertuples():
            tests.append([d.enc("antibiotic", r.antibiotic), PHENO[r.resistance_phenotype],
                          str(r.measurement_sign or ""), str(r.measurement or ""),
                          d.enc("units", r.measurement_units), d.enc("method", r.laboratory_typing_method),
                          d.enc("standard", r.testing_standard)])
        gl = []
        for r in genes[genes.target_acc == t].sort_values(["cls", "symbol"]).itertuples():
            fam = str(r.gene_family or r.symbol)
            gl.append([d.enc("symbol", r.symbol), d.enc("family", fam), d.enc("class", r.cls),
                       d.enc("subclass", r.subclass), 1 if str(r.call_qualifier).upper().startswith("POINT") else 0,
                       fam_cluster.get(fam, (-1, 0))[0], 0 if str(r.call_qualifier) in ("complete", "exact") else 1])
        out.append({"t": t, "bs": m.biosample_acc, "asm": m.asm_acc, "sp": d.enc("species", m.scientific_name),
                    "pds": d.enc("pds", m.PDS_acc), "genome": 1 if str(has_genome.get(t, "1")) in ("1", "True") else 0,
                    "islands": islands_of.get(t, []), "islands3": islands3_of.get(t, []), "tests": tests, "genes": gl})

    # the gene maps of those islands: a slice of genemaps.json (step 33 runs first), with the
    # cluster references re-indexed to a local list so the page carries only what it draws
    gm_slice, n_map_genes = slice_genemaps({n for i in out for n in i["islands"] + i["islands3"]})

    # antibiotic -> subclass tokens, only for antibiotics actually in the data
    targets = {}
    for a in d.as_dict().get("antibiotic", []):
        targets[d.enc("antibiotic", a)] = ANTIBIOTIC_TARGETS.get(a, [])
    unmapped = [a for a in d.as_dict().get("antibiotic", []) if a not in ANTIBIOTIC_TARGETS]

    payload = {"built": time.strftime("%Y-%m-%d"), "absent": 0, "group_isolates": group_isolates,
               "dicts": d.as_dict(), "targets": targets, "isolates": out, "gm": gm_slice}
    mb = dash.write_payload(payload, "antibiograms.json")

    n_tests = sum(len(i["tests"]) for i in out)
    ph = [0, 0, 0]
    for i in out:
        for x in i["tests"]:
            ph[x[1]] += 1
    print("\nchecks")
    print(f"  isolates: {len(out)}; with tier-1/2 islands: {sum(1 for i in out if i['islands'])}; "
          f"with tier-3 windows: {sum(1 for i in out if i['islands3'])}; "
          f"with a genome: {sum(i['genome'] for i in out)}")
    print(f"  tests: {n_tests:,} (S {ph[0]}, I {ph[1]}, R {ph[2]}); antibiotics: {d.size('antibiotic')}")
    print(f"  genotype rows: {sum(len(i['genes']) for i in out):,}; families: {d.size('family')}; "
          f"linked to an analysed cluster: {sum(1 for i in out for g in i['genes'] if g[5] >= 0):,}")
    print(f"  antibiotics without a subclass map: {unmapped if unmapped else 'none'}")
    print(f"  gene maps carried: {len(gm_slice['g'])} islands, {n_map_genes:,} genes, "
          f"{len(gm_slice['cluster']):,} clusters referenced (decoded back and matched)")
    print(f"  size {mb:.2f} MB")
    print(f"\ndone in {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
