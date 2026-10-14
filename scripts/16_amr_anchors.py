"""16 - decide which protein clusters anchor an island, and what KIND of anchor each one is.

THE RULE, as specified: a cluster is AMR-associated if at least one of its members is a known AMR
hit, and then EVERY member of that cluster anchors an island, in every genome. That is deliberate -
it propagates the annotation to homologues the annotation missed, which is the main reason to
cluster at all. This script applies that rule literally and never removes a cluster from it.

WHAT IT ADDS, AND WHY IT IS NOT A FILTER

Cluster-level anchoring has a consequence that has to be measured rather than assumed. `gyrA` carries
the ~3,200 point-mutation hits, but the gene itself is in every genome, so the gyrA cluster anchors
one island per genome - about 6,000 stretches of ordinary chromosome. Chromosomal `ampC`/`blaCMY`
behaves the same way. Those are not wrong, they are simply not mobile resistance islands, and mixing
them with a Tn-borne blaKPC cassette would drown the signal.

So every anchor cluster gets an `anchor_tier`, evaluated in order:

    POINT_ONLY     only point-mutation evidence          -> gyrA: an essential gene with a
                                                            resistance substitution, not cargo
    CORE_AMR       AMR evidence, prevalence >= CORE      -> chromosomal ampC, efflux pumps
    SOFTCORE_AMR   AMR evidence, prevalence >= SOFTCORE
    SHELL_AMR      AMR evidence, prevalence >= SHELL
    MOBILE_AMR     AMR evidence, prevalence <  SHELL     -> blaKPC-2, blaNDM-1, aadA1: the real
                                                            resistance islands

PREVALENCE IS TAKEN AS max(genus-wide, best within a single species). The group is a whole genus:
*C. koseri* is ~10% of the genomes, so a gene core to *C. koseri* sits at genus prevalence 0.10 and
would be called mobile by a genus-wide threshold alone. The species-aware value catches that - it is
what correctly labels the species-specific chromosomal beta-lactamases (blaSED in *C. sedlakii*,
blaCKO in *C. koseri*, blaTEL in *C. farmeri*) and chromosomal qnrB and oqxAB as backbone.

BUT SPECIES PREVALENCE ALONE CAN BE FAKED BY AN OUTBREAK. *C. sedlakii* has 64 genomes spread over
only 14 SNP clusters, so one clonal expansion carrying a plasmid can put a gene in half the species
and make a mobile carbapenemase look chromosomal. A species prevalence therefore only counts when
the cluster is seen in at least MIN_SPECIES_PDS distinct SNP clusters (PDS_acc) of that species -
i.e. it has spread across lineages rather than been inherited once. Same discipline as counting
recurrence in lineages rather than genomes.

CALIBRATION, NOT A GUESSED THRESHOLD. `anchor_prevalence_sweep.tsv` scores candidate thresholds
against a named list of genes that must come out core and must come out mobile (in `groups.json`,
beside the existing `spot_checks`). `sul1` is the hard case: genuinely acquired, yet present in
about half the genomes. The sweep is printed so the threshold is chosen from the data.

Usage: python3 scripts/16_amr_anchors.py [--group G]
"""
import json
import time

import pandas as pd

from amr_common import GROUPS, GROUP, ISLANDS, PANGENOME, read_tsv, write_tsv

CORE_PREV = 0.95          # defaults; the sweep below says whether they are the right ones
SOFTCORE_PREV = 0.90
SHELL_PREV = 0.15
MIN_SPECIES_N = 50        # a species needs this many genomes before its own prevalence counts
MIN_SPECIES_PDS = 3       # ... and the cluster must span this many SNP clusters within it
SWEEP = [0.80, 0.85, 0.90, 0.95, 0.98, 0.99]

# Genes whose tier is not a matter of opinion, matched on the AMR symbols a cluster actually
# carries. Only genes AMRFinderPlus reports can appear here - acrAB/tolC/marA are not in its
# catalogue, so listing them would only ever produce "not found". Overridable in groups.json.
DEFAULT_CALIBRATION = {
    # intrinsic / chromosomal in Citrobacter: must NOT be treated as mobile cargo
    "must_be_core": ["blaCMY", "qnrB", "blaSED", "blaCKO", "blaTEL", "oqxA", "oqxB", "gyrA_T83I"],
    # acquired, plasmid- and transposon-borne: must NOT be treated as backbone
    # NOTE: qnrB2 is deliberately NOT listed. Measured here: at 90% identity the plasmid-borne
    # qnrB alleles (qnrB2, qnrB19 ...) cluster together with the CHROMOSOMAL qnrB that is intrinsic
    # to C. werkmanii and C. freundii - they are >90% identical, so no clustering at this threshold
    # can separate them. The cluster is therefore correctly CORE_AMR for the chromosomal copies and
    # simultaneously contains acquired ones. Islands anchored by it carry `anchor_mixed_qnr` in the
    # QC report; separating them needs allele-level assignment, not clustering.
    "must_be_mobile": ["blaKPC-2", "blaNDM-1", "blaOXA-48", "blaSHV-12", "catA2",
                       "aadA1", "tet(A)", "mph(A)", "dfrA19", "sul1", "sul2", "blaTEM-1"],
}
CORE_TIERS = {"CORE_AMR", "SOFTCORE_AMR", "POINT_ONLY"}
CARGO_TIERS = {"SHELL_AMR", "MOBILE_AMR"}


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def main():
    cfg = GROUPS[GROUP]
    calib = cfg.get("anchor_calibration", DEFAULT_CALIBRATION)
    core_prev = float(cfg.get("core_prevalence", CORE_PREV))
    soft_prev = float(cfg.get("softcore_prevalence", SOFTCORE_PREV))
    shell_prev = float(cfg.get("shell_prevalence", SHELL_PREV))

    ann = read_tsv(PANGENOME / "cluster_annotation.tsv")
    for c in ("n_members", "n_genomes", "n_amr_hits", "n_point_hits", "n_stress_hits",
              "n_virulence_hits", "n_amr_hits_strong"):
        ann[c] = pd.to_numeric(ann[c])
    ann["prevalence"] = pd.to_numeric(ann["prevalence"])
    log(f"{len(ann):,} clusters")

    # ---- per-species prevalence ------------------------------------------------------------------
    gi = read_tsv(PANGENOME / "genome_index.tsv")
    man = read_tsv(ISLANDS.parent / "Genomes" / "manifest" / "genome_manifest.tsv")
    species_of = dict(zip(man["download_acc"], man["scientific_name"].str.split().str[:2].str.join(" ")))
    gi["species"] = gi["download_acc"].map(species_of).fillna("")
    species_n = gi["species"].value_counts()
    big = set(species_n[species_n >= MIN_SPECIES_N].index)
    log(f"{len(big)} species with >= {MIN_SPECIES_N} genomes: "
        f"{', '.join(f'{s} ({species_n[s]:,})' for s in sorted(big, key=lambda s: -species_n[s])[:5])}")

    # Only clusters that carry AMR or point-mutation evidence can ever be anchors, so the
    # species-aware prevalence is computed for those alone (a few hundred, not 208,757).
    evidence = set(ann.loc[(ann["n_amr_hits"] > 0) | (ann["n_point_hits"] > 0), "cluster_id"])
    log(f"{len(evidence):,} clusters carry AMR or point evidence")

    import gzip
    sp_of_acc = dict(zip(gi["download_acc"], gi["species"]))
    pds_of_acc = dict(zip(man["download_acc"], man["PDS_acc"]))
    species_pds = {}                      # species -> set of PDS present in the dataset
    for acc, sp in sp_of_acc.items():
        pds = pds_of_acc.get(acc, "")
        if sp and pds:
            species_pds.setdefault(sp, set()).add(pds)

    per_species = {}                      # cluster -> species -> genome count
    per_species_pds = {}                  # cluster -> species -> set of PDS
    cur, local = None, set()

    def flush(acc):
        sp = sp_of_acc.get(acc, "")
        if sp not in big:
            return
        pds = pds_of_acc.get(acc, "")
        for cid in local:
            per_species.setdefault(cid, {})
            per_species[cid][sp] = per_species[cid].get(sp, 0) + 1
            if pds:
                per_species_pds.setdefault(cid, {}).setdefault(sp, set()).add(pds)

    t0 = time.time()
    with gzip.open(PANGENOME / "gene_cluster_map.tsv.gz", "rt") as fh:
        for line in fh:
            gene_id, cid = line.rstrip("\n").split("\t")
            acc = gene_id.split("|", 1)[0]
            if acc != cur:
                if cur is not None:
                    flush(cur)
                cur, local = acc, set()
            if cid and cid in evidence:
                local.add(cid)
        if cur is not None:
            flush(cur)
    log(f"per-species prevalence computed in {time.time()-t0:.0f}s")

    best_sp, best_sp_name, best_sp_pds = {}, {}, {}
    for cid, d in per_species.items():
        options = []
        for sp, n in d.items():
            n_pds = len(per_species_pds.get(cid, {}).get(sp, ()))
            # a species prevalence only counts if the cluster spans several lineages of that
            # species - otherwise one outbreak would make a plasmid gene look chromosomal
            if n_pds >= MIN_SPECIES_PDS or len(species_pds.get(sp, ())) < MIN_SPECIES_PDS:
                options.append((n / species_n[sp], sp, n_pds))
        if options:
            frac, sp, n_pds = max(options)
            best_sp[cid], best_sp_name[cid], best_sp_pds[cid] = frac, sp, n_pds

    ann["prevalence_species_max"] = ann["cluster_id"].map(best_sp).fillna(0).round(6)
    ann["prevalence_species"] = ann["cluster_id"].map(best_sp_name).fillna("")
    ann["prevalence_species_pds"] = ann["cluster_id"].map(best_sp_pds).fillna(0).astype(int)
    ann["prevalence_eff"] = ann[["prevalence", "prevalence_species_max"]].max(axis=1)

    # ---- anchor tier -----------------------------------------------------------------------------
    def tier(r, core=core_prev, soft=soft_prev, shell=shell_prev):
        if r.n_amr_hits == 0 and r.n_point_hits > 0:
            return "POINT_ONLY"
        if r.n_amr_hits == 0:
            return ("STRESS_ONLY" if r.n_stress_hits else
                    "VIRULENCE_ONLY" if r.n_virulence_hits else "BACKBONE")
        p = r.prevalence_eff
        return ("CORE_AMR" if p >= core else "SOFTCORE_AMR" if p >= soft
                else "SHELL_AMR" if p >= shell else "MOBILE_AMR")

    ann["anchor_tier"] = [tier(r) for r in ann.itertuples()]
    ann["is_anchor"] = ann["anchor_tier"].isin(
        ["POINT_ONLY", "CORE_AMR", "SOFTCORE_AMR", "SHELL_AMR", "MOBILE_AMR"]).astype(int)
    ann["is_mobile_anchor"] = ann["anchor_tier"].isin(["MOBILE_AMR", "SHELL_AMR"]).astype(int)
    ann["backbone_risk"] = ann["prevalence_eff"].ge(soft_prev).astype(int)
    ann["expected_islands"] = ann["n_members"] * ann["is_anchor"]
    write_tsv(ann, "cluster_annotation.tsv", folder=PANGENOME)

    anchors = ann[ann["is_anchor"] == 1]
    write_tsv(anchors[["cluster_id", "anchor_tier", "prevalence", "prevalence_species_max",
                       "prevalence_species", "prevalence_eff", "n_members", "n_genomes",
                       "n_amr_hits", "n_point_hits", "evidence_tier", "propagation_ratio",
                       "amr_symbols", "amr_gene_families", "amr_classes", "product_top"]],
              "anchor_clusters.tsv", folder=ISLANDS)

    # ---- threshold sweep --------------------------------------------------------------------------
    # Scored on the symbols a cluster actually carries. A calibration gene that appears in NO
    # cluster is reported as "not found", never as a failure: absence says the catalogue does not
    # report that gene here, which is not evidence about the threshold.
    must_core = [s for s in calib["must_be_core"]]
    must_mobile = [s for s in calib["must_be_mobile"]]

    def symbols_of(r):
        return {x for x in (str(r.amr_symbols).split(";") + str(r.point_symbols).split(";")) if x}

    amr_like = ann[(ann["n_amr_hits"] > 0) | (ann["n_point_hits"] > 0)].copy()
    sym_index = {}
    for r in amr_like.itertuples():
        for sym in symbols_of(r):
            sym_index.setdefault(sym, []).append(r.cluster_id)

    rows = []
    for thr in SWEEP:
        tiers = dict(zip(amr_like["cluster_id"],
                         [tier(r, core=thr, soft=max(thr - 0.05, 0)) for r in amr_like.itertuples()]))
        ok = fail = missing = 0
        failures = []
        for sym, want in [(s, CORE_TIERS) for s in must_core] + \
                         [(s, CARGO_TIERS) for s in must_mobile]:
            cids = sym_index.get(sym)
            if not cids:
                missing += 1
                continue
            got = {tiers[c] for c in cids}
            if got & want:
                ok += 1
            else:
                fail += 1
                failures.append(f"{sym}->{'/'.join(sorted(got))}")
        counts = pd.Series(list(tiers.values())).value_counts()
        rows.append({"core_prevalence": thr,
                     "n_core_amr": int(counts.get("CORE_AMR", 0)),
                     "n_softcore_amr": int(counts.get("SOFTCORE_AMR", 0)),
                     "n_shell_amr": int(counts.get("SHELL_AMR", 0)),
                     "n_mobile_amr": int(counts.get("MOBILE_AMR", 0)),
                     "calibration_ok": ok, "calibration_failed": fail,
                     "calibration_not_found": missing,
                     "failures": ";".join(failures[:6])})
    sweep = pd.DataFrame(rows)
    write_tsv(sweep, "anchor_prevalence_sweep.tsv", folder=ISLANDS)

    # ---- QC -----------------------------------------------------------------------------------------
    print("\n===== QC 16_amr_anchors =====")
    print(f"thresholds in use: core={core_prev} softcore={soft_prev} shell={shell_prev}")
    print("\nanchor tiers:")
    for t, n in ann["anchor_tier"].value_counts().items():
        sub = ann[ann["anchor_tier"] == t]
        print(f"  {t:15s} {n:7,} clusters  {sub['n_members'].sum():10,} CDS")
    print(f"\nanchor clusters: {len(anchors):,}")
    print(f"CDS that will anchor an island: {int(anchors['n_members'].sum()):,}")
    print(f"  of which POINT_ONLY (gyrA and friends): "
          f"{int(ann.loc[ann.anchor_tier=='POINT_ONLY','n_members'].sum()):,}")
    print(f"  of which CORE_AMR/SOFTCORE_AMR (backbone): "
          f"{int(ann.loc[ann.anchor_tier.isin(['CORE_AMR','SOFTCORE_AMR']),'n_members'].sum()):,}")
    print(f"  of which SHELL_AMR/MOBILE_AMR (cargo): "
          f"{int(ann.loc[ann.anchor_tier.isin(['SHELL_AMR','MOBILE_AMR']),'n_members'].sum()):,}")

    print("\nthreshold sweep (calibration_failed should be 0):")
    print(sweep.to_string(index=False))

    print("\nthe 20 most prevalent AMR clusters (the ones that decide the threshold):")
    top = amr_like.nlargest(20, "prevalence_eff")
    for r in top.itertuples():
        print(f"  {r.cluster_id} prev={r.prevalence:.3f} sp={r.prevalence_species_max:.3f} "
              f"({r.prevalence_species[:22]:22s}) {r.anchor_tier:13s} "
              f"n={r.n_members:6,} hits={r.n_amr_hits:5,} "
              f"{(r.amr_symbols or r.point_symbols)[:34]:34s} {r.product_top[:34]}")
    log("done")


if __name__ == "__main__":
    main()
