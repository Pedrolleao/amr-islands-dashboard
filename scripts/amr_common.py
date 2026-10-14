"""Shared paths, organism-group settings and parsing helpers for the AMR annotation scripts.

Every script accepts  --group <Organism_group>  (or the environment variable AMR_GROUP);
the groups and their settings are listed in scripts/groups.json.
"""
import json
import os
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]                # latest_snps/
GROUPS = {k: v for k, v in json.loads((Path(__file__).with_name("groups.json")).read_text()).items()
          if not k.startswith("_")}


def _selected_group():
    if "--group" in sys.argv:
        i = sys.argv.index("--group")
        name = sys.argv[i + 1]
        del sys.argv[i:i + 2]
        return name
    return os.environ.get("AMR_GROUP", next(iter(GROUPS)))


GROUP = _selected_group()
if GROUP not in GROUPS:
    raise SystemExit(f"unknown group {GROUP!r}; known groups: {', '.join(GROUPS)} (scripts/groups.json)")
os.environ["AMR_GROUP"] = GROUP                           # inherited by scripts started from run_all.py
_cfg = GROUPS[GROUP]

DATA = ROOT / GROUP
PDG = _cfg["pdg"]                                         # Pathogen Detection release (prefix of file names)
TAXGROUP = _cfg["taxgroup"]                               # whitelisted_taxa value for point mutations
ANTIBIOGRAM_QUERY = _cfg["antibiogram_query"]

# Organism-specific reporting settings of 05_profiles_and_clusters.py. The defaults are the
# Enterobacterales ones used for the Citrobacter groups; other organisms override them in groups.json.
KEY_SUBCLASSES = _cfg.get("key_subclasses", ["CARBAPENEM", "CEPHALOSPORIN", "TIGECYCLINE"])
# subclass singled out as its own gene list ("which isolate carries which gene of top concern")
MARKER_SUBCLASS = _cfg.get("marker_subclass", "CARBAPENEM")
MARKER_COLUMN = _cfg.get("marker_column", "carbapenemase_genes")
# symbol -> (class, subclass) assertions; genes absent from a group are skipped
SPOT_CHECKS = _cfg.get("spot_checks", {
    "blaKPC-2": ["BETA-LACTAM", "CARBAPENEM"], "blaNDM-1": ["BETA-LACTAM", "CARBAPENEM"],
    "sul1": ["SULFONAMIDE", "SULFONAMIDE"], "gyrA_T83I": ["QUINOLONE", "QUINOLONE"],
    "mcr-9.1": ["COLISTIN", "COLISTIN"], "tet(A)": ["TETRACYCLINE", "TETRACYCLINE"]})

METADATA_TSV = DATA / "AMR" / f"{PDG}.amr.metadata.tsv"
ALL_ISOLATES_TSV = DATA / "Clusters" / f"{PDG}.reference_target.all_isolates.tsv"
MICROBIGGE_TSV = (ROOT / _cfg["microbigge_tsv"]).resolve()

REF_BASE_URL = "https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest/"
REF_FILES = ["version.txt", "ReferenceGeneCatalog.txt", "fam.tsv", "AMRProt.fa"]
REF_ROOT = ROOT / "reference"

OUT = DATA / "AMR_annotation"
GENOMES = DATA / "Genomes"                     # assemblies downloaded with NCBI datasets
TOOLS = ROOT / "tools"                          # NCBI datasets / dataformat binaries
PANGENOME = DATA / "Pangenome"                  # gene index + protein clusters (scripts 11-15)
ISLANDS = DATA / "Islands"                      # AMR islands and their figures (scripts 16-22)
VARIATION = DATA / "Variation"                  # per-position variation of AMR clusters (script 23)
STRUCTURES = DATA / "Structures"                # ESMFold structures + py3Dmol views (scripts 24-25)
DOMAINS = DATA / "Domains"                      # Pfam domains + Foldseek structure hits (scripts 27-28)
DASHBOARD = DATA / "Dashboard"                  # dashboard payloads + built pages (scripts 30-35)

# Reference databases that live outside the project, installed for other work on this machine.
PFAM_HMM = Path("/home/dmin/Program/interproscan-5.76-107.0/data/pfam/38.0/pfam_a.hmm")
PFAM_DAT = PFAM_HMM.with_name("pfam_a.dat")
FOLDSEEK = Path("/home/dmin/miniconda3/envs/foldseek/bin/foldseek")
# Structure databases already on this machine. afdb_sp is AlphaFold Swiss-Prot v6 (2025_03),
# 590,183 models whose headers carry a protein NAME, so a hit reads as a function and not just an id.
_REFS = Path("/home/dmin/Database/NCBI_Archaea_02-2025/Defense_finder/ArchaeaHQ-DefenseFinder"
             "/Discovery/Explore/validation/refs")
FOLDSEEK_DBS = {
    "pdb": Path("/home/dmin/.local/share/archaeadefense/foldseek/pdb"),   # 324,204 experimental
    "afdb_sp": _REFS / "afdb_sp",                                         # 590,183 predicted, named
}
FOLDSEEK_PDB = FOLDSEEK_DBS["pdb"]     # kept: earlier runs referred to it by this name


def _selected_set():
    """--set {amr,stress,virulence,all}: which group of protein clusters steps 23-28 work on.

    `all` is every cluster of the pangenome, for the domain annotation of steps 27-28 only; the
    per-position and structure steps would be absurd on 200,000 clusters.
    """
    if "--set" in sys.argv:
        i = sys.argv.index("--set")
        name = sys.argv[i + 1]
        del sys.argv[i:i + 2]
    else:
        name = os.environ.get("AMR_SET", "amr")
    if name not in ("amr", "stress", "virulence", "all"):
        raise SystemExit(f"unknown --set {name!r}; use amr, stress, virulence or all")
    return name


SET = _selected_set()
os.environ["AMR_SET"] = SET


def set_dir(base):
    """Output folder of the current cluster set: Variation/ for amr, Variation_stress/ for stress ..."""
    return base if SET == "amr" else base.with_name(f"{base.name}_{SET}")

# Prefix of the protein-cluster ids (CFR-000001). NOT to be confused with the Pathogen Detection
# SNP clusters (PDS_acc) handled by 05_profiles_and_clusters.py - those are a different thing.
CLUSTER_PREFIX = _cfg.get("cluster_prefix", "".join(w[:1] for w in GROUP.split("_")).upper())

# The ONE colour dimension used by every island table, figure and (later) the browser.
# A gene gets exactly one role; precedence is AMR anchor > stress > virulence > MGE > rest.
DRAW_ROLE_COLOURS = {
    "amr_anchor_mobile": "#c0392b",    # acquired resistance gene - the point of the figure
    "amr_anchor_core": "#d98880",      # AMR-annotated but backbone (chromosomal ampC ...)
    "amr_anchor_point": "#7b241c",     # resistance by point mutation in an essential gene (gyrA)
    "stress": "#2471a3",               # metal / biocide resistance
    "virulence": "#7d3c98",
    "mge": "#f1c40f",                  # transposase, integrase, IS element (pale when pseudogene)
    "mge_pseudo": "#f9e79f",
    "hypothetical": "#f2f3f4",
    "other": "#bdc3c7",
}

# metadata column -> element type as used by AMRFinderPlus / MicroBIGG-E
GENOTYPE_COLS = {
    "AMR_genotypes": "AMR",
    "stress_genotypes": "STRESS",
    "virulence_genotypes": "VIRULENCE",
}
COUNT_COLS = {
    "AMR": "number_amr_genes",
    "STRESS": "number_stress_genes",
    "VIRULENCE": "number_virulence_genes",
}


def read_tsv(path, quoted=False, **kw):
    """Read a TSV with every column as string and no NA guessing ('NULL' -> '').

    quoted=True for the Pathogen Detection metadata, which wraps list fields in double quotes.
    """
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False,
                     quoting=0 if quoted else 3, encoding="utf-8-sig", **kw)
    return df.replace("NULL", "")


def write_tsv(df, name, folder=None):
    """Write a TSV into `folder` (default: the group's AMR_annotation/)."""
    folder = OUT if folder is None else folder
    folder.mkdir(parents=True, exist_ok=True)
    path = folder / name
    df.to_csv(path, sep="\t", index=False)
    print(f"wrote {path.relative_to(ROOT)}  ({len(df):,} rows x {df.shape[1]} cols)")
    return path


def ref_dir():
    """Directory of the pinned AMRFinderPlus database (latest one downloaded)."""
    dirs = sorted(REF_ROOT.glob("AMRFinderPlus_*"))
    if not dirs:
        raise SystemExit("No reference found - run 01_fetch_reference.py first")
    return dirs[-1]


def db_version():
    return (ref_dir() / "version.txt").read_text().strip()


def load_metadata():
    return read_tsv(METADATA_TSV, quoted=True)


def split_genotypes(value):
    """'blaCMY=PARTIAL,sul1' -> [('blaCMY', 'PARTIAL'), ('sul1', 'complete')]"""
    out = []
    for token in value.split(","):
        token = token.strip()
        if not token:
            continue
        symbol, _, qualifier = token.partition("=")
        out.append((symbol, qualifier or "complete"))
    return out


def explode_genotypes(md):
    """One row per isolate x element from the three genotype columns of the metadata."""
    records = []
    for row in md.itertuples(index=False):
        core = {s for s, _ in split_genotypes(row.AMR_genotypes_core)}
        for col, etype in GENOTYPE_COLS.items():
            for symbol, qualifier in split_genotypes(getattr(row, col)):
                records.append((row.target_acc, etype, symbol, qualifier,
                                int(etype == "AMR" and symbol in core)))
    return pd.DataFrame(records, columns=["target_acc", "element_type", "symbol",
                                          "call_qualifier", "is_core"])


def load_microbigge():
    mb = read_tsv(MICROBIGGE_TSV)
    mb.columns = (mb.columns.str.lstrip("#").str.strip().str.lower()
                  .str.replace("% ", "pct_").str.replace(" ", "_"))
    mb = mb.rename(columns={"isolate": "target_acc", "biosample": "biosample_acc",
                            "type": "element_type", "element_symbol": "symbol",
                            "pct_coverage_of_reference": "pct_coverage",
                            "pct_identity_to_reference": "pct_identity"})
    # an export can be broader than the group (e.g. whole genus): keep the group's isolates only
    in_group = mb["target_acc"].isin(set(load_metadata()["target_acc"]))
    if not in_group.all():
        print(f"MicroBIGG-E export: kept {int(in_group.sum()):,} of {len(mb):,} rows "
              f"({mb.loc[in_group, 'target_acc'].nunique():,} isolates of {GROUP})")
    return mb[in_group].reset_index(drop=True)


def split_classes(value):
    """'AMINOGLYCOSIDE/QUINOLONE' -> ['AMINOGLYCOSIDE', 'QUINOLONE']"""
    return [c for c in value.split("/") if c]
