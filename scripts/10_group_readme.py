"""Write <Organism_group>/README.md from the tables that were actually produced for that group.

Everything in the text is computed from the group's own files, so re-running it after a new
Pathogen Detection release keeps the numbers honest. Nothing is hard-coded per organism:
what differs between organisms comes from scripts/groups.json.

    python3 scripts/10_group_readme.py --group <Organism_group>
"""
import json

import pandas as pd

from amr_common import (ANTIBIOGRAM_QUERY, DATA, GENOMES, GROUP, KEY_SUBCLASSES, MARKER_COLUMN,
                        MARKER_SUBCLASS, MICROBIGGE_TSV, OUT, PDG, ROOT, TAXGROUP, db_version,
                        load_metadata, read_tsv, split_classes)

# a class carried by at least this share of isolates is flagged as probably intrinsic / near-universal
INTRINSIC_PCT = 90.0
TOP_GENES = 4


def opt(path, **kw):
    return read_tsv(path, **kw) if path.exists() else None


def fmt(n):
    return f"{n:,}"


def pct(n, total):
    return f"{100 * n / total:.1f}%" if total else "n/a"


def label_map(amr, column):
    """profile column suffix -> the class string as it is written in the tables.

    05_profiles_and_clusters.py lowercases and replaces ' ' and '-' by '_' when it builds the
    column names, which is not reversible ('BETA-LACTAM' -> 'beta_lactam'), so map it back from
    the values actually present.
    """
    values = {c for v in amr[column] for c in split_classes(v)}
    return {v.lower().replace(" ", "_").replace("-", "_"): v for v in values}


def genes_behind(amr, value, column):
    """Top gene symbols responsible for a class / subclass, with the number of isolates."""
    m = amr[amr[column].map(lambda v: value in split_classes(v))]
    top = (m.drop_duplicates(["target_acc", "symbol"]).groupby("symbol")["target_acc"].size()
           .sort_values(ascending=False).head(TOP_GENES))
    return ", ".join(f"`{s}` ({fmt(n)})" for s, n in top.items())


def main():
    md = load_metadata()
    long = read_tsv(OUT / "isolate_gene_long.tsv")
    profile = read_tsv(OUT / "isolate_resistance_profile.tsv")
    clusters = read_tsv(OUT / "cluster_amr_summary.tsv")
    hits = opt(OUT / "microbigge_hits.tsv")
    manifest = opt(GENOMES / "manifest" / "genome_manifest.tsv")
    seq_index = opt(DATA / "AMR_sequences" / "sequence_index.tsv")
    ab = opt(DATA / "Antibiograms" / "antibiogram_full.tsv")

    amr = long[long["element_type"] == "AMR"]
    n_iso = len(md)
    n_asm = int((md["asm_acc"] != "").sum())
    n_clustered = int((profile["PDS_acc"] != "").sum())
    n_downloaded = len(list((GENOMES / "assemblies").glob("GC*"))) if (GENOMES / "assemblies").exists() else 0
    missing_file = GENOMES / "manifest" / "download_missing.txt"
    missing = [x for x in missing_file.read_text().split() if x] if missing_file.exists() else []

    species = md["scientific_name"].str.split().str[:2].str.join(" ").value_counts()
    class_cols = [c for c in profile.columns if c.startswith("class_")]
    prevalence = (profile[class_cols].astype(int).mean() * 100).sort_values(ascending=False)

    L = []
    a = L.append
    a(f"# {GROUP.replace('_', ' ')} – AMR annotation of Pathogen Detection group {PDG}\n")
    a(f"Built by `latest_snps/scripts/run_all.py --group {GROUP}`. Every number below was computed from the "
      f"files in this folder; re-run `scripts/10_group_readme.py --group {GROUP}` after a rebuild.\n")
    a(f"Reference catalog pinned at AMRFinderPlus **{db_version()}**. "
      f"MicroBIGG-E export: `{MICROBIGGE_TSV.name}`.\n")

    a("## The dataset\n")
    a("| | count | share |")
    a("|---|---|---|")
    a(f"| isolates in the group | {fmt(n_iso)} | 100% |")
    a(f"| with a public assembly (`asm_acc`) | {fmt(n_asm)} | {pct(n_asm, n_iso)} |")
    a(f"| genomes downloaded here | {fmt(n_downloaded)} | {pct(n_downloaded, n_iso)} |")
    if missing:
        a(f"| assemblies NCBI no longer serves | {fmt(len(missing))} | {', '.join(missing[:3])}"
          f"{' …' if len(missing) > 3 else ''} |")
    if hits is not None:
        n_mb = hits["target_acc"].nunique()
        a(f"| covered by MicroBIGG-E (gene copies) | {fmt(n_mb)} | {pct(n_mb, n_iso)} |")
    a(f"| in an SNP cluster | {fmt(n_clustered)} | {pct(n_clustered, n_iso)} |")
    a(f"| SNP clusters | {fmt(len(clusters))} | |")
    a(f"| carrying >=1 AMR gene | {fmt(amr['target_acc'].nunique())} | "
      f"{pct(amr['target_acc'].nunique(), n_iso)} |")
    n_ab = ab["target_acc"].nunique() if ab is not None and len(ab) else 0
    a(f"| with a measured antibiogram | {fmt(n_ab)} | {pct(n_ab, n_iso)} |")
    a("")

    a("**The group is not one species.** Composition of `scientific_name`:\n")
    a("| species | isolates |")
    a("|---|---|")
    for name, n in species.head(5).items():
        a(f"| {name} | {fmt(n)} ({pct(n, n_iso)}) |")
    if len(species) > 5:
        n_rest = len(species) - 5
        a(f"| {n_rest} further {'taxon' if n_rest == 1 else 'taxa'} | {fmt(int(species[5:].sum()))} |")
    a("\nFilter on `scientific_name` before reporting anything as a species-level figure.\n")

    a("## Resistance overview\n")
    a("Share of isolates carrying at least one gene of each class, and the genes behind it.\n")
    a("| class | isolates | main genes |")
    a("|---|---|---|")
    classes = label_map(amr, "class")
    for col, p in prevalence[prevalence > 0].head(12).items():
        value = classes.get(col[len("class_"):], col[len("class_"):].replace("_", " ").upper())
        a(f"| {value} | {p:.1f}% | {genes_behind(amr, value, 'class')} |")
    a("")

    intrinsic = [classes.get(c[len("class_"):], c[len("class_"):]) for c, p in prevalence.items()
                 if p >= INTRINSIC_PCT]
    if intrinsic:
        a(f"> **Intrinsic, not acquired.** {', '.join(intrinsic)} sit above {INTRINSIC_PCT:.0f}% of isolates. "
          "A class carried by nearly every isolate is part of the species' own chromosome, not a resistance "
          "event; exclude these before comparing isolates or clusters.\n")

    a("## Organism-specific settings (`scripts/groups.json`)\n")
    a("| setting | value | meaning |")
    a("|---|---|---|")
    a(f"| `pdg` | `{PDG}` | Pathogen Detection release, prefix of the file names |")
    a(f"| `taxgroup` | `{TAXGROUP}` | `whitelisted_taxa` used to pick point-mutation rows |")
    a(f"| `key_subclasses` | {', '.join(f'`{s}`' for s in KEY_SUBCLASSES)} | get their own 0/1 column |")
    if MARKER_SUBCLASS:
        a(f"| `marker_subclass` | `{MARKER_SUBCLASS}` | genes of top concern, listed per isolate in "
          f"`{MARKER_COLUMN}` |")
    else:
        a("| `marker_subclass` | *(none)* | no gene of top concern in this organism, so no marker column |")
    a(f"| `antibiogram_query` | `{ANTIBIOGRAM_QUERY}` | BioSample search for measured susceptibility |")
    a("")

    n_point = int((amr["call_qualifier"] == "POINT").sum())
    n_hmm = int((amr["call_qualifier"] == "HMM").sum())
    a("## Read this before interpreting\n")
    if n_point:
        a(f"- **Point mutations**: {fmt(n_point)} calls. Resistance-conferring substitutions are screened for "
          "this organism, so they appear alongside acquired genes.")
    else:
        a("- **No point mutations are screened for this organism**, so mutational resistance (for example in "
          "*gyrA*) is invisible here. That is a blind spot of the dataset, not an absence of resistance.")
    a(f"- **Family-level calls**: {fmt(n_hmm)} of {fmt(len(amr))} AMR calls are profile (`HMM`) hits, which "
      "name a gene family rather than an allele. They can flag a housekeeping homologue as if it were a "
      "resistance gene, so treat a family-level call as a lead, not a result.")
    if n_ab:
        a(f"- **Measured phenotypes** exist for {fmt(n_ab)} isolates only ({pct(n_ab, n_iso)}). Enough to spot "
          "check a genotype, too few to validate the profiles.")
    else:
        a("- **No measured antibiograms exist** for this group in NCBI BioSample, so nothing here can be "
          "checked against a laboratory result. The profiles are gene content only.")
    a(f"- **Gene absence is not susceptibility.** {fmt(n_iso - n_asm)} isolates have no public assembly, and "
      "every call depends on what the pinned catalog knows.\n")

    a("## Files\n")
    a("```")
    a("AMR/ Metadata/ Clusters/ SNP_trees/ Exceptions/   original NCBI download - never modified")
    a("AMR_annotation/    gene -> antibiotic class, per-isolate profiles, per-cluster summaries")
    if manifest is not None:
        a("Genomes/           assemblies + proteins + annotation, and the manifest of what was downloaded")
    if seq_index is not None:
        a("AMR_sequences/     protein / DNA sequences per gene family, for alignment")
    if (DATA / "Antibiograms").exists():
        a(f"Antibiograms/      measured susceptibility from NCBI BioSample"
          f"{f' ({fmt(n_ab)} isolates)' if n_ab else ' - searched, none exists for this group'}")
    a("```\n")
    a("Column-by-column documentation of every table: [`../docs/DATA_DOCUMENTATION.md`]"
      "(../docs/DATA_DOCUMENTATION.md). How the pipeline was built and how to add another organism: "
      "[`../docs/WORKPLAN_new_species.md`](../docs/WORKPLAN_new_species.md).\n")

    a("## Rebuilding\n")
    a("```bash")
    a(f"python3 scripts/run_all.py --group {GROUP}                 # steps 01-05, the tables")
    a(f"python3 scripts/06_genome_manifest.py --group {GROUP}      # decide which accession per isolate")
    a(f"python3 scripts/07_download_genomes.py --group {GROUP}     # long, resumable")
    a(f"python3 scripts/08_extract_amr_sequences.py --group {GROUP}")
    a(f"python3 scripts/09_fetch_antibiograms.py --group {GROUP}")
    a(f"python3 scripts/10_group_readme.py --group {GROUP}         # this file")
    a("```")

    path = DATA / "README.md"
    path.write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"wrote {path.relative_to(ROOT)}  ({len(L)} lines)")


if __name__ == "__main__":
    main()
