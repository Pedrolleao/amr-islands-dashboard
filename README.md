# AMR annotation project – folder guide

> **Live dashboard: <https://pedrolleao.github.io/amr-islands-dashboard/>**
> Interactive, self-contained pages over the *Citrobacter freundii* group (7,287 isolates):
> [Islands](https://pedrolleao.github.io/amr-islands-dashboard/Citrobacter_freundii/Dashboard/islands.html) ·
> [Gene families](https://pedrolleao.github.io/amr-islands-dashboard/Citrobacter_freundii/Dashboard/clusters.html) ·
> [Antibiograms](https://pedrolleao.github.io/amr-islands-dashboard/Citrobacter_freundii/Dashboard/antibiograms.html).
> No install needed; each page is a few MB and takes a moment to load. Source and build steps: [`Citrobacter_freundii/Dashboard/README.md`](Citrobacter_freundii/Dashboard/README.md).

One folder per NCBI Pathogen Detection organism group; the same pipeline builds all of them.

```
latest_snps/
  Citrobacter_freundii/            NCBI Pathogen Detection group PDG000000039.1314 (7,287 isolates, whole genus)
  Citrobacter_portucalensis/       NCBI Pathogen Detection group PDG000000147.9 (581 isolates) – same sub-folders, see its README
  Clostridioides_difficile/        NCBI Pathogen Detection group PDG000000045.850 (39,708 isolates)
  Clostridium_botulinum/           NCBI Pathogen Detection group PDG000000059.230 (1,676 isolates)
  Clostridium_perfringens/         NCBI Pathogen Detection group PDG000000061.366 (5,706 isolates)
    AMR/ Metadata/ Clusters/ SNP_trees/ Exceptions/   original NCBI download – never modified
    AMR_annotation/                gene -> antibiotic class tables, resistance profiles, cluster summaries
    Genomes/                       genome + proteins + annotation of every isolate with a public assembly
    AMR_sequences/                 protein / DNA sequences of the AMR, stress and virulence genes, per gene family
    Antibiograms/                  lab-measured susceptibility tables (MIC / disk) from NCBI BioSample, where any exist
  reference/AMRFinderPlus_<version>/   pinned NCBI reference: gene catalog, gene hierarchy, reference proteins
  docs/                            full documentation + workplan for adding another organism
  scripts/                         pipeline (numbered in running order)
  tools/                           NCBI `datasets` / `dataformat` programs
../microbigge-C_freundii.tsv                               MicroBIGG-E export (input, freundii group)
../microbigge-Clostridioides_difficile.tsv                 MicroBIGG-E export (input)
../microbigge-Clostridium_botulinum.tsv                    MicroBIGG-E export (input; genus-wide, shared with perfringens)
../microbigge-Clostridium_perfringens.tsv                  MicroBIGG-E export (input; identical file, genus-wide)
../Citrobacter_portucalensis/microbigge-C_portucalensis.tsv  MicroBIGG-E export (input, genus-wide; filtered to the portucalensis group)
```

Each generated folder has its own `README.md` describing files and columns.

Workplan for the stress / virulence cluster analyses: [`docs/WORKPLAN_stress_virulence.md`](docs/WORKPLAN_stress_virulence.md)

Full documentation: [`docs/DATA_DOCUMENTATION.md`](docs/DATA_DOCUMENTATION.md) (what we have, how it links, checks,
limitations) · [`docs/WORKPLAN_new_species.md`](docs/WORKPLAN_new_species.md) (how it was done, how to repeat it for
another species).

## Pipeline

| step | script | makes |
|---|---|---|
| 1 | `01_fetch_reference.py` | `reference/` |
| 2–5 | `02_symbol_annotation.py` … `05_profiles_and_clusters.py` | `AMR_annotation/` |
| 6 | `06_genome_manifest.py` | `Genomes/manifest/` (which accession to download per isolate) |
| 7 | `07_download_genomes.py` | `Genomes/assemblies/` (~60 GB, resumable) |
| 8 | `08_extract_amr_sequences.py` | `AMR_sequences/` |
| 9 | `09_fetch_antibiograms.py` | `Antibiograms/` (needs only steps 1–5) |
| 10 | `10_group_readme.py` | `<Organism_group>/README.md`, written from the tables that exist |

Organism groups and their settings are listed in `scripts/groups.json`; every script takes
`--group <Organism_group>` (default: `Citrobacter_freundii`). Besides the release id and the input
paths, a group sets what to report for that organism: `key_subclasses`, the `marker_subclass` singled
out per isolate (carbapenemases for Enterobacterales, `null` where the organism has no such gene) and
`spot_checks`. Omitted keys fall back to the Enterobacterales defaults.

`python3 scripts/run_all.py [--group G]` runs steps 1–5; `python3 scripts/run_all.py --sequences` also runs 6–8.
Step 9 is quick and independent: `python3 scripts/09_fetch_antibiograms.py`.
`python3 scripts/run_islands.py` runs the island pipeline (11-23) and `python3 scripts/run_dashboard.py`
builds the three dashboard pages in `<group>/Dashboard/` from the finished analysis (30-36).
Requirements: Python 3 with pandas; `tools/datasets` from
https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/ .
