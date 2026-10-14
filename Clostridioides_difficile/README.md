# Clostridioides difficile – AMR annotation of Pathogen Detection group PDG000000045.850

Built by `latest_snps/scripts/run_all.py --group Clostridioides_difficile`. Every number below was computed from the files in this folder; re-run `scripts/10_group_readme.py --group Clostridioides_difficile` after a rebuild.

Reference catalog pinned at AMRFinderPlus **2026-08-07.1**. MicroBIGG-E export: `microbigge-Clostridioides_difficile.tsv`.

## The dataset

| | count | share |
|---|---|---|
| isolates in the group | 39,708 | 100% |
| with a public assembly (`asm_acc`) | 31,930 | 80.4% |
| genomes downloaded here | 31,900 | 80.3% |
| assemblies NCBI no longer serves | 30 | GCA_061433475.1, GCA_061433495.1, GCA_061433535.1 … |
| covered by MicroBIGG-E (gene copies) | 31,758 | 80.0% |
| in an SNP cluster | 37,807 | 95.2% |
| SNP clusters | 864 | |
| carrying >=1 AMR gene | 39,649 | 99.9% |
| with a measured antibiogram | 105 | 0.3% |

**The group is not one species.** Composition of `scientific_name`:

| species | isolates |
|---|---|
| Clostridioides difficile | 39,536 (99.6%) |
| Paraclostridium sordellii | 71 (0.2%) |
| Paraclostridium bifermentans | 43 (0.1%) |
| Intestinibacter bartlettii | 31 (0.1%) |
| Peptacetobacter hiranonis | 8 (0.0%) |
| 6 further taxa | 19 |

Filter on `scientific_name` before reporting anything as a species-level figure.

## Resistance overview

Share of isolates carrying at least one gene of each class, and the genes behind it.

| class | isolates | main genes |
|---|---|---|
| LINCOSAMIDE | 99.7% | `cplR` (39,574), `erm(B)` (6,994), `cfr(C)` (956), `lsa` (797) |
| STREPTOGRAMIN | 99.7% | `cplR` (39,574), `erm(B)` (6,994), `cfr(C)` (956), `lsa` (797) |
| PLEUROMUTILIN | 99.7% | `cplR` (39,574), `cfr(C)` (956), `cfr(B)` (443), `cfr(E)` (205) |
| BETA-LACTAM | 99.5% | `blaAHM` (38,509), `blaCDD` (35,650), `blaR1` (22,440), `mreE_V497L` (8,048) |
| NITROIMIDAZOLE | 99.5% | `nimB-Cd` (39,498), `PnimB_G` (5,143) |
| GLYCOPEPTIDE | 98.7% | `vanT-Cd` (33,357), `vanR-Cd` (33,354), `vanG` (33,352), `vanS-Cd` (33,344) |
| QUINOLONE | 55.3% | `gyrA_T82I` (8,033), `gyrB_I139R` (6,646), `gyrB_S366A` (3,108), `gyrB_S416A` (3,107) |
| MACROLIDE | 19.0% | `erm(B)` (6,994), `erm(52)` (296), `mef(H)` (140), `mef(A)` (121) |
| TETRACYCLINE | 18.5% | `tet(M)` (6,461), `tet(40)` (958), `tetA(P)` (670), `tetB(P)` (664) |
| AMINOGLYCOSIDE | 18.5% | `aac(6')-Ie/aph(2'')-Ia` (3,282), `aadE` (3,261), `ant(6)-Ia` (1,686), `aph(3')-IIIa` (1,070) |
| RIFAMYCIN | 6.3% | `rpoB_R505K` (2,334), `rpoB_H502N` (981), `arr` (35), `rpoB_S550F` (12) |
| PHENICOL | 4.6% | `cfr(C)` (956), `cfr(B)` (443), `catP` (269), `cfr(E)` (205) |

> **Intrinsic, not acquired.** LINCOSAMIDE, STREPTOGRAMIN, PLEUROMUTILIN, BETA-LACTAM, NITROIMIDAZOLE, GLYCOPEPTIDE sit above 90% of isolates. A class carried by nearly every isolate is part of the species' own chromosome, not a resistance event; exclude these before comparing isolates or clusters.

## Organism-specific settings (`scripts/groups.json`)

| setting | value | meaning |
|---|---|---|
| `pdg` | `PDG000000045.850` | Pathogen Detection release, prefix of the file names |
| `taxgroup` | `Clostridioides_difficile` | `whitelisted_taxa` used to pick point-mutation rows |
| `key_subclasses` | `VANCOMYCIN`, `NITROIMIDAZOLE`, `QUINOLONE`, `ERYTHROMYCIN`, `TETRACYCLINE`, `LINEZOLID`, `RIFAMYCIN` | get their own 0/1 column |
| `marker_subclass` | `QUINOLONE` | genes of top concern, listed per isolate in `fluoroquinolone_genes` |
| `antibiogram_query` | `antibiogram[filter] AND Clostridioides difficile[Organism]` | BioSample search for measured susceptibility |

## Read this before interpreting

- **Point mutations**: 56,837 calls. Resistance-conferring substitutions are screened for this organism, so they appear alongside acquired genes.
- **Family-level calls**: 29,197 of 445,545 AMR calls are profile (`HMM`) hits, which name a gene family rather than an allele. They can flag a housekeeping homologue as if it were a resistance gene, so treat a family-level call as a lead, not a result.
- **Measured phenotypes** exist for 105 isolates only (0.3%). Enough to spot check a genotype, too few to validate the profiles.
- **Gene absence is not susceptibility.** 7,778 isolates have no public assembly, and every call depends on what the pinned catalog knows.

## Files

```
AMR/ Metadata/ Clusters/ SNP_trees/ Exceptions/   original NCBI download - never modified
AMR_annotation/    gene -> antibiotic class, per-isolate profiles, per-cluster summaries
Genomes/           assemblies + proteins + annotation, and the manifest of what was downloaded
AMR_sequences/     protein / DNA sequences per gene family, for alignment
Antibiograms/      measured susceptibility from NCBI BioSample (105 isolates)
```

Column-by-column documentation of every table: [`../docs/DATA_DOCUMENTATION.md`](../docs/DATA_DOCUMENTATION.md). How the pipeline was built and how to add another organism: [`../docs/WORKPLAN_new_species.md`](../docs/WORKPLAN_new_species.md).

## Rebuilding

```bash
python3 scripts/run_all.py --group Clostridioides_difficile                 # steps 01-05, the tables
python3 scripts/06_genome_manifest.py --group Clostridioides_difficile      # decide which accession per isolate
python3 scripts/07_download_genomes.py --group Clostridioides_difficile     # long, resumable
python3 scripts/08_extract_amr_sequences.py --group Clostridioides_difficile
python3 scripts/09_fetch_antibiograms.py --group Clostridioides_difficile
python3 scripts/10_group_readme.py --group Clostridioides_difficile         # this file
```
