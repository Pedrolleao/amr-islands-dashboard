# Citrobacter freundii – AMR annotation of Pathogen Detection group PDG000000039.1314

Built by `latest_snps/scripts/run_all.py --group Citrobacter_freundii`. Every number below was computed from the files in this folder; re-run `scripts/10_group_readme.py --group Citrobacter_freundii` after a rebuild.

Reference catalog pinned at AMRFinderPlus **2026-08-07.1**. MicroBIGG-E export: `microbigge-C_freundii.tsv`.

## The dataset

| | count | share |
|---|---|---|
| isolates in the group | 7,287 | 100% |
| with a public assembly (`asm_acc`) | 6,063 | 83.2% |
| genomes downloaded here | 6,063 | 83.2% |
| covered by MicroBIGG-E (gene copies) | 6,063 | 83.2% |
| in an SNP cluster | 3,890 | 53.4% |
| SNP clusters | 925 | |
| carrying >=1 AMR gene | 7,286 | 100.0% |
| with a measured antibiogram | 75 | 1.0% |

**The group is not one species.** Composition of `scientific_name`:

| species | isolates |
|---|---|
| Citrobacter freundii | 5,216 (71.6%) |
| Citrobacter koseri | 665 (9.1%) |
| Citrobacter braakii | 483 (6.6%) |
| Citrobacter amalonaticus | 294 (4.0%) |
| Citrobacter farmeri | 184 (2.5%) |
| 7 further taxa | 445 |

Filter on `scientific_name` before reporting anything as a species-level figure.

## Resistance overview

Share of isolates carrying at least one gene of each class, and the genes behind it.

| class | isolates | main genes |
|---|---|---|
| BETA-LACTAM | 99.7% | `blaTEM-1` (2,900), `blaCMY` (2,254), `blaOXA-1` (1,793), `blaKPC-2` (1,526) |
| QUINOLONE | 79.5% | `gyrA_T83I` (3,788), `aac(6')-Ib-cr5` (2,040), `qnrS1` (850), `qnrB38` (641) |
| AMINOGLYCOSIDE | 70.2% | `aadA1` (2,270), `aac(6')-Ib-cr5` (2,040), `aph(6)-Id` (1,631), `aph(3'')-Ib` (1,578) |
| SULFONAMIDE | 62.0% | `sul1` (4,012), `sul2` (2,227), `sul3` (24) |
| TRIMETHOPRIM | 53.9% | `dfrA1` (1,393), `dfrA12` (1,070), `dfrA14` (894), `dfrA19` (411) |
| PHENICOL | 45.6% | `catB3` (2,022), `catA2` (778), `oqxB` (464), `floR` (391) |
| TETRACYCLINE | 39.6% | `tet(A)` (1,363), `tet(D)` (885), `oqxB` (464), `tet(B)` (360) |
| MACROLIDE | 39.6% | `mph(A)` (2,375), `mrx(A)` (1,940), `msr(E)` (499), `mph(E)` (484) |
| RIFAMYCIN | 25.0% | `arr-3` (1,533), `arr` (223), `arr-2` (128) |
| BLEOMYCIN | 19.5% | `ble` (1,408), `bleO` (26) |
| STREPTOTHRICIN | 14.1% | `sat2` (1,030) |
| FOSFOMYCIN | 9.1% | `fosA7` (296), `fosA7.9` (130), `fosA3` (123), `fos` (62) |

> **Intrinsic, not acquired.** BETA-LACTAM sit above 90% of isolates. A class carried by nearly every isolate is part of the species' own chromosome, not a resistance event; exclude these before comparing isolates or clusters.

## Organism-specific settings (`scripts/groups.json`)

| setting | value | meaning |
|---|---|---|
| `pdg` | `PDG000000039.1314` | Pathogen Detection release, prefix of the file names |
| `taxgroup` | `Citrobacter_freundii` | `whitelisted_taxa` used to pick point-mutation rows |
| `key_subclasses` | `CARBAPENEM`, `CEPHALOSPORIN`, `TIGECYCLINE` | get their own 0/1 column |
| `marker_subclass` | `CARBAPENEM` | genes of top concern, listed per isolate in `carbapenemase_genes` |
| `antibiogram_query` | `antibiogram[filter] AND Citrobacter[Organism]` | BioSample search for measured susceptibility |

## Read this before interpreting

- **Point mutations**: 3,828 calls. Resistance-conferring substitutions are screened for this organism, so they appear alongside acquired genes.
- **Family-level calls**: 1,117 of 80,439 AMR calls are profile (`HMM`) hits, which name a gene family rather than an allele. They can flag a housekeeping homologue as if it were a resistance gene, so treat a family-level call as a lead, not a result.
- **Measured phenotypes** exist for 75 isolates only (1.0%). Enough to spot check a genotype, too few to validate the profiles.
- **Gene absence is not susceptibility.** 1,224 isolates have no public assembly, and every call depends on what the pinned catalog knows.

## Files

```
AMR/ Metadata/ Clusters/ SNP_trees/ Exceptions/   original NCBI download - never modified
AMR_annotation/    gene -> antibiotic class, per-isolate profiles, per-cluster summaries
Genomes/           assemblies + proteins + annotation, and the manifest of what was downloaded
AMR_sequences/     protein / DNA sequences per gene family, for alignment
Antibiograms/      measured susceptibility from NCBI BioSample
```

Column-by-column documentation of every table: [`../docs/DATA_DOCUMENTATION.md`](../docs/DATA_DOCUMENTATION.md). How the pipeline was built and how to add another organism: [`../docs/WORKPLAN_new_species.md`](../docs/WORKPLAN_new_species.md).

## Rebuilding

```bash
python3 scripts/run_all.py --group Citrobacter_freundii                 # steps 01-05, the tables
python3 scripts/06_genome_manifest.py --group Citrobacter_freundii      # decide which accession per isolate
python3 scripts/07_download_genomes.py --group Citrobacter_freundii     # long, resumable
python3 scripts/08_extract_amr_sequences.py --group Citrobacter_freundii
python3 scripts/09_fetch_antibiograms.py --group Citrobacter_freundii
python3 scripts/10_group_readme.py --group Citrobacter_freundii         # this file
```
