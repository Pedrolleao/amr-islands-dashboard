# Citrobacter portucalensis – AMR annotation of Pathogen Detection group PDG000000147.9

Built by `latest_snps/scripts/run_all.py --group Citrobacter_portucalensis`. Every number below was computed from the files in this folder; re-run `scripts/10_group_readme.py --group Citrobacter_portucalensis` after a rebuild.

Reference catalog pinned at AMRFinderPlus **2026-08-07.1**. MicroBIGG-E export: `microbigge-C_portucalensis.tsv`.

## The dataset

| | count | share |
|---|---|---|
| isolates in the group | 581 | 100% |
| with a public assembly (`asm_acc`) | 432 | 74.4% |
| genomes downloaded here | 432 | 74.4% |
| covered by MicroBIGG-E (gene copies) | 430 | 74.0% |
| in an SNP cluster | 212 | 36.5% |
| SNP clusters | 68 | |
| carrying >=1 AMR gene | 579 | 99.7% |
| with a measured antibiogram | 5 | 0.9% |

**The group is not one species.** Composition of `scientific_name`:

| species | isolates |
|---|---|
| Citrobacter portucalensis | 581 (100.0%) |

Filter on `scientific_name` before reporting anything as a species-level figure.

## Resistance overview

Share of isolates carrying at least one gene of each class, and the genes behind it.

| class | isolates | main genes |
|---|---|---|
| BETA-LACTAM | 99.7% | `blaCMY` (256), `blaTEM-1` (219), `blaNDM-1` (138), `blaOXA-1` (101) |
| QUINOLONE | 99.5% | `qnrB9` (206), `aac(6')-Ib-cr5` (122), `qnrB6` (103), `qnrS1` (68) |
| AMINOGLYCOSIDE | 71.1% | `aadA2` (178), `aph(3')-I` (171), `aac(3)-IId` (159), `aph(3'')-Ib` (140) |
| SULFONAMIDE | 57.3% | `sul1` (308), `sul2` (165), `sul3` (6) |
| TRIMETHOPRIM | 49.1% | `dfrA12` (149), `dfrA14` (51), `dfrA1` (37), `dfrA27` (37) |
| MACROLIDE | 45.4% | `mph(A)` (229), `mrx(A)` (229), `mph(E)` (55), `msr(E)` (55) |
| TETRACYCLINE | 38.6% | `tet(A)` (117), `tet(D)` (108), `tet(B)` (13), `tet(X4)` (10) |
| PHENICOL | 36.7% | `catB3` (124), `floR` (60), `catA2` (51), `catA1` (29) |
| BLEOMYCIN | 31.0% | `ble` (172), `bleO` (15) |
| RIFAMYCIN | 21.0% | `arr-3` (109), `arr` (8), `arr-2` (8) |
| STREPTOGRAMIN | 11.2% | `msr(E)` (55), `erm(B)` (7), `msr(A)` (2), `vat` (1) |
| FOSFOMYCIN | 9.0% | `fosA3` (33), `fosA7` (11), `fos` (6), `fosA7.5` (1) |

> **Intrinsic, not acquired.** BETA-LACTAM, QUINOLONE sit above 90% of isolates. A class carried by nearly every isolate is part of the species' own chromosome, not a resistance event; exclude these before comparing isolates or clusters.

## Organism-specific settings (`scripts/groups.json`)

| setting | value | meaning |
|---|---|---|
| `pdg` | `PDG000000147.9` | Pathogen Detection release, prefix of the file names |
| `taxgroup` | `Citrobacter_freundii` | `whitelisted_taxa` used to pick point-mutation rows |
| `key_subclasses` | `CARBAPENEM`, `CEPHALOSPORIN`, `TIGECYCLINE` | get their own 0/1 column |
| `marker_subclass` | `CARBAPENEM` | genes of top concern, listed per isolate in `carbapenemase_genes` |
| `antibiogram_query` | `antibiogram[filter] AND Citrobacter portucalensis[Organism]` | BioSample search for measured susceptibility |

## Read this before interpreting

- **No point mutations are screened for this organism**, so mutational resistance (for example in *gyrA*) is invisible here. That is a blind spot of the dataset, not an absence of resistance.
- **Family-level calls**: 194 of 6,283 AMR calls are profile (`HMM`) hits, which name a gene family rather than an allele. They can flag a housekeeping homologue as if it were a resistance gene, so treat a family-level call as a lead, not a result.
- **Measured phenotypes** exist for 5 isolates only (0.9%). Enough to spot check a genotype, too few to validate the profiles.
- **Gene absence is not susceptibility.** 149 isolates have no public assembly, and every call depends on what the pinned catalog knows.

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
python3 scripts/run_all.py --group Citrobacter_portucalensis                 # steps 01-05, the tables
python3 scripts/06_genome_manifest.py --group Citrobacter_portucalensis      # decide which accession per isolate
python3 scripts/07_download_genomes.py --group Citrobacter_portucalensis     # long, resumable
python3 scripts/08_extract_amr_sequences.py --group Citrobacter_portucalensis
python3 scripts/09_fetch_antibiograms.py --group Citrobacter_portucalensis
python3 scripts/10_group_readme.py --group Citrobacter_portucalensis         # this file
```
