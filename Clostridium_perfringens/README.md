# Clostridium perfringens – AMR annotation of Pathogen Detection group PDG000000061.366

Built by `latest_snps/scripts/run_all.py --group Clostridium_perfringens`. Every number below was computed from the files in this folder; re-run `scripts/10_group_readme.py --group Clostridium_perfringens` after a rebuild.

Reference catalog pinned at AMRFinderPlus **2026-08-07.1**. MicroBIGG-E export: `microbigge-Clostridium_perfringens.tsv`.

## The dataset

| | count | share |
|---|---|---|
| isolates in the group | 5,706 | 100% |
| with a public assembly (`asm_acc`) | 5,465 | 95.8% |
| genomes downloaded here | 5,464 | 95.8% |
| assemblies NCBI no longer serves | 1 | GCA_061415275.1 |
| covered by MicroBIGG-E (gene copies) | 5,393 | 94.5% |
| in an SNP cluster | 3,535 | 62.0% |
| SNP clusters | 815 | |
| carrying >=1 AMR gene | 5,043 | 88.4% |
| with a measured antibiogram | 0 | 0.0% |

**The group is not one species.** Composition of `scientific_name`:

| species | isolates |
|---|---|
| Clostridium perfringens | 5,419 (95.0%) |
| [Clostridium] innocuum | 277 (4.9%) |
| Sarcina ventriculi | 9 (0.2%) |
| Clostridium mediterraneense | 1 (0.0%) |

Filter on `scientific_name` before reporting anything as a species-level figure.

## Resistance overview

Share of isolates carrying at least one gene of each class, and the genes behind it.

| class | isolates | main genes |
|---|---|---|
| TETRACYCLINE | 80.0% | `tetA(P)` (4,386), `tetB(P)` (3,063), `tet(44)` (723), `tet(M)` (72) |
| LINCOSAMIDE | 61.3% | `erm(Q)` (1,533), `cplR` (1,456), `lnu(P)` (1,111), `erm(B)` (159) |
| STREPTOGRAMIN | 52.1% | `erm(Q)` (1,533), `cplR` (1,456), `erm(B)` (159), `vat` (50) |
| MACROLIDE | 30.5% | `erm(Q)` (1,533), `erm(B)` (159), `mph(B)` (38), `erm(T)` (17) |
| PLEUROMUTILIN | 25.7% | `cplR` (1,456), `cfr(C)` (16), `cfr(B)` (2) |
| AMINOGLYCOSIDE | 16.0% | `ant(6)-Ib` (685), `aac(6')-Ie/aph(2'')-Ia` (206), `spw` (60), `ant(6)-Ia` (46) |
| PHENICOL | 3.9% | `optrA` (171), `fexA` (170), `cfr(C)` (16), `catA` (10) |
| OXAZOLIDINONE | 3.2% | `optrA` (171), `cfr(C)` (16), `cfr(B)` (2) |
| GLYCOPEPTIDE | 2.6% | `vanR` (148), `vanB` (1), `vanH-B` (1), `vanR-B` (1) |
| TRIMETHOPRIM | 0.3% | `dfrF` (16), `dfrA12` (1) |
| BETA-LACTAM | 0.2% | `blaTEM` (4), `blaTEM-116` (3), `blaCTX-M-55` (1), `bla` (1) |
| STREPTOTHRICIN | 0.2% | `sat4` (10) |

## Organism-specific settings (`scripts/groups.json`)

| setting | value | meaning |
|---|---|---|
| `pdg` | `PDG000000061.366` | Pathogen Detection release, prefix of the file names |
| `taxgroup` | `Clostridium_perfringens` | `whitelisted_taxa` used to pick point-mutation rows |
| `key_subclasses` | `TETRACYCLINE`, `ERYTHROMYCIN`, `LINCOSAMIDE`, `OXAZOLIDINONE`, `VANCOMYCIN`, `STREPTOMYCIN` | get their own 0/1 column |
| `marker_subclass` | `OXAZOLIDINONE` | genes of top concern, listed per isolate in `oxazolidinone_genes` |
| `antibiogram_query` | `antibiogram[filter] AND Clostridium perfringens[Organism]` | BioSample search for measured susceptibility |

## Read this before interpreting

- **No point mutations are screened for this organism**, so mutational resistance (for example in *gyrA*) is invisible here. That is a blind spot of the dataset, not an absence of resistance.
- **Family-level calls**: 1,747 of 14,709 AMR calls are profile (`HMM`) hits, which name a gene family rather than an allele. They can flag a housekeeping homologue as if it were a resistance gene, so treat a family-level call as a lead, not a result.
- **No measured antibiograms exist** for this group in NCBI BioSample, so nothing here can be checked against a laboratory result. The profiles are gene content only.
- **Gene absence is not susceptibility.** 241 isolates have no public assembly, and every call depends on what the pinned catalog knows.

## Files

```
AMR/ Metadata/ Clusters/ SNP_trees/ Exceptions/   original NCBI download - never modified
AMR_annotation/    gene -> antibiotic class, per-isolate profiles, per-cluster summaries
Genomes/           assemblies + proteins + annotation, and the manifest of what was downloaded
AMR_sequences/     protein / DNA sequences per gene family, for alignment
Antibiograms/      measured susceptibility from NCBI BioSample - searched, none exists for this group
```

Column-by-column documentation of every table: [`../docs/DATA_DOCUMENTATION.md`](../docs/DATA_DOCUMENTATION.md). How the pipeline was built and how to add another organism: [`../docs/WORKPLAN_new_species.md`](../docs/WORKPLAN_new_species.md).

## Rebuilding

```bash
python3 scripts/run_all.py --group Clostridium_perfringens                 # steps 01-05, the tables
python3 scripts/06_genome_manifest.py --group Clostridium_perfringens      # decide which accession per isolate
python3 scripts/07_download_genomes.py --group Clostridium_perfringens     # long, resumable
python3 scripts/08_extract_amr_sequences.py --group Clostridium_perfringens
python3 scripts/09_fetch_antibiograms.py --group Clostridium_perfringens
python3 scripts/10_group_readme.py --group Clostridium_perfringens         # this file
```
