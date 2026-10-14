# Clostridium botulinum – AMR annotation of Pathogen Detection group PDG000000059.230

Built by `latest_snps/scripts/run_all.py --group Clostridium_botulinum`. Every number below was computed from the files in this folder; re-run `scripts/10_group_readme.py --group Clostridium_botulinum` after a rebuild.

Reference catalog pinned at AMRFinderPlus **2026-08-07.1**. MicroBIGG-E export: `microbigge-Clostridium_botulinum.tsv`.

## The dataset

| | count | share |
|---|---|---|
| isolates in the group | 1,676 | 100% |
| with a public assembly (`asm_acc`) | 1,449 | 86.5% |
| genomes downloaded here | 1,449 | 86.5% |
| covered by MicroBIGG-E (gene copies) | 1,430 | 85.3% |
| in an SNP cluster | 1,125 | 67.1% |
| SNP clusters | 183 | |
| carrying >=1 AMR gene | 1,307 | 78.0% |
| with a measured antibiogram | 0 | 0.0% |

**The group is not one species.** Composition of `scientific_name`:

| species | isolates |
|---|---|
| Clostridium botulinum | 1,440 (85.9%) |
| Clostridium sporogenes | 212 (12.6%) |
| Clostridium novyi | 14 (0.8%) |
| Clostridium tepidum | 8 (0.5%) |
| Clostridium uliginosum | 1 (0.1%) |
| 1 further taxon | 1 |

Filter on `scientific_name` before reporting anything as a species-level figure.

## Resistance overview

Share of isolates carrying at least one gene of each class, and the genes behind it.

| class | isolates | main genes |
|---|---|---|
| LINCOSAMIDE | 77.2% | `cplR` (1,125), `cfr` (1,110), `lsa` (639), `lnu(A)` (2) |
| STREPTOGRAMIN | 77.1% | `cplR` (1,125), `cfr` (1,110), `lsa` (639), `vat` (129) |
| PLEUROMUTILIN | 77.1% | `cplR` (1,125), `cfr` (1,110) |
| PHENICOL | 75.7% | `cfr` (1,110), `catA` (802) |
| OXAZOLIDINONE | 66.2% | `cfr` (1,110) |
| FOSFOMYCIN | 30.3% | `fosX` (506), `fosD` (1) |
| BETA-LACTAM | 11.9% | `blaR1` (144), `blaCBP-1` (73), `blaCBP` (71), `bla` (54) |
| AMINOGLYCOSIDE | 3.9% | `ant(9)` (60), `ant(6)-Ia` (4), `spw` (4), `aph(3')-IIa` (2) |
| MACROLIDE | 0.7% | `mef(A)` (11), `msr` (2), `erm(B)` (1) |
| GLYCOPEPTIDE | 0.5% | `vanR` (8) |
| TETRACYCLINE | 0.5% | `tet(M)` (6), `tet` (2) |

## Organism-specific settings (`scripts/groups.json`)

| setting | value | meaning |
|---|---|---|
| `pdg` | `PDG000000059.230` | Pathogen Detection release, prefix of the file names |
| `taxgroup` | `Clostridium_botulinum` | `whitelisted_taxa` used to pick point-mutation rows |
| `key_subclasses` | `CHLORAMPHENICOL`, `LINCOSAMIDE`, `LINEZOLID`, `FOSFOMYCIN`, `STREPTOGRAMIN`, `TETRACYCLINE` | get their own 0/1 column |
| `marker_subclass` | *(none)* | no gene of top concern in this organism, so no marker column |
| `antibiogram_query` | `antibiogram[filter] AND Clostridium botulinum[Organism]` | BioSample search for measured susceptibility |

## Read this before interpreting

- **No point mutations are screened for this organism**, so mutational resistance (for example in *gyrA*) is invisible here. That is a blind spot of the dataset, not an absence of resistance.
- **Family-level calls**: 2,250 of 4,826 AMR calls are profile (`HMM`) hits, which name a gene family rather than an allele. They can flag a housekeeping homologue as if it were a resistance gene, so treat a family-level call as a lead, not a result.
- **No measured antibiograms exist** for this group in NCBI BioSample, so nothing here can be checked against a laboratory result. The profiles are gene content only.
- **Gene absence is not susceptibility.** 227 isolates have no public assembly, and every call depends on what the pinned catalog knows.

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
python3 scripts/run_all.py --group Clostridium_botulinum                 # steps 01-05, the tables
python3 scripts/06_genome_manifest.py --group Clostridium_botulinum      # decide which accession per isolate
python3 scripts/07_download_genomes.py --group Clostridium_botulinum     # long, resumable
python3 scripts/08_extract_amr_sequences.py --group Clostridium_botulinum
python3 scripts/09_fetch_antibiograms.py --group Clostridium_botulinum
python3 scripts/10_group_readme.py --group Clostridium_botulinum         # this file
```
