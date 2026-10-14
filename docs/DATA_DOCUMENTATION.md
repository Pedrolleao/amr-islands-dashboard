# AMR annotation project – data documentation

> The dashboard over all of this is planned in
> [`WORKPLAN_dashboard.md`](WORKPLAN_dashboard.md).

> For what is finished and what is outstanding, see [`PROJECT_STATUS.md`](PROJECT_STATUS.md).

*Status: 21 September 2026. Describes everything under `Internship/latest_snps/`. Sections 1–9 use the
*C. freundii* group as the worked example; *C. portucalensis* (section 10) and the three Clostridia groups
(section 11) have the same structure. Each group folder also has its own generated `README.md`.*

This project takes the NCBI Pathogen Detection data for the *Citrobacter* group and adds, step by step:
which antibiotic each resistance gene acts on, per-isolate and per-cluster resistance profiles, the genomes
and proteins of every isolate with a public assembly, the sequences of all resistance genes, and the
lab-measured antibiograms that exist at NCBI.

---

## 1. The dataset in numbers

| | |
|---|---|
| Source | NCBI Pathogen Detection, organism group **Citrobacter_freundii**, release **PDG000000039.1314** (published 16 Sep 2026) |
| Isolates | **7,287** (`target_acc`, PDT…) from 7,276 BioSamples |
| Species | whole genus: *C. freundii* 5,205 · *C. koseri* 664 · *C. braakii* 483 · *C. amalonaticus* 293 · *C. farmeri* 183 · *C. werkmanii* 182 · *C. youngae* 113 · *C. sedlakii* 68 · others |
| In an SNP cluster | 3,890 isolates in **925** clusters (`PDS_acc`) |
| With a public genome assembly | **6,063** (`asm_acc` filled) – the other 1,224 exist only as raw reads (SRA) |
| Distinct gene symbols | 608 AMR · 71 stress · 100 virulence (all annotated) |
| Gene hits with location (MicroBIGG-E) | 256,528 in the 6,063 assembled isolates |
| Distinct resistance-related protein sequences | 8,533 in 325 gene families |
| Isolates with a lab antibiogram | **75** (57 of them with a genome) |
| AMRFinderPlus reference version used | **2026-08-07.1** |

Headline results (gene-based predictions, not measured resistance): 99.7 % of isolates carry a beta-lactam
gene, 79.5 % a quinolone determinant, 70 % aminoglycoside, 67 % a carbapenemase (`blaKPC-2` 1,526 isolates,
`blaNDM-1` 922, `blaOXA-48` 699, `blaKPC-3` 461, `blaVIM-1` 398), 7 % colistin (`mcr`). The database
over-represents resistant clinical isolates, so these are **not** population prevalences.

---

## 2. Folder map

```
Internship/
  microbigge-C_freundii.tsv          INPUT  MicroBIGG-E export (51 MB)
  Citrobacter_portucalensis/microbigge-C_portucalensis.tsv   INPUT  MicroBIGG-E export, genus-wide (53 MB)
  latest_snps/
    README.md                        short folder guide
    docs/                            this documentation + workplan for other species
    scripts/                         pipeline, numbered in running order (section 6)
    tools/                           NCBI datasets / dataformat programs (39 MB)
    reference/AMRFinderPlus_2026-08-07.1/   pinned NCBI reference (7 MB)
    Citrobacter_portucalensis/       second organism group, same sub-folders (section 10)
    Citrobacter_freundii/
      AMR/ Metadata/ Clusters/ SNP_trees/ Exceptions/ *.descriptor.xml    ORIGINAL NCBI download (never modified)
      AMR_annotation/     118 MB   gene -> antibiotic tables, profiles, cluster summaries
      Genomes/             59 GB   genome + proteins + annotation of 6,063 isolates
      AMR_sequences/      377 MB   sequences of the AMR / stress / virulence genes
      Antibiograms/         1 MB   lab-measured susceptibility (MIC values)
```

Every generated folder has its own `README.md` with file and column details; this document is the overview.

---

## 3. Identifiers – how everything links

| identifier | example | what it is | used to join |
|---|---|---|---|
| `target_acc` | PDT000004922.2 | one isolate in Pathogen Detection | **main key of every table** |
| `biosample_acc` | SAMN02603912 | the sample record (metadata, antibiogram) | antibiograms, NCBI lookups |
| `asm_acc` | GCA_000018045.1 | GenBank genome assembly | genomes (via `genome_manifest.tsv`) |
| `download_acc` | GCF_000018045.1 | the assembly copy actually downloaded (GCA or GCF) | folder name in `Genomes/assemblies/` |
| `PDS_acc` | PDS000050776.9 | SNP cluster | clusters, trees (`SNP_trees/<PDS_acc>.tar.gz`) |
| `symbol` | blaKPC-2 | gene / point-mutation name (AMRFinderPlus) | `symbol_annotation.tsv` |
| `gene_family` | blaKPC | family the symbol belongs to | FASTA file names in `AMR_sequences/` |
| `protein` | WP_004199234.1 | the isolate's own protein accession | `protein.faa`, FASTA headers |
| `protein_md5` / `variant_id` | blaKPC_736ad144 | fingerprint of a protein sequence | identical proteins across isolates |

Typical path: isolate (`target_acc`) → its genes (`isolate_gene_long.tsv`) → where they are and which
protein (`microbigge_hits.tsv`) → the sequence (`AMR_sequences/`) or the whole genome
(`genome_manifest.tsv` → `Genomes/assemblies/<download_acc>/`).

---

## 4. Original NCBI download (`Citrobacter_freundii/`, read-only)

| folder / file | content |
|---|---|
| `AMR/PDG…amr.metadata.tsv` (+ `.xml`) | **main input.** One row per isolate, 67 columns: sample metadata (species, date, country, host, source, BioProject, assembly, SRA run) + AMRFinderPlus results: `AMR_genotypes`, `AMR_genotypes_core`, `stress_genotypes`, `virulence_genotypes`, gene counts, `AST_phenotypes`, AMRFinderPlus / database version. The XML has the same content. |
| `Metadata/PDG…metadata.tsv` (+ `.xml`) | the same table published with the metadata step (near-identical; the AMR copy is the one used). |
| `Clusters/…all_isolates.tsv` | every isolate → its SNP cluster (`PDS_acc`, empty = not clustered), distance to nearest neighbour. |
| `Clusters/…cluster_list.tsv` | members of each cluster (3,890 rows). |
| `Clusters/…SNP_distances.tsv` | pairwise SNP distances within clusters (27,344 pairs). |
| `Clusters/…new_isolates.tsv` | isolates added in this release (6). |
| `SNP_trees/<PDS_acc>.tar.gz` (925) | per cluster: Newick tree, tree PDF, VCF of variable positions, parsimony input. No gene or protein sequences. |
| `Exceptions/…exceptions.tsv` | 1,533 samples NCBI excluded and why (species misidentified 579, missing strain info 548, low coverage 74, …). |
| `*.descriptor.xml` | release stamps (which parts are complete, publication date). |

Format traps: list fields in the metadata TSV are wrapped in double quotes; missing values are `NULL`;
genotype entries can carry a suffix (`=POINT`, `=PARTIAL`, `=PARTIAL_END_OF_CONTIG`, `=HMM`,
`=MISTRANSLATION`); each distinct symbol is listed once per isolate even when several copies exist.

---

## 5. Generated data

### 5.1 `reference/AMRFinderPlus_2026-08-07.1/`
`ReferenceGeneCatalog.txt` (11,456 rows: allele, gene family, product, **class, subclass**, reference protein
accession, literature) · `fam.tsv` (gene hierarchy, gives class/subclass for family-level names such as `tet`,
`blaCMY`) · `AMRProt.fa` (10,078 reference protein sequences) · `version.txt`.

### 5.2 `AMR_annotation/` – which gene, which antibiotic, which isolate
| file | rows | one row per | use it for |
|---|---|---|---|
| `symbol_annotation.tsv` | 779 | gene symbol | lookup symbol → `class`, `subclass`, product, reference protein, `match_level` |
| `isolate_gene_long.tsv` | 243,218 | isolate × element | **main table, all 7,287 isolates** |
| `isolate_resistance_profile.tsv` | 7,287 | isolate | metadata + cluster + 0/1 per antibiotic class + carbapenemase genes (join to anything) |
| `cluster_amr_summary.tsv` | 925 | SNP cluster | % of members per class, common genes, carbapenemases, years, countries |
| `microbigge_hits.tsv` | 256,528 | gene copy | contig, coordinates, method, % identity / coverage, protein accession, AMR genes on the same contig |
| `gene_variant_summary.tsv` | 5,220 | symbol × match signature | exact vs variant vs partial forms of each gene |
| `ast_long.tsv` | 1,122 | isolate × antibiotic | S/I/R summary from the metadata (75 isolates) |
| `qc_class_differences.tsv`, `unmatched_symbols.tsv`, `qc_report.txt` | | | quality control |

Key concepts
- **`class` / `subclass`**: drug class (BETA-LACTAM) and the specific drug group (CARBAPENEM). `A/B` = acts on
  both; profile columns count it for each. These are *predictions from gene content*.
- **`match_level`**: `allele` (exact named allele) · `gene_family` (only family known, e.g. `blaCMY`) ·
  `hierarchy_node`/`hierarchy_symbol` (broader family such as `tet`) · `microbigge_only` (`dfrA51`, absent from
  the current catalog).
- **`subclass_ambiguous = 1`**: alleles of the family differ in subclass (`blaOXA`: carbapenemase or not); the
  family-level value is generic, `subclass_all` lists the options, `subclass_as_called` in
  `microbigge_hits.tsv` is specific to the actual hit.
- **Harmonised vs as-called**: isolates were analysed by NCBI with ≥5 database versions. `class`/`subclass`
  use one pinned version for all isolates; `class_as_called`/`subclass_as_called` keep the original.
  Agreement: class 99.86 %, subclass 99.36 %.
- **`call_qualifier` / `hit_quality`**: complete/exact, variant, POINT, PARTIAL, PARTIAL_END_OF_CONTIG (cut by
  a contig edge – usually a real gene), HMM (family only), MISTRANSLATION / internal_stop (probably broken).
- **Copies**: the metadata lists a symbol once; MicroBIGG-E lists every copy (32,297 extra AMR rows). Gene
  *sets* agree for every shared isolate.

### 5.3 `Genomes/` – 6,063 assemblies, 59 GB
`assemblies/<download_acc>/` holds `<acc>_<name>_genomic.fna` (contigs), `protein.faa` (all proteins),
`genomic.gff` (annotation). `manifest/genome_manifest.tsv` maps `target_acc` → `download_acc` and records
assembly level (Contig 5,034 · Scaffold 648 · Complete 353 · Chromosome 28), size, contigs, and why GCA or GCF
was chosen: the copy MicroBIGG-E used (2,844 RefSeq `GCF_` with `NZ_` contigs / `WP_` proteins; 3,219 GenBank
`GCA_`), so coordinates and accessions can be looked up directly. Two assemblies have no annotation
(`GCA_902363265.2`, `GCA_986695725.1`): genome only.

### 5.4 `AMR_sequences/` – the resistance genes themselves
`proteins/<TYPE>/<gene_family>.faa` (every copy, 238,181 proteins) · `proteins_unique/…` (identical sequences
collapsed, 8,533) · `nucleotide/<TYPE>/<gene_family>.fna` (DNA of all 256,528 hits, including genes without a
protein accession) · `sequence_index.tsv` (one row per hit) · `protein_variants.tsv` (one row per distinct
protein: isolates, clusters, species). TYPE = AMR / STRESS / VIRULENCE. Headers:
`>isolate|protein accession|symbol|hit quality|contig:start-stop(strand)`.
19 proteins (mostly `aadA1/2`) were re-annotated by NCBI and have DNA only.

### 5.5 `Antibiograms/` – measured susceptibility
`antibiogram_full.tsv` (1,168 tests, 75 isolates): antibiotic, S/I/R, **measured value with sign and units**,
method (MIC 1,109 / disk diffusion 59), platform, vendor, standard (CLSI 1,158 / EUCAST 9).
`antibiogram_isolates.tsv` (per isolate; `has_genome`), `antibiogram_outside_dataset.tsv` (30 *Citrobacter*
BioSamples outside this dataset), raw XML. Verified exhaustive: every BioSample of the dataset was searched.

---

## 6. Scripts (`scripts/`)

| script | input → output |
|---|---|
| `groups.json` | **species-specific settings** (release id, MicroBIGG-E file, point-mutation taxon, antibiogram query), one entry per organism group |
| `amr_common.py` | paths, `--group <Organism_group>` selection (default `Citrobacter_freundii`), shared readers |
| `01_fetch_reference.py` | NCBI FTP → `reference/` |
| `02_symbol_annotation.py` | metadata + catalog + hierarchy → `symbol_annotation.tsv` |
| `03_isolate_gene_long.py` | metadata + lookup → `isolate_gene_long.tsv`, `ast_long.tsv` |
| `04_microbigge_hits.py` | MicroBIGG-E export → `microbigge_hits.tsv`, `gene_variant_summary.tsv` |
| `05_profiles_and_clusters.py` | → `isolate_resistance_profile.tsv`, `cluster_amr_summary.tsv` |
| `06_genome_manifest.py` | NCBI assembly summaries → `Genomes/manifest/` |
| `07_download_genomes.py` | NCBI datasets → `Genomes/assemblies/` (resumable, under 1 h here) |
| `08_extract_amr_sequences.py` | genomes + hits → `AMR_sequences/` (~20 min) |
| `09_fetch_antibiograms.py` | NCBI BioSample → `Antibiograms/` (~40 s) |
| `run_all.py` | runs 1–5 (`--sequences` adds 6–8) and writes `AMR_annotation/qc_report.txt` |

Requirements: Python 3 + pandas; internet; `tools/datasets` for steps 6–7.

---

## 7. Quality checks that passed

- Rows per isolate in `isolate_gene_long.tsv` equal NCBI's `number_amr_genes`, `number_stress_genes`,
  `number_virulence_genes` for all 7,287 isolates; 100 % of gene occurrences annotated.
- Gene sets per isolate identical between metadata and MicroBIGG-E (6,062 / 6,062 AMR, 6,063 / 6,063 stress,
  4,844 / 4,844 virulence).
- Spot checks: `blaKPC-2`, `blaNDM-1` → carbapenem; `sul1` → sulfonamide; `gyrA_T83I` → quinolone;
  `mcr-9.1` → colistin; `tet(A)` → tetracycline.
- Genomes: 6,063 / 6,063 downloaded, no empty files, protein file for 6,061 as predicted.
- Sequences: DNA cut for 256,528 / 256,528 hits with the exact expected length; proteins found for
  238,181 / 238,200; translated DNA equals the annotated protein for 98.8 % (differences = genes at contig
  ends or alternative start codons).
- Antibiograms: S/I/R identical to the metadata summary (1,122 combinations, none missing either way).

---

## 8. Limitations to keep in mind

1. **Predicted ≠ measured.** Class / subclass come from gene content. Only 75 isolates (1 %) have lab data.
2. **Sampling bias.** Public databases over-represent resistant, clinical, outbreak isolates and a few
   countries; do not read percentages as prevalence.
3. **Whole genus.** Filter on `scientific_name` for *C. freundii* only. Species labels are submitters' labels
   (NCBI excluded 579 samples as misidentified, see `Exceptions/`).
4. **1,224 isolates without a public assembly** have gene lists and classes but no location, copy number,
   protein or genome (raw reads exist in SRA: `Run` column).
5. **Draft assemblies.** 83 % are contig-level: genes can be split at contig ends
   (`partial_contig_end`), copy numbers can be collapsed, and "same contig" underestimates plasmid linkage.
6. **Family-level calls** (`match_level = gene_family`, `subclass_ambiguous = 1`) have a generic subclass.
7. **Intrinsic genes.** Chromosomal `blaCMY` (AmpC) is a normal part of the *C. freundii* genome: a `blaCMY`
   family gene is found in 98.2 % of *C. freundii* isolates (81.7 % of all isolates), and for 48 % of the
   "cephalosporin" isolates it is the only cephalosporin gene. So the high beta-lactam / cephalosporin figures
   mostly reflect the species, not acquired resistance; look at acquired genes (ESBLs, carbapenemases) instead.
8. **Snapshot.** Everything reflects release 1314 and catalog 2026-08-07.1; NCBI updates both continuously
   (19 proteins were already re-annotated between MicroBIGG-E and the genome download).

---

## 9. Not done yet (ideas)

- Annotate the SNP trees (`SNP_trees/`) with resistance classes / carbapenemases (iTOL or ggtree).
- Genotype vs phenotype comparison for the 75 isolates with antibiograms.
- Plots and statistics (by species, year, country, source).
- Alignments / phylogenies of gene families from `AMR_sequences/proteins_unique/` with `AMRProt.fa` references.
- Plasmid / mobile-element context of carbapenemases in the complete genomes (353).
- Assemble the 1,224 read-only isolates if they matter for a question.

---

## 10. Second organism group: *Citrobacter portucalensis* (PDG000000147.9)

Processed with the same pipeline (`--group Citrobacter_portucalensis`); identical folders, files and columns
under `latest_snps/Citrobacter_portucalensis/`. Details and caveats: that folder's `README.md`.

| | |
|---|---|
| Isolates | 581, all *C. portucalensis*; none shared with the *C. freundii* group (the two groups are disjoint, so tables can be concatenated) |
| In an SNP cluster | 212 in 68 clusters |
| With a genome (downloaded) | 432 (4.2 GB) – 430 in MicroBIGG-E |
| Gene symbols | 259 AMR · 59 stress · 27 virulence, all annotated |
| Gene hits / distinct proteins | 13,447 / 1,484 in 204 gene families |
| Antibiograms | 5 isolates (77 MIC tests), only 1 with a genome |
| Carbapenemase carriers | 343 (59 %): `blaNDM-1` 138, `blaKPC-2` 59, `blaKPC-3` 35, `blaNDM-5` 25, `blaOXA-48` 21, `blaIMP-4` 17 |

Differences that matter when comparing the two groups:
- **No point mutations are reported** for *C. portucalensis* (0 `POINT` calls), whereas `gyrA_T83I` is in half of the
  *C. freundii* group – a screening difference, not biology.
- **Intrinsic genes**: `blaCMY` in 579 / 581 and chromosomal `qnrB` in 576 / 581 isolates; hence ~100 % beta-lactam,
  cephalosporin *and quinolone* classes. Compare acquired genes only.
- The MicroBIGG-E input of this group is a genus-wide export; only this group's 13,447 rows are used.
- NDM dominates the carbapenemases here (KPC in the *C. freundii* group); both collections are biased towards
  clinical, resistant isolates (363 clinical, 73 environmental/other, 145 unstated).

---

## 11. Third genus: the Clostridia (`Clostridioides_difficile`, `Clostridium_perfringens`, `Clostridium_botulinum`)

Same pipeline, same folders and columns, selected with `--group`. These three groups are the reason the
reporting settings moved out of the code and into `scripts/groups.json` (section 9 of the workplan): what
counts as a resistance gene of interest is not the same for a *Clostridium* as for an Enterobacterales.

| | *C. difficile* | *C. perfringens* | *C. botulinum* |
|---|---|---|---|
| Release | PDG000000045.850 | PDG000000061.366 | PDG000000059.230 |
| Isolates | 39,708 | 5,706 | 1,676 |
| Dominant species in the group | *C. difficile* 98.9% | *C. perfringens* 95.0% | *C. botulinum* 70.5% |
| In an SNP cluster | 37,807 in 864 clusters | 3,535 in 815 clusters | 1,125 in 183 clusters |
| Genomes downloaded | 31,900 (228 GB) | 5,464 (32 GB) | 1,449 (9.6 GB) |
| Gene symbols (AMR / stress / virulence) | 186 / 8 / 8 | 87 / 12 / 26 | 27 / 5 / 38 |
| Gene hits with location | 667,433 | 28,821 | 8,066 |
| Protein copies extracted | 508,512 | 12,548 | 6,515 |
| Distinct proteins / gene families | 9,208 / 132 | 933 / 70 | 652 / 61 |
| Measured antibiograms | 105 isolates, 689 tests | none | none |
| `marker_subclass` | `QUINOLONE` | `OXAZOLIDINONE` | *none* |

### Why each group reports what it reports

- ***C. difficile*** is the only one of the three for which NCBI screens **point mutations**, so its profiles
  contain 56,837 `POINT` calls. The marker column is `fluoroquinolone_genes`, carried by 21,979 isolates
  (`gyrA_T82I` 8,035, plus `gyrB` substitutions): it splits the collection and tracks the epidemic lineages,
  which the treatment drugs cannot, because they look universal (next point).
- ***C. perfringens*** reports `oxazolidinone_genes`: transferable `optrA` in 171 isolates. Tetracycline
  determinants (`tetA(P)`, `tetB(P)`) are in 80% of isolates and are the background, not the signal.
- ***C. botulinum*** has **no marker column at all** (`marker_subclass: null`). Its apparent linezolid signal
  is 1,110 family-level `cfr` profile hits, which is a gene family assignment rather than a named allele;
  promoting that to a headline resistance figure would have been wrong. 2,250 of its 4,826 AMR calls are
  `HMM` hits, the highest share of the three.

### Traps specific to these groups

- **Intrinsic genes dominate *C. difficile*.** `nimB-Cd` (nitroimidazole) is in 39,498 of 39,708 isolates,
  the chromosomal `vanG` cluster in ~98.8%, and `blaCDD` / `blaAHM` give ~99% beta-lactam and ~99.6%
  "carbapenem". None of that is acquired resistance, and metronidazole and vancomycin remain the treatment
  drugs. Always list the genes behind a class above 90% before quoting the percentage: each group's
  `README.md` does this automatically.
- **One MicroBIGG-E export serves two groups.** The *C. botulinum* and *C. perfringens* exports are the same
  byte-identical genus-wide file (37,010 rows, 6 species); the loader keeps 8,066 rows for botulinum and
  28,821 for perfringens. That is expected, but check the "kept N of M rows" line.
- **Withdrawn assemblies.** 30 *C. difficile* and 1 *C. perfringens* accessions are in the metadata but are
  no longer served by NCBI `datasets` (all recent `GCA_0614…`). They are listed in
  `Genomes/manifest/download_missing.txt` and in the group README.
- **NCBI FTP refuses recursive `wget`** (`robots.txt`), which silently yields zero files. Fetch the directory
  listing, then download from an explicit file list.

### Antibiograms

Only *C. difficile* has measured data: 105 isolates, 689 tests (505 MIC, 184 disk diffusion). Erythromycin,
tetracycline, clindamycin and moxifloxacin show substantial resistance; vancomycin, metronidazole and
fidaxomicin are almost entirely susceptible. Searching NCBI BioSample exhaustively returns **one** antibiogram
for the whole genus *Clostridium*, and that isolate is not in any of these groups (kept in
`Clostridium_botulinum/Antibiograms/antibiogram_outside_dataset.tsv`). For those two groups the profiles
therefore rest on gene content alone, with nothing to validate them against.

---

## 12. Resistance islands (steps 11-22, `Citrobacter_freundii/Pangenome/` and `Islands/`)

Run so far for *C. freundii* only, with `scripts/run_islands.py --group Citrobacter_freundii`. The group
README and `Islands/README.md` carry the numbers; this section says what the pipeline is and where the
judgement calls sit.

**The question it answers.** Sections 5.2-5.4 say which isolate carries which resistance gene. The island
pipeline asks what each resistance gene is *sitting next to* in the assembly: the cassette array, the
transposon, the plasmid backbone around it.

**How it is built.**

| steps | what they produce |
|---|---|
| 11-12 | a gene index of every CDS in all 6,061 genomes, and every MicroBIGG-E hit mapped onto its coding sequence |
| 13-15 | protein clusters at 90% identity / 80% coverage, both directions (MMseqs2), with per-cluster annotation |
| 16 | the anchor set: clusters carrying resistance evidence, tiered core / shell / mobile |
| 17-19 | islands, the gene table behind them, and recurrent island types |
| 20-22 | per-island diagrams, the F01-F04 summary figures, and 17 QC checks |

**The three decisions that shape every number.**

1. **Anchoring is at cluster level, not hit level.** A gene anchors if its protein cluster contains a known
   resistance gene anywhere in the dataset, so the label carries over to copies the annotation missed. The
   gyrA check shows the effect: 6,055 anchors against ~3,202 hits.
2. **An island is +/-10 genes around each anchor**, merged where windows overlap, never spanning a contig.
3. **Stress and virulence never anchor.** They are coloured in the diagrams and annotated on the clusters,
   but an island is defined by resistance alone. This is why the stress and virulence work of section 13
   changes nothing here.

**What must be read before using the output.** All of it is in `Islands/README.md`; the two that bite
hardest are that **74% of islands run off a contig edge**, so size and composition claims need the clean set
(1,419 tier1 islands with a complete window), and that **a cluster, not a gene label, is the unit of "a
different gene"** - tiering on labels inflated the multi-resistance count 1.4-fold before it was fixed.
One known limit stands: at 90% identity the plasmid-borne and chromosomal `qnrB` alleles fall in a single
cluster, and no threshold separates them.

## 13. Per-cluster variation, structures and views (steps 23-25)

Three steps that take a **set of protein clusters** and describe each one in depth. They run per set:
`--set amr` (the 411 anchor clusters, the default), `--set stress` (328) or `--set virulence` (67).
`scripts/cluster_sets.py` defines the sets: `amr` reads `Islands/anchor_clusters.tsv`, the other two select
on `element_type` in `Pangenome/cluster_annotation.tsv`. Output goes to `Variation[_set]/` and
`Structures[_set]/`, one folder per set, each with its own README.

| step | what it does | output |
|---|---|---|
| 23 | MAFFT alignment of the distinct sequences, then per-position substitution / deletion / insertion shares against a reference | `Variation*/` summary, position table, alignments, one figure per cluster |
| 24 | ESMFold structure of each cluster's representative, locally, no MSA | `Structures*/pdb/` (B-factor = pLDDT) and `pdb_variation/` (B-factor = step-23 variation) |
| 25 | py3Dmol page plus PNG stills per structure, colourable by either B-factor | `Structures*/views/` |

**Two conventions make the numbers comparable.** Variation is counted **per lineage** (a Pathogen Detection
SNP cluster, or the isolate itself outside one) so a clonal outbreak counts once, with the per-copy share
drawn as a second line. And the reference is the most common **full-length** sequence, except where the
catalog defines point-mutation sites: there it is the most common sequence that is wild type at all of them,
because for gyrA the commonest sequence is the resistant I83 mutant.

**Where the three sets stand.**

| set | clusters | with >=2 sequences | structures | median mean pLDDT | below 70 |
|---|---|---|---|---|---|
| amr | 411 | 253 | 411 | 91.4 | 10 |
| stress | 328 | 237 | 328 | 85.8 | 79 |
| virulence | 67 | 40 | 61 folded, 6 deferred | 82.9 | 10 |

**Confidence falls with how well the family is represented in protein databases, not with the quality of the
run.** Beta-lactamases are among the best-covered families anywhere; small periplasmic metal-binding
proteins (`silE`, `pcoE`, `asr`: pLDDT 41-56) and giant NRPS modules are not. A low pLDDT on a disordered
protein is the honest answer. Treat anything below 70 as a lead, not a model.

**Two limits specific to the new sets.** `clbB` (3,206 aa) and `sslE` (1,505 aa) exceed what ESMFold fits on
a 16 GB card and are listed in `Structures_virulence/deferred_for_colabfold.{tsv,faa}` rather than dropped.
Four more went the same way on time: `sigA` (1,287 aa) held the GPU for a full hour without returning, and
`icmF`, `tsh` and `clbN` were dropped with it, being longer.
Step 24 evicts the ESM-2 language model from the card as soon as it has embedded the sequence: it runs once
per protein and the trunk never touches it, so its ~5.7 GB is dead weight through the recycling. This is
bit-identical, since the embeddings are still computed on the GPU with the same weights, and it is the
default (`--no-offload` restores the old behaviour). Measured here: silA (1,048 aa) 328 s -> 131 s at an
unchanged pLDDT of 88.30, sigA (1,287 aa) from over an hour with no result to 562 s, sslE (1,505 aa) 1,311 s.
Activation demand is about 7,200 bytes x length squared, so the no-spill ceiling is ~1,210 aa on this 16 GB
card with a desktop attached; beyond that CUDA on WSL spills into host memory over PCIe and slows down
3-30x without ever raising an error.

Step 24 also enforces a **time cap of one hour per protein** (`--max-seconds`, default 3600): a protein that
passes it is dropped, and so is every protein still queued, because the queue is sorted by length and the
rest are longer. Everything dropped, for either reason, is written to `Structures*/colabfold/` as its own
FASTA file and carries `status=deferred_too_long` or `deferred_too_slow` with a reason in
`structure_summary.tsv`. Fold those elsewhere with ColabFold, put the PDB in `colabfold/returned/`, and run
`scripts/26_ingest_colabfold.py`: it matches the file to its cluster by the `CFR-` id in the name, checks the
residue count against the representative sequence, writes both B-factor variants, and turns the deferred row
into a predicted one.
And cluster `CFR-016271` carries the virulence symbol `astA` although its proteins are 402-aa IS256
transposases: the 114 bp `astA` hit lies inside the transposase CDS and was attached at tier
`T4_overlap_partial`, which step 15 reports without filtering. It is the only cluster of any size in that
position; `Variation_virulence/README.md` documents it.

## 14. The dashboard (steps 30-37, `Citrobacter_freundii/Dashboard/`)

Three self-contained HTML pages over everything above, opened straight from the folder with no
server: `islands.html` (46,921 tier-1 and tier-2 islands with gene maps, a distribution tab and an
aligned view per type, cluster or AMR class), `clusters.html` (the 806 analysed clusters in full, each
with its Cα trace drawn on a canvas, every other cluster searchable) and `antibiograms.html` (the 75
isolates with measured susceptibility beside their genotype, step 37). `python3 scripts/run_dashboard.py` rebuilds
them in half a minute: 31-34 write the five JSON payloads (16.3 MB), 35 splices them and the shared
`scripts/dashboard/core.{css,js}` into the templates in `scripts/dashboard/`, and 36 decodes the
payloads again and compares them with the source tables and PDBs. The folder's `README.md` documents
the pages, the payloads and the codec; `docs/WORKPLAN_dashboard.md` holds the design, the field
schemas and the three packing bugs that the round-trip check exists to catch. Pages link to each
other through the address bar (`clusters.html?cluster=CFR-000010`, `islands.html?island=CFRI046929`).
