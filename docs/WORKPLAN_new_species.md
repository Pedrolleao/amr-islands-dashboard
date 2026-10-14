# Workplan – repeating this analysis for another bacterial species

> **Steps 01-10 only.** Islands are steps 11-22, cluster analyses 23-26, domain annotation
> 27-28. For what exists per organism group and what is missing, see
> [`PROJECT_STATUS.md`](PROJECT_STATUS.md).

*Written 19 September 2026 from the Citrobacter freundii run (PDG000000039.1314) and validated the same day by
repeating everything for Citrobacter portucalensis (PDG000000147.9) – see section 12.*
What was done, in order, why it was done that way, what to change for a new organism, and the traps we hit.
The data produced by the Citrobacter run is described in `DATA_DOCUMENTATION.md`.

---

## 0. Overview

Goal: starting from an NCBI Pathogen Detection organism group, end up with
(1) per-isolate resistance profiles, (2) resistance linked to SNP clusters, (3) gene / protein-level data
(sequences), plus (4) whatever measured antibiograms exist.

```
 A. Pathogen Detection download ──┐
 B. Reference Gene Catalog ───────┼─> C. gene -> antibiotic class for ALL isolates ──> E. profiles + clusters
 D. MicroBIGG-E export ───────────┘         (D adds location, protein, copies for assembled isolates)
                                   F. genomes + proteins (NCBI datasets) ──> G. gene / protein sequences
                                   H. antibiograms (NCBI BioSample)
```

| phase | script | time (Citrobacter, 7,287 isolates) | disk |
|---|---|---|---|
| A, D | manual downloads | minutes | 0.2 GB |
| B, C, E | `01`–`05` (`run_all.py`) | < 1 min | 0.12 GB |
| F | `06`, `07` | under 1 h for 6,063 genomes | 59 GB (≈10 MB per genome: FASTA + proteins + GFF) |
| G | `08` | ≈ 20 min | 0.4 GB |
| H | `09` | 40 s | 1 MB |

Requirements: Python 3 + pandas, internet, the NCBI `datasets` program (single file, no installation), enough
disk for the genomes. No NCBI account or API key was needed.

---

## 1. Before starting – decisions to take

1. **Which organism group?** Pathogen Detection groups are not always one species. "Citrobacter_freundii"
   contains the whole genus (71 % *C. freundii*). Check `scientific_name` early and decide whether to keep all.
   List of groups: https://ftp.ncbi.nlm.nih.gov/pathogen/Results/
2. **How big is it?** Look at the line count of the metadata TSV first.
   - up to ~20,000 isolates: this workplan works as is.
   - *Klebsiella*, *E. coli/Shigella*, *Salmonella*, … (10⁵–10⁶ isolates): tables still work but need more
     memory; **do not download all genomes** (terabytes) – choose a subset (species, country, clusters, carriers
     of a gene) by filtering `genome_manifest.tsv` before step F; get MicroBIGG-E through Google BigQuery
     (`ncbi-pathogen-detect.pdbrowser.microbigge`) instead of the web export.
3. **Which questions?** Profiles only → phases A–C, E are enough (minutes). Gene location / plasmid linkage /
   variants → add D. Sequences, alignments, anything genome-based → add F, G.
4. **Freeze the versions.** Note the PDG release and the AMRFinderPlus database version (the scripts store
   the latter in every table). NCBI updates both continuously; mixing releases causes small mismatches.
   Example: we used release 1314 (16 Sep 2026); three days later `latest_snps/` on the FTP site already
   pointed to 1317. Old releases stay available under their own `PDG….N/` folder.

---

## 2. Phase A – Pathogen Detection download (manual)

From `https://ftp.ncbi.nlm.nih.gov/pathogen/Results/<Organism_group>/latest_snps/` (or a fixed
`PDGxxxxxxxxx.N/` folder) download and keep the structure:

```
<Organism_group>/
  AMR/<PDG>.amr.metadata.tsv            <- the essential file
  Metadata/<PDG>.metadata.tsv
  Clusters/<PDG>.reference_target.all_isolates.tsv, cluster_list.tsv, SNP_distances.tsv, new_isolates.tsv
  SNP_trees/*.tar.gz                    (optional, only for tree work)
  Exceptions/<PDG>.exceptions.tsv
```

Put it under `latest_snps/<Organism_group>/` and never edit these files.

First look (what we did): list the columns, count filled values of `AMR_genotypes`, `AST_phenotypes`,
`asm_acc`, species composition, number of clustered isolates. Findings that shaped the rest:
- gene **symbols** are present for every isolate, but **no antibiotic** and no protein information;
- `AST_phenotypes` (measured) exists for ~1 % of isolates;
- `asm_acc` tells which isolates have a public genome (83 %).

## 3. Phase B/C – gene → antibiotic class (scripts 01–03)

Source: AMRFinderPlus database, `https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest/`
→ `ReferenceGeneCatalog.txt`, `fam.tsv`, `AMRProt.fa`, `version.txt` (script `01`, pinned in
`reference/AMRFinderPlus_<version>/`).

Why the catalog first: it covers **all** isolates (also those without an assembly), is tiny, needs no account,
and reuses NCBI's own AMRFinderPlus calls.

Matching rules that were needed (script `02`), in this order:
1. symbol = `allele` (for point mutations keep the row whose `whitelisted_taxa` is the organism group);
2. symbol = `gene_family` → take class / subclass from the `fam.tsv` node (what AMRFinderPlus reports for a
   family-level call); flag `subclass_ambiguous` when the alleles of the family disagree;
3. symbol = a `fam.tsv` node or gene symbol (broad families: `tet`, `erm`, `aph(3')-I`…), top-most node;
4. last resort: the class MicroBIGG-E reports (needed for 1 symbol renamed / absent in the current catalog);
5. otherwise list in `unmatched_symbols.tsv` – **look at this file for a new species**.

Script `03` explodes the three genotype columns into `isolate_gene_long.tsv` and parses `AST_phenotypes`.

Checks to repeat: rows per isolate = NCBI's `number_*_genes`; % of occurrences annotated (we had 100 %).

## 4. Phase D – MicroBIGG-E (manual export + script 04)

MicroBIGG-E (Microbial Browser for Identification of Genetic and Genomic Elements) has one row **per gene
copy per assembled isolate**.

**Where the tables come from – how our two exports were made**

1. Open the MicroBIGG-E browser with the organism as search term (the text after `#` is the search):
   - https://www.ncbi.nlm.nih.gov/pathogens/microbigge/#Citrobacter%20portucalensis
   - for another organism replace the text after `#` (spaces as `%20`), or type it in the search box.
   - The *C. freundii* export (`microbigge-C_freundii.tsv`) came from the same browser; its exact search term was
     not written down. It contains precisely the 6,063 assembled isolates of the *C. freundii* group, so it was
     a narrower search than the *C. portucalensis* one (see the note on free-text search below).

   Link exactly as used for *C. portucalensis* (it came from an NCBI Insights post, hence the tracking
   parameters, which can be dropped):
   `https://www.ncbi.nlm.nih.gov/pathogens/microbigge/%5d?utm_source=ncbi_insights&utm_medium=referral&utm_campaign=pathogen-detection-ast-20240501#Citrobacter%20portucalensis`
2. Wait for the result table, then use **Download** (top of the table) → all rows, tab-separated (TSV).
3. Save as `Internship/<Organism_group>/microbigge-<short name>.tsv` and put that path in `scripts/groups.json`
   (`microbigge_tsv`). Note the search term and the download date (both Citrobacter files are dated 19 Sep 2026).

Columns the scripts need (the browser's default column set, names as in our exports):
`Scientific name, Protein, BioSample, Isolate, Contig, Start, Stop, Strand, Element symbol, Element name, Type,
Scope, Subtype, Class, Subclass, Method, % Coverage of reference, % Identity to reference`.
If columns were hidden in the browser ("Choose columns"), switch them back on before downloading.

**The search term is free text, not a strict species filter.** `Citrobacter portucalensis` returned 270,625 rows from
6,513 isolates of the *whole genus* (only 430 isolates were *C. portucalensis*). This is harmless – the loader
keeps only isolates that are in the group's Pathogen Detection metadata and prints how many rows were kept –
but check that printed line: if far fewer isolates are kept than there are isolates with `asm_acc`, the
search was too narrow. (A fielded search such as `taxgroup_name:"Citrobacter portucalensis"` should give a
smaller file; not tested here.) A genus-wide export can be reused for every group of that genus.

For very large groups (*Salmonella*, *E. coli*, *Klebsiella*) the browser download is impractical: use Google
BigQuery, table `ncbi-pathogen-detect.pdbrowser.microbigge`.

What it adds: location (contig, coordinates), the isolate's own protein accession, identity / coverage to
the reference, detection method, **copy number**, co-location of genes on a contig.
What it does not cover: isolates without a public assembly (17 % here).

Checks to repeat: isolate set = isolates with `asm_acc`; gene *sets* identical to the metadata
(MicroBIGG-E has more rows because it lists every copy); class agreement with the catalog
(`qc_class_differences.tsv`).

## 5. Phase E – profiles and clusters (script 05)

`isolate_resistance_profile.tsv` (one row per isolate: metadata, cluster, 0/1 per class, key subclasses,
carbapenemase genes, copy-level numbers) and `cluster_amr_summary.tsv` (one row per SNP cluster).

Adapt per species: `KEY_SUBCLASSES` and the "carbapenemase" logic are what mattered for an
Enterobacterales; for e.g. *Staphylococcus* you would flag METHICILLIN / VANCOMYCIN / LINEZOLID instead.
Remember intrinsic genes (chromosomal AmpC `blaCMY` in 98 % of *C. freundii*): every species has some, and
they inflate class percentages.

## 6. Phase F – genomes and proteins (scripts 06–07)

Idea (better than fetching proteins one by one): download each assembly **with its protein FASTA and GFF**
using NCBI `datasets`, then look the proteins of interest up locally.

- Tool: `https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/datasets` → `tools/`.
- `06`: asks NCBI for the assembly summary of every `asm_acc` and decides **GCA or GCF per isolate**:
  MicroBIGG-E reports some isolates on the RefSeq copy (`NZ_` contigs, `WP_` proteins) and others on the
  GenBank copy. Download the copy MicroBIGG-E used, otherwise accessions and contig names will not match.
- `07`: dehydrated packages in batches of 500 + `datasets rehydrate` (NCBI's recommended route for
  thousands of genomes); resumable; second pass retries failures; result in `Genomes/assemblies/<acc>/`.
- Test with 3 accessions first and confirm that MicroBIGG-E proteins and contigs are found in the files
  (we did: 100 %).

Checks: all accessions present, no empty files, `protein.faa` present wherever an annotation was expected.

## 7. Phase G – sequences of the genes of interest (script 08)

Per hit: protein by accession from the isolate's `protein.faa`; DNA cut from the genome by
contig + start/stop/strand (covers genes without a protein accession). Output per gene family:
`proteins/`, `proteins_unique/` (identical sequences collapsed by md5), `nucleotide/`, plus
`sequence_index.tsv` and `protein_variants.tsv`.

Why md5: GenBank protein accessions are unique per assembly, only `WP_` accessions are shared between
identical proteins – so identical proteins can only be recognised from the sequence itself.

Checks: DNA length = stop − start + 1 (100 %); proteins found (99.99 %); translation of the DNA equals the
annotated protein (98.8 %; the rest are contig-end fragments and alternative start codons).

## 8. Phase H – antibiograms (script 09)

NCBI keeps measured AST only as tables attached to BioSamples. Script `09` searches (a) every BioSample of
the dataset by accession with `antibiogram[filter]` and (b) `antibiogram[filter] AND <Genus>[Organism]`,
downloads the records, and writes MIC / disk values, method, platform, standard.
Expect little data (1 % here). More would have to come from the publications behind the big BioProjects
(`bioproject_acc` column).

---

## 9. What to edit for a new species

The scripts are species-independent. Everything species-specific is one entry in **`scripts/groups.json`**:

```json
"Citrobacter_portucalensis": {
  "pdg": "PDG000000147.9",
  "taxgroup": "Citrobacter_freundii",
  "microbigge_tsv": "../Citrobacter_portucalensis/microbigge-C_portucalensis.tsv",
  "antibiogram_query": "antibiogram[filter] AND Citrobacter portucalensis[Organism]"
}
```

| key | meaning |
|---|---|
| entry name | organism group = folder `latest_snps/<name>/` holding the NCBI download |
| `pdg` | release id, the prefix of the NCBI file names |
| `taxgroup` | value of `whitelisted_taxa` in `ReferenceGeneCatalog.txt` used to pick the right point-mutation rows (closest organism that has a set; irrelevant if the group has no `POINT` calls) |
| `microbigge_tsv` | path of the MicroBIGG-E export, relative to `latest_snps/`. It may be broader than the group (whole genus): rows of other isolates are dropped automatically |
| `antibiogram_query` | BioSample search used in addition to the per-accession search |

What to report is organism-specific and is set in the same entry (omit a key to get the Enterobacterales
default used for the Citrobacter groups):

| key | meaning |
|---|---|
| `key_subclasses` | subclasses that get their own 0/1 column in the profile, besides the class columns. Pick them from the data, not from habit: run steps 02-04 first, then look at which subclasses actually occur (`isolate_gene_long.tsv`) |
| `marker_subclass` | the subclass singled out as its own per-isolate gene list, the way carbapenemases are for Enterobacterales. Use `null` when the organism has no such acquired determinant: an invented marker is worse than none. *C. botulinum* is the worked example - its apparent linezolid signal is family-level `cfr` profile hits, so it has no marker |
| `marker_column` | name of that gene list column (`carbapenemase_genes`, `oxazolidinone_genes`, …) |
| `spot_checks` | `symbol -> [class, subclass]` assertions that fail the run loudly if the catalog says something else. Genes absent from a group are skipped. Take the expected values from the catalog, not from memory: `erm(B)` is `LINCOSAMIDE/MACROLIDE/STREPTOGRAMIN`, not `MACROLIDE` |

Choosing these per organism is the one judgement call the pipeline cannot make. For the Clostridia:
*C. difficile* uses `QUINOLONE` as its marker, because the `gyrA` / `gyrB` substitutions split the dataset
and track the epidemic lineages, while vancomycin and metronidazole look universal (see intrinsic genes
below); *C. perfringens* uses `OXAZOLIDINONE`, i.e. transferable `optrA`.

`reference/` and `tools/` are shared between species.

Run order for a new species (`G` = entry name):
```
# download phase A into latest_snps/G/   (wget -r -np -nH --cut-dirs=4 -R "index.html*" <FTP release folder>/)
python3 scripts/run_all.py --group G                 # 01-05  (needs phases A and D in place)
python3 scripts/06_genome_manifest.py --group G      # check the printed numbers, subset here if the group is large
python3 scripts/07_download_genomes.py --group G     # long; can be interrupted and re-run
python3 scripts/08_extract_amr_sequences.py --group G
python3 scripts/09_fetch_antibiograms.py --group G
```
(`export AMR_GROUP=G` once instead of repeating `--group G` also works.) After changing shared scripts, re-run
an existing group and compare checksums of its tables – we did this for *C. freundii*: byte-identical.

---

## 10. Traps we hit (and how they were solved)

| trap | symptom | solution |
|---|---|---|
| Metadata TSV quotes its list fields | ~380 "unknown" symbols starting with `"` | read that file with normal CSV quoting (`read_tsv(..., quoted=True)`); the catalog, `fam.tsv`, MicroBIGG-E and cluster files contain no quotes and are read literally (gene names contain `'`) |
| `NULL` strings | empty-looking values counted as data | converted to empty on reading |
| Metadata lists a symbol once, MicroBIGG-E once per copy | row counts differ by 32,000 | compare gene **sets**; copy number comes only from MicroBIGG-E |
| Family-level symbols (`blaCMY`, `tet`, `bla`) | not found as alleles | cascade allele → family → hierarchy (`fam.tsv`) |
| Same point mutation listed for many taxa | duplicated rows | keep the row of the organism group (`TAXGROUP`) |
| Isolates analysed with different database versions | same gene, different class between isolates | one pinned catalog for the harmonised class; original kept as `*_as_called` |
| A symbol missing from the current catalog (`dfrA51`) | unmatched | fall back to MicroBIGG-E class; report |
| Family-level subclass too generic (`blaOXA`) | carbapenemases missed or over-called | use the as-called subclass of the actual hit; count per (isolate, gene), never per gene name |
| GenBank vs RefSeq copies of an assembly | proteins / contigs "not found" | follow MicroBIGG-E: `NZ_`/`WP_` → download `GCF_`, otherwise `GCA_` |
| GenBank protein accessions unique per genome | "144,000 different proteins" | group by sequence (md5), not by accession |
| NCBI re-annotation after MicroBIGG-E was built | 19 proteins not found | DNA cut from coordinates still available |
| Windows drive (NTFS) under WSL | file names differing only by case or by `(` `'` would collide | sanitised, collision-checked file names |
| Long download in one command | risk of losing everything on interruption | batches + resumable script + log + `download_missing.txt` |
| Group ≠ species | genus-wide numbers reported as species numbers | always report / filter `scientific_name` |
| MicroBIGG-E export broader than the group | *C. portucalensis* export contained 6,513 isolates of the whole genus | loader keeps only isolates present in the group's metadata and prints how many rows were kept |
| Point mutations not screened for every organism | 0 `POINT` calls in *C. portucalensis* while half of *C. freundii* has `gyrA_T83I` | not an error: document it as a blind spot of that group |
| Intrinsic genes differ per species | ~100 % quinolone in *C. portucalensis* (chromosomal `qnrB`) | list the genes behind every class near 100 % before interpreting |

---

## 11. Checklist

- [ ] Organism group chosen, size checked, PDG release noted
- [ ] Phase A files in `latest_snps/<Organism_group>/` (untouched)
- [ ] MicroBIGG-E export saved, **filter and date written down**
- [ ] entry added to `scripts/groups.json`; `KEY_SUBCLASSES` / spot checks reviewed
- [ ] `run_all.py` finished; `qc_report.txt` read; `unmatched_symbols.tsv` empty or explained
- [ ] Gene sets identical between metadata and MicroBIGG-E
- [ ] Species composition and intrinsic genes of the species understood before interpreting percentages
- [ ] Genome subset decided (all, or filter) and disk space checked (≈10 MB per genome)
- [ ] Download complete (`download_missing.txt` empty), 3-genome lookup test passed
- [ ] Sequence extraction QC read (`extraction.log`)
- [ ] Antibiograms fetched; number of isolates with both genome and antibiogram noted
- [ ] READMEs / documentation updated with the new numbers

---

## 12. Second run: *Citrobacter portucalensis* (validation of this workplan)

Done on 19 Sep 2026 following the steps above; the whole run (download of the release, steps 01–09, 432 genomes) took about 10 minutes.

| | *C. freundii* group | *C. portucalensis* group |
|---|---|---|
| release | PDG000000039.1314 | PDG000000147.9 |
| isolates / with genome | 7,287 / 6,063 | 581 / 432 |
| SNP clusters | 925 | 68 |
| symbols unmatched | 0 | 0 |
| gene sets metadata = MicroBIGG-E | all | all (430 / 430) |
| genomes downloaded / size / time | 6,063 / 59 GB / < 1 h | 432 / 4.2 GB / few minutes |
| proteins found | 238,181 / 238,200 | 13,059 / 13,067 |
| distinct proteins / families | 8,533 / 325 | 1,484 / 204 |
| isolates with antibiogram (with genome) | 75 (57) | 5 (1) |
| carbapenemase carriers | 67 % (top: `blaKPC-2`) | 59 % (top: `blaNDM-1`) |

What the second run taught us (now built into the scripts): settings moved to `groups.json` + `--group`;
the MicroBIGG-E export can be genus-wide; spot-check genes may be absent; a group can have zero antibiograms or
zero point mutations; the file for records outside the dataset got a species-neutral name
(`antibiogram_outside_dataset.tsv`).
