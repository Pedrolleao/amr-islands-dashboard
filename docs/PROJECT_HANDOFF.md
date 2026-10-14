# Project handoff - AMR comparative genomics from NCBI Pathogen Detection

*Written 29 September 2026. One self-contained summary of the whole project, for someone (or an AI
assistant) who has none of the files. It covers what was done, in what order, why each choice was made,
what came out, and what to reuse in a new project. The detailed docs it condenses are listed at the end.*

---

## 1. The project in one paragraph

An internship analysis of antimicrobial resistance (AMR) in bacterial genomes from **NCBI Pathogen
Detection**. Starting from one organism group (*Citrobacter freundii* group, release PDG000000039.1314:
7,287 isolates, 6,063 with public genomes), a numbered Python pipeline (steps 01-37) goes from NCBI's gene
lists to (1) per-isolate resistance profiles, (2) resistance linked to SNP clusters (outbreak lineages),
(3) a gene- and protein-level study: pangenome clustering, resistance islands, per-position variation,
predicted 3D structures, Pfam domains and structural fold matches, and finally (4) a three-page offline
HTML dashboard over everything. Steps 01-10 were repeated for four more groups (*C. portucalensis*,
*C. difficile*, *C. perfringens*, *C. botulinum*) to prove the pipeline is organism-independent.

**Original goals (19 Sep 2026)** and status:

| goal | status |
|---|---|
| 1. Per-isolate resistance profiles | done, all five groups |
| 2. Link resistance to SNP clusters and trees | clusters done; putting resistance *onto the trees* not done |
| 3. Gene / protein-level study | done for *C. freundii* (islands, variation, structures, domains, folds, dashboard) |

---

## 2. Inputs (all public, no account or API key needed)

| input | where | what it gives |
|---|---|---|
| Pathogen Detection release | `ftp.ncbi.nlm.nih.gov/pathogen/Results/<Organism_group>/<PDG release>/` | `AMR/*.amr.metadata.tsv` (one row per isolate: metadata + AMRFinderPlus gene symbols, `AST_phenotypes`, `asm_acc`), `Clusters/` (SNP cluster membership), `SNP_trees/*.tar.gz` (per-cluster Newick + VCF), `Exceptions/` |
| AMRFinderPlus reference catalog | `ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest/` | `ReferenceGeneCatalog.txt` (gene -> class/subclass), `fam.tsv` (gene hierarchy), `AMRProt.fa` |
| MicroBIGG-E export | `ncbi.nlm.nih.gov/pathogens/microbigge/#<organism>` -> Download TSV (BigQuery `ncbi-pathogen-detect.pdbrowser.microbigge` for huge groups) | one row per gene *copy* per assembled isolate: contig, coordinates, strand, protein accession, % identity/coverage, method |
| Genomes | NCBI `datasets` CLI (single binary) | genome FASTA + `protein.faa` + GFF per assembly (~10 MB each; 59 GB for 6,063) |
| Antibiograms | NCBI BioSample, `antibiogram[filter]` | measured MIC / disk values (only ~1% of isolates have any) |

**Key identifiers:** `target_acc` (PDT..., one isolate - the main key of every table), `biosample_acc`,
`asm_acc` (GCA_), `download_acc` (GCA_ or GCF_ actually downloaded), `PDS_acc` (SNP cluster),
`symbol` (gene name), `protein_md5` (identity of a protein sequence).

---

## 3. The pipeline, step by step

Python 3 + pandas throughout; every script takes `--group <Organism_group>`; all organism-specific
settings live in one `groups.json` entry. Runners: `run_all.py` (01-10), `run_islands.py` (11-22),
`run_dashboard.py` (30-37); steps 23-29 run per cluster set with `--set amr|stress|virulence`.

### Phase 1 - Annotation and profiles (steps 01-05, under 1 minute)

| step | what | method / key rule |
|---|---|---|
| 01 | Fetch and **pin** the AMRFinderPlus catalog | stored under `reference/AMRFinderPlus_<version>/`; version written into every table |
| 02 | Map every gene symbol to antibiotic class/subclass | cascade: exact allele -> gene family (via `fam.tsv`) -> broader hierarchy node -> MicroBIGG-E's class as last resort -> `unmatched_symbols.tsv`. Point mutations: keep the catalog row whose `whitelisted_taxa` matches the organism. Flag `subclass_ambiguous` when a family's alleles disagree (e.g. `blaOXA`: carbapenemase or not) |
| 03 | Explode genotypes to a long table (isolate x gene), parse `AST_phenotypes` | check: rows per isolate = NCBI's `number_*_genes` |
| 04 | Load MicroBIGG-E | adds location, protein accession, copy number; keep only isolates in the group's metadata (exports are often genus-wide) |
| 05 | Per-isolate profile (0/1 per class + key subclasses + marker gene list) and per-SNP-cluster summary | `key_subclasses` and `marker_subclass` chosen per organism *from the data* (carbapenemases for Enterobacterales, fluoroquinolone for *C. difficile*, `optrA` for *C. perfringens*, none for *C. botulinum*) |

### Phase 2 - Genomes and sequences (steps 06-10)

| step | what | key rule |
|---|---|---|
| 06 | Genome manifest: choose GCA or GCF per isolate | download **the copy MicroBIGG-E used** (`NZ_`/`WP_` -> GCF, else GCA), otherwise accessions and contigs will not match |
| 07 | Download genomes | dehydrated packages in batches of 500 + `datasets rehydrate`; resumable; retry pass; test on 3 accessions first |
| 08 | Extract AMR / stress / virulence sequences per gene family | protein by accession from `protein.faa`; DNA cut by contig+coordinates; identical proteins collapsed **by md5 of the sequence** (GenBank accessions are unique per assembly) |
| 09 | Fetch antibiograms from BioSample | search every BioSample of the dataset + a genus query |
| 10 | Write a README per group from the tables that exist | automatically lists the genes behind any class above 90% (intrinsic-gene warning) |

### Phase 3 - Pangenome and resistance islands (steps 11-22, *C. freundii* only)

| step | what |
|---|---|
| 11-12 | Index every CDS of all genomes; map each MicroBIGG-E hit onto its CDS (tiered matching; weak overlaps kept but flagged, never propagate labels) |
| 13-15 | **Protein clustering with MMseqs2 at 90% identity / 80% coverage, both directions** -> 208,757 clusters; per-cluster annotation (product, prevalence, element type AMR/STRESS/VIRULENCE, core/shell/cloud) |
| 16 | **Anchors**: clusters containing any known resistance gene (409 clusters), tiered CORE / SOFTCORE / SHELL / MOBILE by prevalence (threshold calibrated with a prevalence sweep) |
| 17 | **Islands**: union of +/-10-gene windows around every anchor on a contig, overlaps merged, never crossing a contig |
| 18-19 | Island gene table; recurrent **island types** (same ordered arrangement of informative clusters) |
| 20-21 | Per-island gene-map figures; summary figures F01-F04 |
| 22 | 17 QC checks (e.g. gyrA propagation, known controls land in the right tier, truncation rate recorded) |

Island tiers: `tier1_multi_acquired` (>=2 distinct acquired resistance clusters), `tier2_single_acquired`,
`tier3_core_only` (backbone control - later shown as "chromosomal window", never "island").

### Phase 4 - Deep per-cluster analysis (steps 23-29), run for three sets

Sets: **amr** = 411 anchor clusters, **stress** = 328 (metal/biocide/acid), **virulence** = 67.

| step | what | method |
|---|---|---|
| 23 | Per-position variation | MAFFT alignment of distinct sequences; substitution / deletion / insertion share per position against a reference; **counted per lineage** (SNP cluster, or isolate if unclustered) so clonal outbreaks count once. Reference = most common full-length sequence, except where the catalog defines point-mutation sites: then the most common *wild-type* sequence (for gyrA the commonest sequence is the resistant T83I mutant) |
| 24 | 3D structure of each cluster representative | **ESMFold** locally on a 16 GB GPU (protocol in section 6); output PDBs with B-factor = pLDDT and a second copy with B-factor = variation |
| 25 | Interactive views | py3Dmol HTML page + PNG stills per structure |
| 26 | Ingest external ColabFold models | for proteins too long for local folding; matched by cluster id, residue count checked |
| 27 | Pfam domains | `hmmsearch` (faster than `hmmscan`) with Pfam 38.0 gathering thresholds + clan competition; also run over the **whole pangenome** (208,757 reps, 11.6 h) |
| 28 | Fold matches | **Foldseek** against PDB + AlphaFold Swiss-Prot; same fold = TM-score >= 0.5; model pLDDT carried alongside |
| 29 | Fold a giant modular protein by domain | groups Pfam domains into modules under a max length, cuts in linkers, 20 aa overlap (used for colibactin `clbB`, 3,206 aa) |

### Phase 5 - Dashboard (steps 30-37)

Three self-contained offline HTML pages (no server, no JS libraries), data packed into inline JSON:

- **islands.html** - every island with an SVG gene map, filters (species, lineage, AMR class, tier),
  architecture types, and *aligned views* of all islands of a type / cluster / AMR class with a
  conservation table.
- **clusters.html** - one record per analysed cluster: identity, anchor tier, variation profile, Pfam
  architecture, islands, and a canvas 3D Cα-trace viewer with confidence/variation colouring and ranked
  Foldseek hits; all 208,757 clusters searchable.
- **antibiograms.html** - measured susceptibility of the 75 tested isolates beside their genotype
  (carriers vs non-carriers per antibiotic and gene family).

Steps 30-34, 37 write payloads with a custom codec (base64 varint + zigzag + delta coding, ~16 MB total);
step 35 splices them into templates; **step 36 decodes every payload again and compares it with the
source tables** (the only thing that caught three silent encoding bugs).

---

## 4. Final output - the numbers (*C. freundii* group)

| item | value |
|---|---|
| isolates / with genome / in SNP clusters | 7,287 / 6,063 / 3,890 in 925 clusters |
| species in the "C. freundii" group | whole genus: *C. freundii* 71%, *C. koseri*, *C. braakii*, ... (17 species) |
| gene symbols annotated | 608 AMR, 71 stress, 100 virulence - 100% matched |
| gene copies located (MicroBIGG-E) | 256,528 |
| headline genotypes (predicted, biased sample) | 99.7% beta-lactam, 79.5% quinolone, 67% carbapenemase (`blaKPC-2` 1,526, `blaNDM-1` 922, `blaOXA-48` 699), 7% `mcr` |
| protein clusters | 208,757; 409 resistance anchors |
| islands | 66,005 (17,665 tier 1 / 29,256 tier 2 / 19,084 tier 3); 12,338 recurrent types; 74% touch a contig edge; clean set = 1,419 complete tier-1 islands |
| structures | 806/806 clusters folded (median pLDDT: amr 91.4, stress 85.8, virulence 82.9); 5 via ColabFold, `clbB` in 3 modules |
| Pfam named | amr 401/411, stress 309/328, virulence 60/67; whole pangenome 60.8% (99.4% of core vs 59.9% of cloud) |
| fold matched | amr 410, stress 264, virulence 63 |
| antibiograms | 75 isolates, 1,168 tests (MIC mostly, CLSI) |
| dashboard | 3 pages, rebuilt in ~1 min, verified by round-trip decoding + screenshots |

**Folder layout per group:** original NCBI download (never modified) + `AMR_annotation/`, `Genomes/`,
`AMR_sequences/`, `Antibiograms/`, `Pangenome/`, `Islands/`, `Variation[_set]/`, `Structures[_set]/`,
`Domains[_set]/`, `Dashboard/` - each with its own README of files and columns.

---

## 5. Design decisions that shape every number

1. **Pin database versions** (AMRFinderPlus catalog, Pfam, AlphaFold DB); keep NCBI's original call in
   `*_as_called` columns next to the harmonised one.
2. **Count per lineage, not per genome** - public databases are full of outbreak clones.
3. **Anchoring at cluster level, not hit level** - the resistance label carries to copies the annotation
   missed (gyrA: 6,055 anchors vs ~3,202 hits).
4. **A protein cluster, not a gene label, is the unit of "a different gene"** - tiering on labels
   inflated multi-resistance 1.4x.
5. **Stress and virulence genes never define islands** - they are annotated but only resistance anchors.
6. **Nothing is silently missing** - every skipped/deferred item has a row and a reason.
7. **Predicted is not measured, and sampling is biased** - percentages are never presented as prevalence.
8. **Intrinsic genes are background, not signal** - always list the genes behind any class above 90%
   (chromosomal AmpC `blaCMY` in 98% of *C. freundii*; `nimB`, `vanG` in *C. difficile*).

---

## 6. Traps already paid for (reuse these lessons)

**Data**
| trap | fix |
|---|---|
| Metadata TSV wraps list fields in double quotes; missing = `NULL` | parse with CSV quoting on; convert `NULL` to empty. Other NCBI files contain `'` in gene names and must be read *without* quoting |
| Metadata lists each symbol once, MicroBIGG-E once per copy | compare gene *sets*; copy number only from MicroBIGG-E |
| Family-level symbols (`blaCMY`, `tet`) | allele -> family -> hierarchy cascade |
| GenBank vs RefSeq copy of an assembly | download the copy MicroBIGG-E used |
| GenBank protein accessions are per-assembly | group identical proteins by sequence md5 |
| Organism group is not one species | report/filter `scientific_name` always |
| MicroBIGG-E free-text search returns the whole genus | loader filters to the group; check the "kept N of M rows" line |
| Point mutations screened only for some organisms | zero POINT calls is a blind spot, not biology |
| NCBI FTP refuses recursive `wget` (robots.txt) - silently gets nothing | fetch an explicit file list |
| Windows/NTFS drive under WSL | file names are case-insensitive: sanitise and collision-check |
| Long downloads | batches, resumable, logged, missing-list file |

**Structure prediction (ESMFold on a 16 GB consumer GPU under WSL)**
- Load, then `model.esm = model.esm.half()`, **then** `.cuda()` (reverse order fills the card).
- **Evict the ESM-2 language model after it embeds the sequence** - it is not used again during
  recycling; bit-identical results, ~5.7 GB freed, 2.5x faster on long proteins.
- Too little VRAM on WSL **never raises OOM** - CUDA silently spills to host RAM and runs 3-30x slower,
  looking like a hang. Memory ~7,200 bytes x length^2; ~1,200 aa is the no-spill ceiling on 16 GB.
- One-hour watchdog per protein (`os._exit` from a thread, since CUDA calls cannot be interrupted);
  queue sorted by length so everything after a timeout is deferred too.
- Deferred proteins go to one FASTA each for ColabFold, and come back through an ingest step.
- **Fold giant modular proteins by domain**, not as one chain: clbB in 3 modules scored higher (pLDDT
  80-89) than any whole-chain ColabFold model of shorter proteins.
- Low pLDDT on small disordered metal-binding proteins is the honest answer, not a failure; TM-score is
  the wrong test for them - read the hit's name before dismissing a weak score.

**Dashboard / code**
- Three encoding bugs were invisible to "no errors": varint flag bit collided with the continuation bit;
  rounding *deltas* accumulates error (quantise positions first, then delta); a string column sorted as
  text (`10` before `2`). Fix: a round-trip verifier **plus looking at screenshots** (Playwright).
- Generate the JS decoder from the Python encoder so they cannot drift.
- After refactoring shared scripts, re-run an existing group and confirm **byte-identical** output.
- Assert every text replacement when patching templates; never edit built output.

**Tooling**: many tools live in conda envs or private folders and are not on PATH - search
`~/miniconda3/envs/*/bin/`, search databases by file-type marker (`*.dbtype`, `*.hmm`), and grep old
project logs before concluding something is missing.

---

## 7. Open items and caveats

- Other four groups: only steps 01-10 done; islands/structures/dashboard are one command each
  (*C. difficile* is 5.3x larger).
- Resistance on the SNP trees (goal 2) untouched.
- KOfam / KEGG pathways not run (installed).
- `run_islands.py` never run end to end in one pass.
- Cluster `CFR-016271` is labelled `astA` but is an IS256 transposase (label from a 114 bp hit inside the
  CDS) - documented, not removed.
- The `qnrB` cluster merges plasmid and chromosomal alleles at 90% identity; needs allele-level assignment.

---

## 8. How to reuse this in a new project

1. **Scope first**: pick the group, check its size (>20k isolates -> subset genomes, use BigQuery), check
   species composition, note and pin all release/database versions.
2. **Cheap answers before expensive ones**: gene lists + catalog give profiles for every isolate in
   minutes; download genomes only for the questions that need location or sequence.
3. **Put every organism-specific choice in one config entry** (release id, input paths, point-mutation
   taxon, key subclasses, marker subclass, spot-check assertions that fail loudly).
4. **Choose what to report from the data** after the first annotation pass, never from habit; identify
   intrinsic genes before quoting any percentage.
5. **Build QC into every step** (counts that must match upstream, spot checks, controls with known
   answers) and write a README next to every output folder.
6. **Pangenome-level questions**: cluster proteins (MMseqs2 90/80), anchor on clusters, define genomic
   neighbourhoods with fixed gene windows, track contig-edge truncation, keep a backbone control tier.
7. **Protein-level questions**: variation per lineage, local structure prediction with a memory protocol
   and an external fallback, then two independent annotations (Pfam domains + Foldseek folds).
8. **Presentation**: offline self-contained pages with a verified data codec; one page per question,
   honest captions ("80 of N drawn, all counted").
9. **Keep a living status file** (done / missing / caveats) and update it whenever a step finishes.

---

## 9. Tools and versions used

Python 3 + pandas; NCBI `datasets`/`dataformat` v2; AMRFinderPlus catalog 2026-08-07.1; MMseqs2; MAFFT;
ESMFold (`facebook/esmfold_v1`, PyTorch, RTX 5070 Ti 16 GB); ColabFold (Google Colab) for 5 long proteins;
py3Dmol; HMMER 3.4 with Pfam 38.0 (from InterProScan 5.76); Foldseek 10 against PDB and AlphaFold
Swiss-Prot v6; Playwright for page screenshots. Machine: WSL2, 32 threads, 64 GB RAM.

## 10. The detailed documents this condenses

`README.md` (folder guide) · `docs/PROJECT_STATUS.md` (done / missing) · `docs/DATA_DOCUMENTATION.md`
(every file and column, QC, limitations) · `docs/WORKPLAN_new_species.md` (steps 01-10 for a new group,
traps, checklist) · `docs/WORKPLAN_stress_virulence.md` (steps 23-29) · `docs/WORKPLAN_dashboard.md`
(dashboard design, schemas, verification) · per-folder `README.md` files.
