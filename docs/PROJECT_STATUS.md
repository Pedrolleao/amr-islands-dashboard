# Project status - what is done and what is not

*Written 23 September 2026, updated 24 September (dashboard complete).* One place to see the whole project. The detail lives in the per-folder
READMEs and in `DATA_DOCUMENTATION.md`; this file says what exists, for which organism, and what is left.

## The three original goals

Stated 19 September 2026, for the NCBI Pathogen Detection Citrobacter group:

1. **Per-isolate resistance profiles** - done, for all five groups.
2. **Linking resistance to SNP clusters and trees** - *partly*. Profiles are joined to SNP clusters and
   every variation figure counts by lineage rather than by genome, so clonal expansion never inflates a
   number. What is **not** done is putting resistance onto the trees themselves (`SNP_trees/*.tar.gz`).
3. **Gene and protein-level study** - done and then some, for *C. freundii*: islands, per-position
   variation, predicted structures, domain architecture and fold matches.

## The pipeline

| steps | what they do | runner |
|---|---|---|
| 01-10 | metadata, resistance profiles, genome download, AMR sequences, antibiograms, group README | `run_all.py` |
| 11-22 | gene index, protein clusters, anchors, resistance islands, figures, 17 QC checks | `run_islands.py` |
| 23-26 | per-position variation, ESMFold structures, py3Dmol views, ColabFold ingest | per set |
| 27-29 | Pfam domains, Foldseek folds, folding a too-long protein by module | per set |
| 30-37 | dashboard payloads, the three pages, and a round-trip check of the packing | `run_dashboard.py` |

Steps 23-28 take `--set amr|stress|virulence` (27 also takes `all`), and write to
`Variation[_set]/`, `Structures[_set]/`, `Domains[_set]/`.

## Where each organism group stands

| group | isolates | genomes | steps 01-10 | islands 11-22 | clusters 23-28 |
|---|---|---|---|---|---|
| *Citrobacter freundii* | 7,287 | 6,063 | done | done | done |
| *Citrobacter portucalensis* | 581 | 432 | done | **not started** | **not started** |
| *Clostridioides difficile* | 39,708 | 31,900 | done | **not started** | **not started** |
| *Clostridium perfringens* | 5,706 | 5,464 | done | **not started** | **not started** |
| *Clostridium botulinum* | 1,676 | 1,449 | done | **not started** | **not started** |

## Citrobacter freundii in detail

**Islands (steps 11-22).** 208,757 protein clusters from 6,061 genomes; 409 carry resistance evidence and
anchor. 66,005 islands: 17,665 `tier1_multi_acquired`, 29,256 `tier2_single_acquired`, 19,084
`tier3_core_only`. 12,338 recurrent island types. 74% of islands run off a contig edge, so the clean set
for any claim about size or composition is the 1,419 tier-1 islands with a complete window. 17/17 QC
checks pass.

**The three cluster sets (steps 23-28).**

| set | clusters | variation figures | structures | viewer pages | Pfam named | fold matched |
|---|---|---|---|---|---|---|
| amr | 411 | 253 | 411 | 411 | 401 | 410 |
| stress | 328 | 237 | 328 | 328 | 309 | 264 |
| virulence | 67 | 40 | 67 (5 ColabFold, `clbB` in 3 modules) | 69 | 60 | 63 |

**`clbB` was folded in three modules on 24 September** (`scripts/29_fold_by_domain.py`), so every one of
the 67 virulence clusters now has a structure. Its cluster row reads `folded_by_domain` and no whole-chain
PDB exists, by design.

**The five ColabFold structures were taken in on 23 September** (`sigA`, `icmF`, `tsh`, `clbN`, `sslE`),
each verified against its representative sequence by residue count, confidence 70.6-87.4, all five matching
the expected fold. Only `clbB` has no structure.

**Pfam over the whole pangenome** (`--set all`) finished on 24 September, 11.6 hours: all 208,757
representatives, 168,110 domains retained, 126,850 clusters (60.8%) named. Output and its README are in
`Domains_all/`. The gradient is the story: 99.4% of core clusters are named and 59.9% of cloud ones, and
the 81,907 unnamed are short (median 116 aa), rare (43% in one genome) and mostly annotated only as
hypothetical.

## The dashboard (steps 30-37)

Three self-contained pages over the finished analysis, **built and verified on 24 September**:
`Citrobacter_freundii/Dashboard/{islands,clusters,antibiograms}.html`, opened straight from the folder.
Islands with gene maps, architecture types and aligned views per type, cluster or AMR class; one
record per analysed cluster with its variation, Pfam architecture, predicted fold and Foldseek hits;
and the 75 measured antibiograms beside each isolate's genotype. Tier-3 backbone windows are carried
but hidden behind a filter. Design, schemas, what was built and the checks that pass are in
[`WORKPLAN_dashboard.md`](WORKPLAN_dashboard.md); the folder's own `README.md` explains the pages,
the payloads and the codec. `python3 scripts/run_dashboard.py` rebuilds everything in about a minute.

| payload | size | contents |
|---|---|---|
| `islands.json` | 8.0 MB | 66,005 islands of all three tiers (tier 3 hidden by default on the page), their architecture types |
| `genemaps.json` | 6.0 MB | 695,638 genes, packed to about 7.5 characters each, with each cluster's AMR class |
| `cluster_index.json` | 4.6 MB | all 208,757 clusters, searchable |
| `structures.json` | 1.7 MB | 808 backbone traces, 210,971 residues |
| `clusters.json` | 0.9 MB | the 806 analysed clusters in full |
| `antibiograms.json` | 0.2 MB | 75 tested isolates: calls, measurements, genotype, islands and their gene maps |

`36_verify_payloads.py` decodes the packing again and compares it against the source tables; it
passes. Three encoder bugs have been caught this way, none of which crashed (the last: genes sorted
by position as text, `10` before `2`). Templates are in `scripts/dashboard/`, one copy for every group.

## What is missing, in the order I would do it

1. **The other four organism groups.** Islands and cluster analyses have only ever been run for
   *C. freundii*. Each is one command. *C. difficile* is 5.3x the size of *C. freundii* and is the only one
   likely to strain the machine.
2. **Goal 2: resistance on the SNP trees.** Untouched. The trees are in `<group>/SNP_trees/*.tar.gz`.
3. **The interactive island browser** is superseded by the dashboard, which is built.
4. **KOfam.** Installed on this machine and would add KEGG pathway membership. A different question from
   domain architecture, not a substitute for it. See `WORKPLAN_stress_virulence.md` section 9.
5. **Re-run `run_islands.py` end to end.** Steps 11-22 were run one at a time, so there is no
   `qc_islands_report.txt` from a single clean pass.

## Two open caveats that are not bugs

- **`CFR-016271` is labelled `astA` but is an IS256 transposase.** Its 402-aa proteins cannot be the 39-aa
  EAST1 enterotoxin. The label came from a 114 bp hit falling inside the transposase CDS, attached at tier
  `T4_overlap_partial` with `propagates_label = 0`; step 15 reports weak links by design and
  `cluster_sets.py` selects on `element_type`. Structure (pLDDT 93 over 402 aa) and Pfam
  (`Transposase_mut`, 94% coverage) both agree it is a transposase. It is documented, not removed:
  whether to drop it from the virulence set is a judgement call that has not been made.
- **The `qnrB` cluster is two populations.** At 90% identity the plasmid-borne and chromosomal alleles fall
  in one cluster and no threshold separates them; it needs allele-level assignment.

## Reproducibility notes

- The AMRFinderPlus catalog is pinned at `2026-08-07.1`, Pfam at 38.0 via the InterProScan path, AlphaFold
  Swiss-Prot at v6 (2025_03).
- `Structures/structure_summary.tsv` (the amr set) predates the `status` column; read it with a default.
- The 61 virulence structures have an empty `ptm` and `device = recovered_from_pdb`, because their rows
  were rebuilt from the PDB files after a crash. Re-run step 24 with `--redo` if pTM is ever needed.
