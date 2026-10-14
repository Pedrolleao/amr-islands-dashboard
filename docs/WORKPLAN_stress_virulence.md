# Workplan – repeating the cluster analysis for the stress and virulence clusters

> The dashboard over all of this is planned in
> [`WORKPLAN_dashboard.md`](WORKPLAN_dashboard.md).

> **Status of the whole project, including what is still missing, is in
> [`PROJECT_STATUS.md`](PROJECT_STATUS.md).** Sections 0-7 below are the original plan for the stress and
> virulence sets, which is now **complete** apart from the six proteins listed in section 8. Sections 8
> and 9 were added on 23 September and are the outstanding work.

*Written 22 September 2026, from the resistance run on the Citrobacter freundii group (PDG000000039.1314).*

What we did for the **411 resistance-associated protein clusters** – per-position variation (step 23),
ESMFold structures (step 24), py3Dmol views (step 25) – repeated for the **stress** clusters and then for
the **virulence** clusters. Nothing in the island pipeline changes: stress and virulence stay non-anchors,
as decided when the islands were designed, and keep their own colours in the island figures.

Order: **stress first, virulence second.** Stress is the larger and simpler set, it validates the code
change, and it is the group that actually travels with resistance on the same transposons.

---

## 0. What the three groups look like

From `Pangenome/cluster_annotation.tsv` (`element_type`), for the *C. freundii* group:

| group | clusters | proteins | clusters with >=2 distinct sequences | residues to fold |
|---|---|---|---|---|
| AMR (done) | 411 | 100,352 | 253 | 100,337 |
| **STRESS** | 328 | 156,576 | 237 | 69,883 |
| **VIRULENCE** | 67 | 11,792 | 40 | 39,823 |

Stress is metal and biocide resistance: arsenic (`arsC`, `arsD`, `arsA`), mercury (`merA`, `merP`, `merE`,
`merD`, `merR`), acid survival (`asr`), quaternary ammonium efflux (`qacE`).
Virulence is invasins and siderophores: `sinH`, yersiniabactin (`ybtP`, `ybtQ`), aerobactin (`iucABCD`),
salmochelin (`iroC`), enterotoxin `senB`, plus the colibactin cluster.

The pinned AMRFinderPlus catalog covers both groups (259 STRESS and 1,034 VIRULENCE entries in
`ReferenceGeneCatalog.txt`, with their proteins in `AMRProt.fa`), so the "closest catalog protein" naming in
step 23 works unchanged.

---

## 1. The one code change (do this once, before the stress run)

Steps 23–25 currently read the resistance cluster list from `Islands/anchor_clusters.tsv`. They need to take
a **cluster set** instead. Everything else in them is group-agnostic.

**1a. New helper `scripts/cluster_sets.py`** – returns the cluster list for a set name, in the column shape
the later scripts already expect:

| set | source | rows |
|---|---|---|
| `amr` | `Islands/anchor_clusters.tsv` (unchanged, so the resistance output stays reproducible) | 411 |
| `stress` | `cluster_annotation.tsv` where `element_type == "STRESS"` | 328 |
| `virulence` | `cluster_annotation.tsv` where `element_type == "VIRULENCE"` | 67 |

For the stress and virulence sets, map the columns onto the names steps 23–25 use:
`stress_symbols` / `virulence_symbols` -> `amr_gene_families` and `amr_symbols`; `element_type` ->
`amr_classes`; `anchor_tier` -> the string `NOT_AN_ANCHOR`; `point_symbols` -> empty (neither group has
catalog point mutations, so the wild-type-reference rule of step 23 is a no-op and must not crash on empty).
Keep `cluster_id`, `rep_md5`, `n_seqs`, `n_members`, `prevalence`, `product_top` as they are.

**1b. Add `--set {amr,stress,virulence}` to steps 23, 24 and 25**, default `amr`. It selects the cluster
list and the output folder:

| step | output for `amr` (existing, untouched) | output for another set |
|---|---|---|
| 23 | `Variation/` | `Variation_stress/`, `Variation_virulence/` |
| 24 | `Structures/` | `Structures_stress/`, `Structures_virulence/` |
| 25 | `Structures/views/` | `Structures_stress/views/`, … |

Separate folders, not a shared one with a column: each folder keeps its own README, and the resistance
output stays exactly as it is now.

**1c. Re-run the resistance set once after the change and compare checksums** of
`Variation/cluster_variation_summary.tsv` and `Structures/structure_summary.tsv`. They must be identical –
the same check we ran when the pipeline was first parametrised by organism group. Do not start the stress
run until this passes.

---

## 2. Stress – the run

```bash
cd latest_snps
python3 scripts/23_cluster_variation.py --group Citrobacter_freundii --set stress
python3 scripts/24_structure_prediction.py --group Citrobacter_freundii --set stress
tools/viz_env/bin/python scripts/25_structure_views.py --group Citrobacter_freundii --set stress
```

| step | what it does | expected scale | time |
|---|---|---|---|
| 23 | MAFFT alignment + per-position variation for the 237 stress clusters with >=2 sequences | ~237 figures, ~65,000 position rows | ~4 min |
| 24 | ESMFold structure of all 328 representatives | 328 PDB pairs, ~70,000 residues | ~25 min on the GPU |
| 25 | py3Dmol page + stills | 328 pages, ~565 stills | ~6 min |

Notes specific to stress:

- **Two long proteins**: the SilA copper/silver efflux permease (1,047 aa) and the ClpK heat-shock ATPase
  (947 aa). Both are inside the range that worked for OqxB (1,051 aa), so the existing fallback
  (GPU chunk 16 -> chunk 4 -> CPU) handles them. No change needed.
- **Variation will be high**: the mer and ars operons are on plasmids and transposons, so expect more
  polymorphism than in the chromosomal resistance clusters. The largest stress cluster has 296 distinct
  sequences, more than any resistance cluster except blaCMY.
- **No catalog point mutations**, so no figure carries a catalog-site marker, and every reference is simply
  the most common full-length sequence.

## 3. Virulence – the run

Same three commands with `--set virulence`.

| step | expected scale | time |
|---|---|---|
| 23 | ~40 figures (40 clusters have >=2 sequences) | ~2 min |
| 24 | 67 structures, ~40,000 residues | see the length problem below |
| 25 | 67 pages, ~107 stills | ~2 min |

**The length problem, and the decision it needs.** Virulence proteins are much longer than resistance ones:
median 456 aa, and **10 clusters above 900 aa**, including

| cluster | length | protein |
|---|---|---|
| CFR-014064 | 3,205 aa | colibactin PKS/NRPS ClbB |
| CFR-080248 | 1,504 aa | lipoprotein metalloprotease SslE |
| CFR-014197 | 1,455 aa | colibactin NRPS ClbN |
| CFR-115901 | 1,377 aa | hemoglobin-binding autotransporter Hbp |

ESMFold memory grows with the square of the length. OqxB at 1,051 aa already peaked near 14 GB of the
16 GB card at chunk 16. Plan for three bands:

- **up to ~1,100 aa** (57 of the 67): the current settings work.
- **1,100–1,500 aa** (3 clusters): force chunk 4, and expect a CPU fallback taking tens of minutes each.
- **above 1,500 aa** (1 cluster, ClbB at 3,205 aa): do **not** fold whole on this machine. Options, in the
  order I would try them: (a) skip it and record why in the summary; (b) fold the domains separately after
  splitting on the NRPS/PKS module boundaries; (c) leave it for a machine with more GPU memory. A 3,205-aa
  single-chain ESMFold run needs roughly 5 GB per pair tensor and would not finish in reasonable time on the
  CPU here.

Whichever is chosen, the cluster must appear in `structure_summary.tsv` with an explicit reason, never be
silently missing – the same rule that caught the missing dfrG row in the resistance run.

**Decided, 23 September 2026.** Option (a), with a rule that covers both kinds of failure. Step 24 defers a
protein when it is longer than `--max-len` (1,500 aa) **or** when it passes the `--max-seconds` cap, one hour
per protein, and defers everything still queued with it since the queue is sorted by length. `clbB` and
`sslE` were deferred on length before any folding started. Every deferred protein is written to
`Structures*/colabfold/<name>.fasta`, one file per protein, so each can be run through ColabFold separately.
The returned PDB goes in `colabfold/returned/` and `scripts/26_ingest_colabfold.py` files it: it checks the
residue count against the representative sequence, writes `pdb/` and `pdb_variation/`, and rewrites the
summary row as `predicted` with `source=colabfold`. Nothing is ever silently missing, which was the rule.

---

## 4. Checks to repeat, per set

Copied from the resistance run; all of these are printed by the scripts themselves.

**Step 23**
- protein copies per cluster equal `n_members` in `cluster_annotation.tsv` (was 253/253)
- every reference is unchanged by the alignment
- flat figures are explained in their subtitle (fragments / insertions / length outside the reference)

**Step 24**
- structures written = clusters in the set, failures listed by name
- residue count of each PDB equals the representative's length
- median mean pLDDT, count below 70 and below 50; fragment vs non-fragment medians
- **expect lower confidence than the resistance set.** Beta-lactamases are among the best-represented
  families in any protein database; metal-efflux membrane proteins and giant NRPS modules are not. A median
  pLDDT in the 70s for virulence would be normal, not a failure.

**Step 25**
- pages written = clusters in the set; stills = clusters + clusters with variation
- no still under 5 kB (a blank render)

**A biological sanity check per set**, equivalent to gyrA T83I and KPC H272Y in the resistance run:
- stress: the `merA` cluster should be tightly conserved in its NADPH-binding motif, while `qacE` /
  `qacEdelta1` should split into two length classes – the truncated `qacEdelta1` is the one carried by class
  1 integrons, and it should appear as a separate cluster or as a fragment flag.
- virulence: the aerobactin genes `iucA`–`iucD` should be nearly identical across lineages (one plasmid
  lineage), while `sinH` should be the most variable – it is a surface-exposed invasin under immune
  selection.

---

## 5. Traps already paid for – do not rediscover them

| trap | what happens | what to do |
|---|---|---|
| ESMFold loaded with `.cuda()` before `.half()` | fp32 model fills the 16 GB card, anything over ~600 aa spills into WSL system memory and looks hung (gyrA: >10 min vs 92 s) | load, then `model.esm = model.esm.half()`, then `.cuda()` – 8.4 GB |
| proteins over ~900 aa | out of memory at the default chunk size | chunk 16 or 4; the script already falls back automatically |
| `pkill -f <script>` from a tool call | also matches the calling shell and kills the command | kill by PID |
| `--redo` on one cluster | overwrote the whole summary table | fixed: the merge now keeps rows the run did not touch. Verify the row count after any partial re-run |
| py3Dmol `_make_html()` | pulls 3Dmol from a CDN and hides the viewer in a global `viewer_<digits>` created inside a promise | rewrite the URL to `views/lib/3Dmol-min.js`, scrape the uuid, poll for the global |
| Playwright Chromium | the pip package expects a newer build than the cached one | launch with `executable_path=` the cached Chromium + swiftshader flags |
| system python is externally managed (PEP 668) | `pip install` refuses | use `tools/viz_env/bin/python` for step 25 |
| a common sequence may not be the sensible reference | for gyrA the commonest sequence was the resistant I83 mutant | only matters where catalog point mutations exist, so not for stress or virulence |

---

## 6. Documentation to write at the end of each set

Mirror `Structures/README.md` and `Variation/README.md`:

- what was folded and from which representative;
- the confidence table, and which clusters fell below 70;
- the group-specific caveats. For stress and virulence, two matter most: **these are genus-wide numbers**
  (the group is not one species, so background variation partly reflects species differences), and
  **virulence screening by NCBI is far less thorough for *Citrobacter* than for *Salmonella* or *E. coli***,
  so 67 clusters understates what is present.
- update `docs/DATA_DOCUMENTATION.md`, which still stops at step 10 and covers neither the islands
  (11–22) nor the cluster analyses (23–25).

---

## 7. Total cost

| | stress | virulence |
|---|---|---|
| step 23 | ~4 min | ~2 min |
| step 24 | ~25 min GPU | ~20 min GPU + the long-protein decision |
| step 25 | ~6 min | ~2 min |
| disk | ~120 MB structures + ~380 MB views | ~40 MB + ~90 MB |

Plus the one-off code change in section 1, and the reproducibility re-run of the resistance set.

---

## 8. ClbB, folded by domain - DONE 24 September 2026

*Agreed 23 September 2026.* `CFR-014064__clbB`, 3,206 aa, is the one structure still missing. ColabFold on
Google Colab crashes on it repeatedly, which is expected: AlphaFold2 memory grows with the square of the
length and this exceeds what a Colab runtime holds. Our own measurement puts ESMFold's activations alone at
about 74 GB for this sequence.

**The decision is option (b) of section 3: fold it in modules, not as one chain.** ClbB is a hybrid
NRPS/PKS built from discrete catalytic domains joined by flexible linkers. A single-chain model is both the
hardest thing to compute and the least trustworthy to interpret, because the arrangement of the modules is
genuinely dynamic; the confident part of any successful whole-chain prediction would be the individual
domains anyway.

**Done on 24 September** with `scripts/29_fold_by_domain.py`, in 506 s. The 14 domains grouped into three
modules of 1,103 / 1,127 / 1,016 aa, cutting at the NRPS/PKS junction and at the condensing-to-reductive
transition without being told to. Confidences 89.2, 87.4 and 80.0, weighted mean 85.7, which is higher
than any of the five whole-chain ColabFold models. All three match a known fold at TM 0.94-0.95. All
3,206 residues covered. Details in `Citrobacter_freundii/Structures_virulence/README.md`.

**The steps that were followed, kept for the next long protein:**

1. Annotate the sequence with the NRPS/PKS domain models (condensation, adenylation, carrier protein,
   ketosynthase, acyltransferase, ketoreductase) using HMMER, which is installed. This falls out of the
   domain-annotation step of section 9 and should not be done separately.
2. Cut on the domain boundaries HMMER reports, keeping the linkers with the preceding domain and allowing a
   small overlap so no residue is lost.
3. Fold each module with step 24. Each lands in the 400-600 aa range, which the language-model eviction
   protocol handles in well under a minute at high confidence.
4. Record ClbB in `structure_summary.tsv` as folded **by domain** rather than as a single chain, with the
   boundaries in the row, and keep the per-domain PDBs under a name that marks them as parts. It must never
   read as a whole-chain model, and it must never read as missing.

If a whole-chain model is ever wanted anyway, the levers that make it survivable on Colab are, in order of
effect: `--max-msa 32:64`, then `--num-models 1 --num-recycle 1`, no Amber relaxation,
`TF_FORCE_UNIFIED_MEMORY=1` with `XLA_PYTHON_CLIENT_MEM_FRACTION=4.0`, on an A100 runtime.

---

## 9. Domain and fold annotation of every cluster representative (steps 27-28)

*Built 23 September 2026.* Two independent answers to "what is this protein", for the representative of
each protein cluster. Neither replaces the other and both are cheap, so both are run.

**Everything needed was already on this machine** and nothing was downloaded: Pfam 38.0 comes pressed
from the InterProScan 5.76 installation (`/home/dmin/Program/interproscan-5.76-107.0/data/pfam/38.0/`),
which pins the version by path, and Foldseek 10 with its PDB database sits in
`~/.local/share/archaeadefense/foldseek/`. HMMER 3.4 is on the PATH.

**Step 27, `scripts/27_pfam_domains.py --set {amr,stress,virulence,all}`** runs `hmmsearch` (not
`hmmscan` - same answers, an order of magnitude faster) with each family's own gathering threshold, then
applies Pfam clan competition so homologous families do not count one region twice. Output per set in
`Domains[_<set>]/`: `pfam_domains.tsv` (one row per domain, with coordinates) and `cluster_domains.tsv`
(the ordered architecture per cluster).

**Step 28, `scripts/28_foldseek_hits.py --set {amr,stress,virulence}`** searches the step-24 structures
against the PDB. It carries each model's own pLDDT into the output, because a fold match against a bad
model means nothing. `same fold` is TM-score >= 0.5.

**What the first run showed, and why both methods are needed.**

| set | clusters | Pfam names | fold matches | named by neither |
|---|---|---|---|---|
| amr | 411 | 401 | 408 | 0 |
| stress | 328 | 309 | 228 | 18 |
| virulence | 67 | 60 | 52 | 3 |

In the resistance set the two are almost interchangeable, and every one of the 10 clusters Pfam missed
was rescued by structure. In the stress set they are not: 18 clusters have neither a domain nor a fold
match, and they are the same small periplasmic metal-binding proteins that folded at pLDDT 41-54 -
`silE`, `pcoE`, `silF`, `asr`. Three independent methods agreeing that there is nothing to see is itself
the finding.

**Measured cost** (this machine, 32 threads): about 138 s of fixed database streaming plus 0.84 ms per
residue. So a set of a few hundred clusters costs 3-4 minutes, and the whole pangenome, 208,757
representatives and 47.5 M residues, costs about 11 hours. That is an overnight job, not the two to three
days first guessed, so there is no need to restrict it to non-singleton clusters.

**AlphaFold added, 23 September 2026.** The AlphaFold database WAS already here, in another project's
reference folder rather than under any name a search for "alphafold" would find:
`.../ArchaeaHQ-DefenseFinder/Discovery/Explore/validation/refs/afdb_sp`, Swiss-Prot v6 (2025_03),
590,183 models with taxonomy. Step 28 now searches PDB and AlphaFold together (`--db`), and the best hit
is taken across both. There is no full AFDB50 or proteome database on the machine, confirmed by scanning
every file over 500 MB on both drives.

| set | fold matches, PDB only | with AlphaFold added | named by nothing |
|---|---|---|---|
| amr | 408 | 410 | 0 |
| stress | 228 | 264 | 13 |
| virulence | 52 | 55 | 3 |

It paid off most where it was needed most. Of the 19 stress clusters Pfam could not name, PDB rescued 1
and AlphaFold rescued 6, naming `merD` (TM 0.94), `pcoD` (0.97), `arsC` (0.88) and `terW` (0.81).

**One lesson to carry forward: TM-score is the wrong test for a disordered protein.** The `silE` and
`pcoE` clusters come back against the reviewed *silver-binding protein* at TM 0.28-0.45, below the
same-fold line. The identification is almost certainly right; the score is low because neither the query
nor the target has a stable shape to superpose. Always read the target's name before dismissing a weak
score.

**Still to add:** KOfam, which is installed (`exec_annotation` 1.3.0 in the `ewaste-pipeline` env, 28,160
KOs with curated thresholds) and would add KEGG pathway membership - a different question from domain
architecture, not a substitute.
