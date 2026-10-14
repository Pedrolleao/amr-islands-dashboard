# Workplan – the Citrobacter freundii resistance dashboard

*Written 24 September 2026. The data layer, the front end and every addition below were built and verified that same day. Sections 1-4 are the design and the schemas, section 5 what was built and in what order, section 6 the checks that pass.*

Three self-contained HTML pages over the finished analysis, so that answering a question about one
gene family stops meaning opening five folders.

## 1. What was decided, and why

| question | decision |
|---|---|
| cluster scope | The 806 analysed clusters get full detail. The other ~208,000 stay searchable with product, prevalence and Pfam architecture, so a gene in an island always has a name. |
| delivery | Three pages sharing one build script: `islands`, `clusters` and `antibiograms`. A separate `structures` page was built first and then folded into the cluster record the same day: it had no question of its own, and half of it duplicated the record's sections. `antibiograms` was added as a page of its own at the user's request: **they prefer a new page to inflating an existing one.** |
| islands | Tier 1 and tier 2 by default, 46,921 islands. The 19,084 tier-3 backbone-only windows were excluded at first (they exist to calibrate the tiering); the user asked for them, so they are carried too, hidden behind the tier filter (`?tier=2`) and labelled *chromosomal window*, never *island*. Every count, the distribution tab and the aligned views follow the filter. |
| structures | Inline packed CA traces drawn on a canvas, as a full-width panel of the cluster record. No dependency on the 1.2 GB of py3Dmol pages from step 25. |

Pages are separate files with no server, so **cross-page links carry state in the query string**
(`clusters.html?cluster=CFR-000010`), and anything two pages both need is duplicated in each.

## 2. The two references

**Layout and data organisation** come from the ArchaeaHQ DefenseFinder dashboard. Read the templates,
`…/ArchaeaHQ-DefenseFinder/Discovery/Dashboard/src/{islands,atlas,structures}_template.html`
(1,901 / 834 / 876 lines). **Never open the built HTML** next to them: each is one 8–11 MB line.
Its build script, `…/Discovery/scripts/23_build_pages.py`, is 38 lines and is the model for ours.

It uses **no JavaScript libraries at all**. Gene maps are hand-generated SVG polygons, every chart is
a list of `div`s with percentage widths, and the 3D viewer is a bespoke canvas CA-trace renderer with
Catmull-Rom splining and painter's-algorithm depth sorting (lines ~1544–1660 of the islands
template). The per-residue ribbon is generated SVG (`foldRibbon`, ~line 1697).

**Visual design** comes from the Simple brand book,
`~/Protocols/protocols_website/design/Simple/design-system_Simple.html` (`:root` at lines 10–25) and
`motion-components-simple.html` (lines 10–48 hold the superset token block to copy verbatim).

```css
--bg:#0E0E0E  --bg-2:#161616  --highlight:#1F1F1F
--fg:#E5E5E5  --fg-mute:#888888  --fg-dim:#555555
--border:#2A2A2A  --border-strong:#3A3A3A
--c-teal:#6DC7C8  --c-green:#61CA6E  --c-coral:#CE6A68  --c-purple:#B579CF  --c-mustard:#CDC37C
```

JetBrains Mono for everything, 12px/1.55, letter-spacing .01em. **No border radius, no shadows, no
gradients** — the only curve in the system is a 6px status dot. The signature component is `.frame`,
a 1px rectangle whose `[BRACKETED_LABEL]` sits on the top border and cuts it. Labels, buttons and
badges are uppercase with 0.08–0.22em tracking; prose stays lowercase. **Colour is categorical
metadata only** — thin outlines, small dots, or text, never fills. There is a single dark theme and
no light mode; the `Light/` folder is empty.

The source defines **no tabs, no nav and no form inputs**. Design them by analogy from its segmented
toggle (`motion-components-simple.html` lines 82–95): one hairline strip, 1px internal dividers,
active cell filled with `--highlight`. Inputs: transparent background, 1px `--border-strong`,
uppercase mono label, focus shown by a border-colour change to an accent, never a glow.

## 3. The data layer (scripts 30–34, 36, 37)

All take `--group`, so the other four organism groups work later. Output to
`Citrobacter_freundii/Dashboard/`. They do **not** use the `--set` machinery, because they merge all
three cluster sets into one payload.

| script | produces | measured |
|---|---|---|
| `30_dashboard_common.py` | the codec: base64-varint, zigzag, quality strings, `Dicts`, id helpers, and `js_decoder()` | `python3 scripts/30_dashboard_common.py` self-tests it |
| `31_cluster_payload.py` | `clusters.json`, `cluster_index.json` | 0.92 MB, 4.57 MB |
| `32_island_payload.py` | `islands.json` | 7.95 MB |
| `33_island_genemap_payload.py` | `genemaps.json` | 6.04 MB |
| `34_structure_payload.py` | `structures.json` | 1.66 MB |
| `37_antibiogram_payload.py` | `antibiograms.json` | 0.16 MB |
| `36_verify_payloads.py` | decodes them again and compares against the source TSVs and PDBs | exits non-zero on any mismatch |

Counts that should reproduce exactly: 66,005 islands (17,665 tier 1, 29,256 tier 2, 19,084 tier 3),
695,638 genes, 806 analysed clusters of which 530 have variation, 808 chains over 210,971 residues,
439 of them carrying variation; 75 tested isolates with 1,156 antibiogram calls. Before tier 3 was
added the island payloads held 46,921 islands and 350,064 genes; `islands.html` grew from 9.3 to
14.1 MB with it. Built pages: `islands.html` 14.1 MB, `clusters.html` 7.2 MB, `antibiograms.html`
0.2 MB.

### Three bugs already paid for — do not reintroduce them

1. **The gene flags must be a varint, not a raw character.** The pseudogene bit is 32, which is also
   the varint continuation bit, so a raw character makes the whole stream undecodable. It costs one
   extra character on the 10.7% of genes that are pseudogenes.
2. **Quantise positions, then take deltas — never round the deltas.** Rounding deltas lets the error
   accumulate along an island; it reached 46 bp by the end of one. Quantising first bounds every
   coordinate to within 5 bp, always.

3. **Sort `pos_in_island` as an integer.** `read_tsv` returns every column as a string, so the
   genes of an island came out `0, 1, 10, 11, …, 2, 3`: every coordinate was right, but the map
   clipped at the wrong end and every delta after the first was a long jump. Both step 33 and the
   verifier sorted the same way, so the round-trip check passed; step 36 now also checks that the
   decoded genes are in contig order.

All three were caught only by looking (36, or the built page) and none would have surfaced as a
crash, only as a drawing subtly wrong. Run step 36 after any change to the encoders.

## 4. Payload schemas — what the front end must decode

**The codec** (`30_dashboard_common.py`): 5 bits per character over
`ABC…XYZabc…xyz0…9-_`, bit 32 continues. `zig(n)` is zigzag; quality strings are one character per
residue, `chr(35 + round(value * 0.9))`, so a 0–100 value decodes as `(code - 35) / 0.9`.
`js_decoder()` emits the JavaScript half — **generate it, never hand-write it**, so the two cannot
drift. Cluster and island ids are dense integer ranges, so payloads store the bare integer:
`CFR-{n:06d}` and `CFRI{n:06d}`.

**`genemaps.json`** — `{bp:10, roles[8], flag_label:16, flag_pseudo:32, cluster[], clabel[],
cclass[], csub[], dicts{label[], amrclass[], amrsubclass[]}, g{islandNum: packed}}`. `cclass[i]` and
`csub[i]` are the AMR class and subclass sets of `cluster[i]` (indices into the dicts, −1 for none),
added for the class view. Decode `g[n]` as a flat varint stream: first value is the
island's quantised start; then per gene four values, `zigzag(gapFromPreviousEnd)`, `width`,
`cs`, `flags`, plus `labelIndex+1` only when `flags & 16`. Multiply positions by `bp`.
From `cs`: `cluster = (cs >> 1) - 1` (−1 means none), `strand = cs & 1`. From `flags`: role in the
low 4 bits, bit 16 label override, bit 32 pseudogene. A gene's label is `clabel[cluster]` unless it
carries its own. Roles, in the order the array gives them: `amr_anchor_core`, `amr_anchor_mobile`,
`amr_anchor_point`, `hypothetical`, `mge`, `other`, `stress`, `virulence`.

**`islands.json`** — `{fields[30], type_fields[13], dicts{}, ids[], isl[], types{}, cooc[]}`.
`ids[i]` is the island number of `isl[i]`; the kept set is not dense, so the parallel array is
needed. Field order is in `fields`. `types` is keyed by type number, `cooc` rows are
`[clusterA, clusterB, geneAIdx, geneBIdx, nPds]`.

**`clusters.json`** — `{fields[31], dicts{}, pfam[[acc,name,desc]], recs{}, var{}, dom{}, fold{},
isl{}}`, all keyed by cluster number. `var[n]` is `{ref, lin, cop, sites[], top[[pos,fracx1000]]}`
where `lin` and `cop` are quality strings over reference positions. `dom[n]` is
`[[pfamIdx, start, end]]`. `fold[n]` is `[[dbIdx, targetIdx, descIdx, tmx1000, fidentx1000,
sameFold]]`. `isl[n]` is up to 60 island numbers, tier 1/2 first; the totals are the `n_islands` (tier 1 and 2)
and `n_islands_t3` fields.

**`cluster_index.json`** — a structure of arrays over all 208,757 clusters. **Position is the id**:
slot `i` describes `CFR-{i+1:06d}`. Arrays `nm` (members), `pv` (prevalence per mille), `pi`
(product index), `ai` (architecture index), `ci` (class index), `st` (set: 0 amr, 1 stress,
2 virulence, −1 not analysed). Build the search client-side over the ~21,000 dictionary strings.

**`antibiograms.json`** — `{absent, group_isolates, dicts{antibiotic, units, method, standard,
symbol, family, class, subclass, species, pds}, targets{abxIdx: [subclass tokens]}, isolates[]}`. An
isolate is `{t, bs, asm, sp, pds, genome, islands[], islands3[], tests[[abx, pheno 0 S/1 I/2 R, sign, value,
units, method, standard]], genes[[symbol, family, class, subclass, point, cluster or −1, partial]]}`.
`targets` is the hand-written antibiotic-to-subclass map in step 37; it only decides which genotype
rows the page flags as *expected*. `gm` is a slice of `genemaps.json` for the tested isolates'
islands and tier-3 windows, re-indexed to a local cluster list and decoded back as a check,
so the page can draw them. Written empty, with `absent: 1`, when step 09 never ran.

**`structures.json`** — `{dicts{}, chains[], recs{}}`. A chain is `{c, m, n, name, seq, xyz, pl,
var?, set, src, frag}`; `xyz` is the varint-zigzag of deltas at 0.1 Å, x/y/z interleaved; `pl` and
`var` are quality strings; `var` is **absent** for the 369 chains that are all zero. `recs[n]` is
`{chains:[idx]}` plus, for clbB only, `{whole:3206, by_domain:1, modules:[[start,end,chainIdx,
domainsIdx]]}`.

## 5. The front end, as built (24 September 2026)

The templates live in **`scripts/dashboard/`**, not in `Dashboard/src/` as first planned: they are
code, one copy serves every organism group, and `35_build_pages.py` writes the built pages into
`<group>/Dashboard/`. The build inlines `core.css`, then the codec from `js_decoder()`, a `META`
object (group, prefixes, build date) and `core.js`, then each `<!--DATA:__VAR:file.json-->` marker.
Every length in the stylesheet and templates is 1.25x the brand book's (base 15px): the user judged
the 12px original right only at 125% browser zoom, so that is now 100%. Scrollbars are styled in the
page's own greys.

| file | what it holds |
|---|---|
| `core.css` | the Simple token block, `.frame`, `.btn`, `.badge`, `.data-table`, `.kv`; tabs, `.seg`, inputs, `.rail`, `.bars`, `.tip` and scrollbars designed by analogy; the gene-role and confidence-bucket variables |
| `core.js` | `esc/fmt/pct`, dictionary lookup, `QS` query-string state and the cross-page `link` helpers, the role table and `roleChip()` (chip wording identical to the legend: *acquired resistance*, never *mobile*, for `amr_anchor_mobile`), bucket functions, `wireTabs`, `wireSeg`, `barsHTML`, tooltips |
| `islands_template.html` | tabs `[ISLANDS] [DISTRIBUTION] [HOW TO READ]`. The islands tab is the `330px \| 1fr \| 400px` grid, collapsing at 1600 and 1100 px: rail with search, species, lineage, AMR class and subclass selects, the tier segment `[TIER 1+2] [1] [2] [3] [ALL]`, integron and MGE toggles; centre the fit-or-detail SVG gene map, then the gene table; right the record and the architecture-type panel. **Aligned views** replace the centre: every island of a *type* (`?type=N`), of a *cluster* (`?cluster=X`, aligned on it), or of an *AMR class or subclass* (`?amrclass=`, `?subclass=`, aligned on whichever gene of the class each island carries), flipped so the anchor points right, capped at 80 rows in id order, with a **conservation table** over all of them: islands carrying each cluster, median offset from the anchor, anchor share in the class view; hover highlights, click filters. All of it follows the tier filter, and a cluster filter says how many isolates, lineages, species and types the list spans. Distribution tab: eight bar-list charts over the tier filter, clickable into the views. Buttons say *80 of N drawn, all counted* when a cap applies. |
| `clusters_template.html` | segmented set filter with counts; pangenome-wide search (dictionary scan, then the parallel arrays); AMR class, subclass, pangenome class and tier selects, variation and fold toggles; sortable sticky-header table with *amr class* and *subclass* columns, capped at 400 rows; the record: header with badges and *islands · N* / *windows · N* buttons, a `columns:4 375px` masonry (identity, anchor tier, per-position variation, Pfam architecture, islands) whose sections are only added when their data exists, and last the full-width **structure panel**: the Catmull-Rom tube renderer with depth shading, confidence-or-variation toggle, module segments for clbB, the ribbon of confidence, variation and domains along the whole chain, and beside it the source, confidence and module facts with the ranked Foldseek hits beneath. The canvas is rebuilt per record by `mountStage()`. `?cluster=`, `?set=`, `?amr=`, `?subclass=` |
| `antibiograms_template.html` | tabs `[GENOTYPE × PHENOTYPE] [ISOLATES] [COVERAGE] [HOW TO READ]`. Pick an antibiotic: every gene family carried by at least the chosen minimum of tested isolates, with its resistance rate among carriers against non-carriers, lineage counts, the difference in points and an *expected* flag from the step-37 subclass map; pick a family for the same across antibiotics; both, or a row, for the isolate lists. Isolates tab, `330px \| 1fr \| 400px`: rail, the isolate's calls with measured values and the genes acting on each drug, its genotype by class (green dot = drawn below, hover picks the arrow out, links into the clusters page), and its islands and tier-3 windows drawn one row each on one scale. `?antibiotic=`, `?family=`, `?isolate=` |
| `35_build_pages.py`, `run_dashboard.py`, `Dashboard/README.md` | build (three pages, removes a stale `structures.html`), runner (asserts upstream tables and names the step to rerun; `--pages-only` for template edits), folder documentation |

**Order of the additions, all on 24 September, all at the user's request:** the 1.25x scale; the
type view with its conservation table; the same view for a cluster filter; AMR class and subclass
columns and filters on the clusters page (step 31 gained `amr_subclasses`); the structures page
merged into the cluster record and dropped, then the structure panel moved below the cards with the
fold matches beside it; the AMR class view on the islands page (step 33 gained `cclass`/`csub`); the
scrollbars; the antibiograms page (step 37); its drawn islands; tier 3 everywhere, hidden by default
(steps 31-33 and 37 carry it, step 31 counts it separately as `n_islands_t3`).

Decisions made while building, beyond section 1: green marks catalog sites, same-fold matches and
*drawn in an island* and nothing else; the pLDDT buckets are captioned everywhere they are drawn; a
cluster outside the three sets opens on the clusters page with an identity section only and a badge
saying so; the spin button is disabled under `prefers-reduced-motion`; a type is the same
arrangement of the same informative genes (step 19's `arch_core_key`), so the cluster and class
views are the wider questions and the pages say so.

**Light theme, added 29 September 2026 at the user's request.** The brand book has none, so it was
designed here: a sun/moon button in the nav strip of every page (icon shows the theme it switches to)
sets `data-theme="light"` on `<html>`, and `:root[data-theme="light"]` in `core.css` redefines every
token under the same name - greys inverted, each accent darkened to 4.7:1 or better as text on `--bg`
and `--bg-2`, role tints lowered to 20-24%. Dark stays the default. The choice is kept in
`localStorage` (`amr-dashboard-theme`, guarded, since storage can be blocked on `file://`) and a
one-line script that step 35 puts in `<head>` applies it before first paint. The structure canvas
reads resolved colours, so its outline (`--edge-rgb`, `--edge-k`) and the module veil (`--veil`) are
tokens too, and it redraws on the `themechange` event. Checked with
`~/tools/shotter/theme-check.mjs <Dashboard dir> <out dir>`: each page dark, toggled, reloaded.

### Gene colours — the role convention must survive the move to five accents

Identity is carried by the stroke; fills are tints capped at 30% so a 15px arrow still reads.

| role | accent | treatment |
|---|---|---|
| `amr_anchor_mobile` acquired resistance | coral | 32% tint, 1.5px stroke |
| `amr_anchor_core` chromosomal backbone | coral | 12% tint, 1px stroke |
| `amr_anchor_point` point mutation | coral | no fill, 1.5px stroke, centre dot |
| `stress` | teal | 30% tint, 1px stroke |
| `virulence` | purple | 30% tint, 1px stroke |
| `mge` mobile element | mustard | 30% tint, 1px stroke |
| `mge` + pseudogene bit | mustard | 10% tint, dashed stroke |
| `other` / `hypothetical` | grey | `--bg-2` fill, `--border` stroke |

**Green is deliberately not a gene role.** It is reserved for pass and conserved states, so the
legend cannot be misread as having eight categories.

Two conflicts between the references, already resolved. The brand book forbids colour as fill, but a
pure-outline arrow is illegible at small widths, so tint plus solid stroke is the compromise. And
AlphaFold's canonical confidence colours are replaced by the Simple accents for consistency, with the
bucket thresholds captioned so domain readers are not misled: coral below 50, mustard 50–70, teal
70–90, `--fg` above 90.

## 6. Verification — passes as of 24 September 2026

Checked with `~/tools/shotter/dashcheck.mjs` (Playwright; screenshots read back; note that
`shot.mjs` cannot take a `file://` URL with a query string, so page states are driven from a script):

1. `run_dashboard.py` writes the six payloads and three pages in about a minute and
   `36_verify_payloads.py` passes: gene coordinates within 5 bp, contig order, cluster ids, roles,
   pseudogene flags, AMR classes of all 10,248 referenced clusters, chain coordinates within 0.06 Å,
   confidence within 0.6.
2. Each page opens from `file://` with no console errors; the only external hosts are
   `fonts.googleapis.com` and `fonts.gstatic.com`.
3. Biology: `blaCMY` at 80.5% prevalence with five class C beta-lactamase PDB matches at TM 1.00,
   and its 7,630 tier-3 windows aligned show the same backbone (lipocalin, sugE, ecnA, efp) in over
   85% of genomes; `CFRI046929` draws its two pseudogene transposases dashed beside `blaSHV-12`;
   `gyrA` marks T83 and D87 as catalog sites with T83 the most variable position; `clbB` renders as
   three modules with their true ranges; `blaKPC` carriers are meropenem-resistant in 4 of 6 against
   4 of 57 non-carriers; type 3 is 1,077 qacEΔ1–sul1 islands; carbapenem islands align on 22 anchor
   clusters led by blaKPC-2 at 43%.
4. Degenerate cases: a cluster without variation has no variation section, fragments and the five
   ColabFold structures are badged, an index-only cluster opens with its identity section alone, a
   tier-3 id opens with the *chromosomal window* badge while the default rail still shows
   46,921 of 66,005.
5. Caps: the island rail shows 300 with a note; a search for "protein" lists 400 of 2,187 with a
   note; aligned views draw 80 and count all; filtering resets the cap.
6. Keyboard: arrow keys move between tabs, rail rows and table rows are focusable, focus is shown
   by a teal border. Reduced motion disables spin. Single column with no horizontal scroll at 800 px.

Two lessons from the day, for whoever edits next: a text replacement that does not match fails
silently, so assert every one; and "no console errors" proved nothing three times (genes out of
order, a handler never replaced, a stale variable), so look at the screenshot.
