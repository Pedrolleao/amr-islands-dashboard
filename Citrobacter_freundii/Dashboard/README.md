# Dashboard – three self-contained pages over the finished analysis

Open any of the `.html` files straight from the folder. There is no server, no build step at
view time and no network use beyond the Google Fonts link (JetBrains Mono); everything each page
shows is spliced into the file, which is why they are megabytes each.

| page | what it answers | data inlined |
|---|---|---|
| `islands.html` | where a resistance gene sits: 46,921 tier-1 and tier-2 islands, gene maps in contig coordinates, architecture types, aligned views per type, cluster or AMR class, and a distribution tab; the 19,084 tier-3 backbone windows around core genes behind the tier filter | `islands.json`, `genemaps.json` |
| `clusters.html` | one record per analysed gene family (806: 411 resistance, 328 stress, 67 virulence) with its predicted fold drawn as a Cα trace, coloured by confidence or by where the family varies, plus variation, Pfam architecture, fold matches and islands; every other cluster of the pangenome searchable by product, architecture, AMR class and subclass | `clusters.json`, `cluster_index.json`, `structures.json` |
| `antibiograms.html` | the measured susceptibility NCBI holds for the group (75 isolates, 1,156 tests) beside each isolate's resistance genotype: resistance rates among carriers and non-carriers of every gene family per antibiotic, and one isolate at a time | `antibiograms.json` |

The pages link to each other through the address bar, so a link opens on a record:
`clusters.html?cluster=CFR-000010`, `islands.html?island=CFRI046929`,
`islands.html?cluster=CFR-000010` (every island carrying a cluster, aligned on it, with a conservation table),
`islands.html?amrclass=BETA-LACTAM` or `?subclass=CARBAPENEM` (every island with a gene of that class, aligned on that gene),
`islands.html?cluster=CFR-000010&tier=2` (the tier-3 backbone windows around a core gene; `tier=` takes `12` by default, `0`, `1`, `2`, or empty for all),
`islands.html?type=3` (every island of an architecture type, aligned, with a conservation table),
`clusters.html?set=stress`, `clusters.html?subclass=CARBAPENEM`, `islands.html?tab=distribution`,
`antibiograms.html?antibiotic=meropenem&family=blaKPC`, `antibiograms.html?isolate=PDT000050772.2`.
Tier-3 windows, the chromosomal context around a core resistance gene, are hidden by default and shown by the tier filter; they are labelled *chromosomal window*, not *island*.

## Rebuilding

```
python3 scripts/run_dashboard.py [--group G]            # payloads 31-34 and 37, pages 35, round-trip check 36
python3 scripts/run_dashboard.py --pages-only           # only 35, the loop when editing a template
```

The runner first checks that every upstream table exists and names the step to rerun if one is
missing. Its output is kept in `qc_dashboard_report.txt`. The page templates are
`scripts/dashboard/{islands,clusters,antibiograms}_template.html` with the shared `core.css` and
`core.js`; one copy serves every organism group. **Never edit a built page**: it is one
multi-megabyte line and the edit is lost on the next build.

## The payloads

| file | size | contents |
|---|---|---|
| `islands.json` | 8.0 MB | positional records for every island of all three tiers (`fields` gives the column order), their architecture types, the family co-occurrence edges, and the dictionaries every index points into |
| `genemaps.json` | 6.0 MB | the 695,638 genes of those islands, packed (see below), plus the cluster list, labels, AMR classes and subclasses they reference |
| `clusters.json` | 0.9 MB | the 806 analysed clusters: record, per-position variation, Pfam domains, Foldseek hits, up to 60 island ids each |
| `cluster_index.json` | 4.6 MB | all 208,757 clusters as parallel arrays, position = id: members, prevalence, product, architecture, class, set |
| `structures.json` | 1.7 MB | 808 Cα chains (clbB is three module chains) with sequence, coordinates, confidence and variation strings |
| `antibiograms.json` | 0.05 MB | the tested isolates with their calls, measurements, genotype and islands, and the antibiotic-to-subclass map that flags expected genes |

Cluster and island ids are dense integer ranges, so a payload stores the bare integer and the page
formats `CFR-000010` / `CFRI046929` back. Field semantics are documented in
`docs/WORKPLAN_dashboard.md` section 4.

## The codec (`scripts/30_dashboard_common.py`)

Three techniques keep the total at 16 MB: dictionary encoding of every repeated string, positional
records with the field order declared once, and base64-varint packing of the two high-cardinality
layers. The varint uses 5 bits per character over `A-Za-z0-9-_` with bit 32 as the continuation
flag; signed values are zigzagged; per-residue 0–100 values are one character each,
`chr(35 + round(v * 0.9))`.

A gene is four varints in island order: zigzag(gap from the previous gene's end), width, cluster
reference with the strand folded into the low bit, and flags (role in the low 4 bits, 16 = own
label follows, 32 = pseudogene). Coordinates are quantised to 10 bp *before* the deltas are taken,
so every position is within 5 bp of the annotation and the error never accumulates. A chain is the
zigzag deltas between consecutive Cα atoms at 0.1 Å, x/y/z interleaved.

The JavaScript half of the codec is generated by `js_decoder()` and spliced into every page by
step 35, so the two halves cannot drift. `36_verify_payloads.py` decodes the payloads again in
Python and compares them with the source tables and PDB files; run it after touching any encoder,
because a codec bug does not crash, it draws something subtly wrong. Three such bugs have been
caught this way: a raw flag character colliding with the continuation bit, deltas rounded instead
of positions, and genes sorted by position as text (`10` before `2`).

## Reading the colours

Colour is categorical metadata, never decoration. Gene roles in the maps: coral is resistance
(solid fill acquired, faint fill chromosomal, outline with a centre dot for a point-mutation
target), teal stress, purple virulence, mustard mobile element (dashed and pale when a pseudogene),
grey everything else. Green is reserved for pass and conserved states (catalog sites, same-fold
matches) and is never a gene role. Model confidence buckets: coral below 50, mustard 50–70, teal
70–90, white above 90; below 50 the model is not evidence. Variation buckets: none, under 5% of
lineages (teal), 5–20% (mustard), over 20% (coral).
