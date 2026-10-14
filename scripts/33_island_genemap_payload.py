"""33 - the packed gene maps: every gene of every kept island, small enough to inline.

350,064 genes across 46,921 islands. Written as ordinary JSON - four integers and a string per gene -
this is about 7 MB before the cluster ids are even spelled out. Packed as below it is about 3 MB, and
the browser decodes only the island it is drawing.

THE PACKING, per gene, in `pos_in_island` order:

  gap    varint(zigzag(q(start) - q(previous_end)))          where q(x) = round(x / 10). Genes in an
                                                              operon abut or overlap, so the gap from
                                                              the last one is usually a single
                                                              character; an absolute coordinate never
                                                              is. The rounding is applied to the
                                                              POSITIONS and the deltas taken after,
                                                              never the other way round - rounding
                                                              deltas lets the error accumulate along
                                                              the island, and it reached 46 bp by the
                                                              end of one. This way every coordinate
                                                              is within 5 bp of the truth, always.
  width  varint(q(end) - q(start))
  cs     varint((cluster_index + 1) * 2 + strand_bit)         strand folded into the cluster
                                                              reference, so it costs no extra field;
                                                              0 and 1 mean "no cluster", strand -/+
  flags  varint: the draw role in the low 4 bits, bit 16 set when this gene's label differs from
                                                              its cluster's usual one, bit 32 set
                                                              for a pseudogene (a disrupted
                                                              transposase is drawn pale and dashed).
                                                              A varint and not a raw character: the
                                                              pseudogene bit is 32, which IS the
                                                              varint continuation bit, so a raw
                                                              character would make the stream
                                                              undecodable. It costs one extra
                                                              character on 10.7% of genes.
  label  varint(label_index), PRESENT ONLY when that bit is set

Labels are the reason this fits. Most clusters carry one label everywhere they appear, so the label
is stored once per cluster and only the exceptions - and the genes that have no cluster at all,
mostly pseudogenes - spend bytes on their own.

Coordinates are stored at 10 bp resolution. A gene map is drawn a few hundred pixels wide across tens
of kilobases, so 10 bp is far below one pixel.

Usage: python3 scripts/33_island_genemap_payload.py [--group G]
"""
import importlib.util
import os
import time
from collections import Counter, defaultdict

from amr_common import ISLANDS, PANGENOME, read_tsv

_spec = importlib.util.spec_from_file_location(
    "dash", os.path.join(os.path.dirname(os.path.abspath(__file__)), "30_dashboard_common.py"))
dash = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(dash)

KEEP_TIERS = ("tier1_multi_acquired", "tier2_single_acquired", "tier3_core_only")
BP = 10                      # coordinate resolution
LABEL_OVERRIDE = 16          # flag bit: this gene carries its own label
PSEUDO = 32                  # flag bit: pseudogene, drawn pale and dashed
COLS = ["island_id", "pos_in_island", "cluster_id", "start", "end", "strand", "draw_role",
        "is_pseudo", "label_text"]


def main():
    t0 = time.time()
    isl = read_tsv(ISLANDS / "islands.tsv", usecols=["island_id", "island_tier"])
    keep = set(isl.loc[isl.island_tier.isin(KEEP_TIERS), "island_id"])
    print(f"{len(keep):,} islands kept; reading the gene table", flush=True)

    g = read_tsv(ISLANDS / "island_gene_table.tsv", usecols=COLS)
    print(f"  gene table: {len(g):,} rows ({time.time() - t0:.0f} s)", flush=True)
    g = g[g.island_id.isin(keep)]
    # read_tsv gives strings, and "10" sorts before "2" as a string: the order must be numeric, or
    # the genes come out 0, 1, 10, 11, ... and every delta after the first is a long jump backwards
    g["pos_in_island"] = g.pos_in_island.astype(int)
    g = g.sort_values(["island_id", "pos_in_island"])
    print(f"  genes in kept islands: {len(g):,}", flush=True)

    g["cluster_id"] = g.cluster_id.fillna("")
    g["label_text"] = g.label_text.fillna("")
    g["draw_role"] = g.draw_role.fillna("other")

    roles = sorted(g.draw_role.unique())
    if len(roles) > 15:
        raise SystemExit(f"{len(roles)} draw roles will not fit the flag byte: {roles}")
    role_idx = {r: i for i, r in enumerate(roles)}

    # the usual label of each cluster: what a gene gets for free when it does not disagree
    label_counts = defaultdict(Counter)
    for cid, lab in zip(g.cluster_id, g.label_text):
        if cid:
            label_counts[cid][lab] += 1
    clusters = sorted(label_counts)
    cluster_idx = {c: i for i, c in enumerate(clusters)}
    default_label = {c: label_counts[c].most_common(1)[0][0] for c in clusters}
    varying = sum(1 for c in clusters if len(label_counts[c]) > 1)

    # the AMR class and subclass of each referenced cluster, so the page can align islands on
    # "whichever gene of this class the island carries" - the class view of the islands page
    ann = read_tsv(PANGENOME / "cluster_annotation.tsv",
                   usecols=["cluster_id", "amr_classes", "amr_subclasses"])
    ann = ann[ann.cluster_id.isin(cluster_idx)]
    cls_of = dict(zip(ann.cluster_id, ann.amr_classes.fillna("")))
    sub_of = dict(zip(ann.cluster_id, ann.amr_subclasses.fillna("")))

    d = dash.Dicts()
    cclass = [d.enc("amrclass", cls_of.get(c, "")) for c in clusters]
    csub = [d.enc("amrsubclass", sub_of.get(c, "")) for c in clusters]
    packed, overrides, n_genes, n_pseudo = {}, 0, 0, 0
    for island, grp in g.groupby("island_id", sort=False):
        parts, prev_end, first = [], None, None
        for r in grp.itertuples():
            qs, qe = round(int(r.start) / BP), round(int(r.end) / BP)
            if first is None:
                first = qs
                prev_end = qs
            parts.append(dash.vint(dash.zig(qs - prev_end)))
            parts.append(dash.vint(max(0, qe - qs)))
            ci = cluster_idx[r.cluster_id] + 1 if r.cluster_id else 0
            strand = 1 if str(r.strand) in ("1", "+", "+1") else 0
            parts.append(dash.vint(ci * 2 + strand))
            lab = str(r.label_text)
            override = (not r.cluster_id) or (lab != default_label.get(r.cluster_id, ""))
            flag = role_idx[r.draw_role] | (LABEL_OVERRIDE if override else 0)
            if str(r.is_pseudo) in ("1", "True", "1.0"):
                flag |= PSEUDO
                n_pseudo += 1
            parts.append(dash.vint(flag))
            if override:
                parts.append(dash.vint(d.enc("label", lab) + 1))   # +1 so "" is encodable as 0
                overrides += 1
            prev_end = qe
            n_genes += 1
        packed[dash.island_num(island)] = dash.vint(first) + "".join(parts)

    payload = {
        "built": time.strftime("%Y-%m-%d"), "bp": BP, "roles": roles,
        "flag_label": LABEL_OVERRIDE, "flag_pseudo": PSEUDO,
        "cluster": [dash.cluster_num(c) for c in clusters],
        "clabel": [default_label[c] for c in clusters],
        "cclass": cclass, "csub": csub,
        "dicts": d.as_dict(), "g": packed,
    }
    mb = dash.write_payload(payload, "genemaps.json")

    chars = sum(len(v) for v in packed.values())
    print("\nchecks")
    print(f"  islands packed: {len(packed):,}; genes: {n_genes:,} "
          f"({n_genes / max(1, len(packed)):.2f} per island)")
    print(f"  packed characters: {chars:,} ({chars / max(1, n_genes):.2f} per gene)")
    print(f"  label overrides: {overrides:,} of {n_genes:,} ({100 * overrides / n_genes:.1f}%); "
          f"{varying:,} of {len(clusters):,} clusters have more than one label")
    print(f"  draw roles ({len(roles)}): " + ", ".join(roles))
    print(f"  pseudogenes flagged: {n_pseudo:,} ({100 * n_pseudo / n_genes:.1f}%)")
    print(f"  clusters referenced: {len(clusters):,}; distinct override labels: {d.size('label'):,}")
    print(f"  clusters with an AMR class: {sum(1 for v in cclass if v >= 0):,} "
          f"({d.size('amrclass')} distinct class sets, {d.size('amrsubclass')} subclass sets)")
    print(f"  size {mb:.2f} MB")
    print(f"\ndone in {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
