"""The three groups of protein clusters that steps 23-25 can run on.

  amr        the 411 resistance anchors of Islands/anchor_clusters.tsv - unchanged, so the resistance
             output stays reproducible
  stress     metal / biocide resistance: element_type == STRESS in cluster_annotation.tsv
  virulence  element_type == VIRULENCE

Every set is returned in the column shape the resistance run already uses, so steps 23-25 need no
knowledge of which group they are looking at:
  cluster_id, rep_md5, n_seqs, n_members, prevalence, product_top,
  anchor_tier, amr_symbols, amr_gene_families, amr_classes, point_symbols
For stress and virulence the symbol columns are filled from their own columns and `point_symbols` is
empty: neither group has catalog point mutations, so the wild-type-reference rule of step 23 is a no-op.
"""
from amr_common import ISLANDS, PANGENOME, SET, read_tsv

COLUMNS = ["cluster_id", "rep_md5", "n_seqs", "n_members", "prevalence", "product_top",
           "anchor_tier", "amr_symbols", "amr_gene_families", "amr_classes", "point_symbols"]

ANNOTATION_COLS = ["cluster_id", "rep_md5", "n_seqs", "n_members", "prevalence", "product_top",
                   "element_type", "stress_symbols", "virulence_symbols", "aa_len_mean"]


def load_cluster_set(name=None):
    name = name or SET
    if name == "all":
        # every cluster of the pangenome - for the domain annotation of steps 27-28 only
        d = read_tsv(PANGENOME / "cluster_annotation.tsv", usecols=ANNOTATION_COLS + ["anchor_tier"])
        sym = d.stress_symbols.fillna("").where(d.stress_symbols.notna(), "")
        d["amr_symbols"] = sym.where(sym != "", d.virulence_symbols.fillna(""))
        d["amr_gene_families"] = d.amr_symbols.fillna("").str.split(";").str[0]
        d["amr_classes"] = d.element_type.fillna("")
        d["anchor_tier"] = d.anchor_tier.fillna("NOT_AN_ANCHOR")
        d["point_symbols"] = ""
    elif name == "amr":
        anchors = read_tsv(ISLANDS / "anchor_clusters.tsv")
        ann = read_tsv(PANGENOME / "cluster_annotation.tsv",
                       usecols=["cluster_id", "rep_md5", "n_seqs", "point_symbols"])
        d = anchors.merge(ann, on="cluster_id")
    else:
        d = read_tsv(PANGENOME / "cluster_annotation.tsv", usecols=ANNOTATION_COLS)
        d = d[d.element_type == name.upper()].copy()
        symbols = d.stress_symbols if name == "stress" else d.virulence_symbols
        d["amr_symbols"] = symbols
        # the "family" label used in figure titles and file names: the first symbol of the cluster
        d["amr_gene_families"] = symbols.fillna("").str.split(";").str[0]
        d["amr_classes"] = d.element_type
        d["anchor_tier"] = "NOT_AN_ANCHOR"       # stress and virulence never anchor an island
        d["point_symbols"] = ""
    d["n_seqs"] = d.n_seqs.astype(int)
    d["n_members"] = d.n_members.astype(int)
    return d[COLUMNS].sort_values("cluster_id").reset_index(drop=True)
