"""Build symbol -> antibiotic class lookup from the Reference Gene Catalog (+ fam.tsv hierarchy).

One row per distinct (symbol, element_type) seen in the metadata genotype columns.
Match cascade (recorded in match_level):
  allele -> gene_family -> hierarchy_node -> hierarchy_symbol -> microbigge_only -> unmatched
"""
import pandas as pd

from amr_common import (TAXGROUP, db_version, explode_genotypes, load_metadata,
                        load_microbigge, read_tsv, ref_dir, split_classes, write_tsv)

COLUMNS = ["symbol", "element_type", "subtype", "scope", "gene_family", "product_name",
           "class", "subclass", "subclass_all", "subclass_ambiguous",
           "refseq_protein_accession", "pubmed_reference", "match_level", "db_version",
           "n_isolates"]


def consensus(rows, col):
    """Most common non-empty value of a column ('' if none)."""
    vals = rows[col][rows[col] != ""]
    return vals.mode().iat[0] if len(vals) else ""


def union_subclasses(rows):
    return "/".join(sorted({s for v in rows["subclass"] for s in split_classes(v)}))


def n_pairs(rows):
    pairs = rows[["class", "subclass"]].drop_duplicates()
    return len(pairs[(pairs["class"] != "") | (pairs["subclass"] != "")])


def annotate(symbol, etype, cat, fam, mb_cls):
    rec = dict.fromkeys(COLUMNS, "")
    rec.update(symbol=symbol, element_type=etype, subclass_ambiguous=0)

    typed = cat[cat["type"] == etype]
    alleles = typed[typed["allele"] == symbol]
    if alleles["subtype"].str.startswith("POINT").any():
        own = alleles[alleles["whitelisted_taxa"] == TAXGROUP]
        alleles = own if len(own) else alleles
    family = typed[typed["gene_family"] == symbol]
    node = fam[fam["#node_id"] == symbol]
    by_symbol = fam[(fam["gene_symbol"] == symbol) & (fam["type"] == etype)]
    # several nodes can share a gene symbol (e.g. eae intimin types): keep the top-most one(s)
    by_symbol = by_symbol[~by_symbol["parent_node_id"].isin(by_symbol["#node_id"])]

    if len(alleles):
        r = alleles.iloc[0]
        rec.update(match_level="allele", subtype=r["subtype"], scope=r["scope"],
                   gene_family=r["gene_family"], product_name=r["product_name"],
                   refseq_protein_accession=r["refseq_protein_accession"],
                   pubmed_reference=r["pubmed_reference"],
                   subclass_all=union_subclasses(alleles),
                   subclass_ambiguous=int(n_pairs(alleles) > 1))
        rec["class"], rec["subclass"] = r["class"], r["subclass"]
    elif len(family):
        # family-level call: AMRFinderPlus reports the hierarchy node's class/subclass
        src = node if len(node) else by_symbol if len(by_symbol) else family
        rec.update(match_level="gene_family", subtype=consensus(family, "subtype"),
                   scope=consensus(family, "scope"), gene_family=symbol,
                   product_name=(consensus(src, "family_name") if "family_name" in src
                                 else consensus(family, "product_name")),
                   subclass_all=union_subclasses(family),
                   subclass_ambiguous=int(n_pairs(family) > 1))
        rec["class"], rec["subclass"] = consensus(src, "class"), consensus(src, "subclass")
    elif len(node) or len(by_symbol):
        src = node if len(node) else by_symbol
        rec.update(match_level="hierarchy_node" if len(node) else "hierarchy_symbol",
                   subtype=consensus(src, "subtype"), gene_family=symbol,
                   product_name=consensus(src, "family_name"),
                   subclass_all=union_subclasses(src),
                   subclass_ambiguous=int(n_pairs(src) > 1))
        rec["class"], rec["subclass"] = consensus(src, "class"), consensus(src, "subclass")
    elif (symbol, etype) in mb_cls.index:
        r = mb_cls.loc[(symbol, etype)]
        rec.update(match_level="microbigge_only", subtype=r["subtype"], scope=r["scope"],
                   product_name=r["element_name"], subclass_all=r["subclass"])
        rec["class"], rec["subclass"] = r["class"], r["subclass"]
    else:
        rec["match_level"] = "unmatched"
    return rec


def main():
    cat = read_tsv(ref_dir() / "ReferenceGeneCatalog.txt")
    fam = read_tsv(ref_dir() / "fam.tsv")
    version = db_version()

    long = explode_genotypes(load_metadata())
    counts = (long.groupby(["symbol", "element_type"])["target_acc"].nunique()
              .rename("n_isolates").reset_index())

    # most frequent MicroBIGG-E annotation per symbol, used only as a last resort
    mb = load_microbigge()
    mb_cls = (mb.groupby(["symbol", "element_type"])
              [["subtype", "scope", "element_name", "class", "subclass"]]
              .agg(lambda s: s.mode().iat[0]))

    recs = [annotate(r.symbol, r.element_type, cat, fam, mb_cls)
            for r in counts.itertuples(index=False)]
    ann = pd.DataFrame(recs, columns=COLUMNS)
    ann["db_version"] = version
    ann["n_isolates"] = counts["n_isolates"].values
    ann = ann.sort_values(["element_type", "symbol"]).reset_index(drop=True)

    write_tsv(ann, "symbol_annotation.tsv")
    write_tsv(ann.loc[ann["match_level"] == "unmatched", ["symbol", "element_type", "n_isolates"]],
              "unmatched_symbols.tsv")

    print("\nmatch levels:")
    print(ann.groupby(["element_type", "match_level"]).size().unstack(fill_value=0))
    no_class = ann[(ann["element_type"] == "AMR") & (ann["class"] == "")]
    print(f"\nAMR symbols without a class: {len(no_class)}  {no_class['symbol'].tolist()}")
    print(f"family-level symbols with >1 subclass among alleles: {int(ann['subclass_ambiguous'].sum())}")


if __name__ == "__main__":
    main()
