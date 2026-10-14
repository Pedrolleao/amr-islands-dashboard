"""Catalog join for all isolates: one row per isolate x AMR/stress/virulence element.

Also parses AST_phenotypes (lab-measured susceptibility) into a long table.
"""
import pandas as pd

from amr_common import (COUNT_COLS, OUT, explode_genotypes, load_metadata, load_microbigge,
                        read_tsv, write_tsv)

ANNOTATION_COLS = ["subtype", "scope", "gene_family", "product_name", "class", "subclass",
                   "subclass_all", "subclass_ambiguous", "refseq_protein_accession",
                   "match_level", "db_version"]


def parse_ast(md):
    records = []
    for row in md[md["AST_phenotypes"] != ""].itertuples(index=False):
        for token in row.AST_phenotypes.split(","):
            antibiotic, _, result = token.strip().rpartition("=")
            if antibiotic:
                records.append((row.target_acc, antibiotic, result))
    return pd.DataFrame(records, columns=["target_acc", "antibiotic", "result"]).drop_duplicates()


def main():
    md = load_metadata()
    ann = read_tsv(OUT / "symbol_annotation.tsv")
    in_mb = set(load_microbigge()["target_acc"])

    long = explode_genotypes(md)
    long = long.merge(md[["target_acc", "biosample_acc", "scientific_name"]], on="target_acc")
    long = long.merge(ann[["symbol", "element_type"] + ANNOTATION_COLS],
                      on=["symbol", "element_type"], how="left", validate="many_to_one")
    long["in_microbigge"] = long["target_acc"].isin(in_mb).astype(int)
    long = long[["target_acc", "biosample_acc", "scientific_name", "element_type", "symbol",
                 "call_qualifier", "is_core"] + ANNOTATION_COLS + ["in_microbigge"]]
    write_tsv(long, "isolate_gene_long.tsv")
    write_tsv(parse_ast(md), "ast_long.tsv")

    # QC: element counts per isolate must equal NCBI's own number_*_genes columns
    n = long.groupby(["target_acc", "element_type"]).size().unstack(fill_value=0)
    n = n.reindex(md["target_acc"], fill_value=0)
    for etype, col in COUNT_COLS.items():
        expected = pd.to_numeric(md.set_index("target_acc")[col]).fillna(0).astype(int)
        bad = int((n[etype] != expected).sum())
        print(f"QC {etype}: isolates whose row count != {col}: {bad}")
    annotated = (long["match_level"] != "unmatched").mean()
    print(f"QC annotated element occurrences: {annotated:.2%}")
    print(f"isolates with >=1 element: {long['target_acc'].nunique():,} / {len(md):,}")


if __name__ == "__main__":
    main()
