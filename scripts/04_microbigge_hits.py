"""Enrich the MicroBIGG-E per-hit table with the harmonised catalog annotation,
hit quality, SNP cluster and AMR gene co-location per contig; summarise protein variants.
"""
import pandas as pd

from amr_common import (ALL_ISOLATES_TSV, OUT, explode_genotypes, load_metadata,
                        load_microbigge, read_tsv, write_tsv)

HIT_QUALITY = {
    "EXACTP": "exact", "EXACTX": "exact", "ALLELEP": "exact", "ALLELEX": "exact",
    "BLASTP": "variant", "BLASTX": "variant",
    "POINTP": "point", "POINTX": "point",
    "PARTIALP": "partial", "PARTIALX": "partial",
    "PARTIAL_CONTIG_ENDP": "partial_contig_end", "PARTIAL_CONTIG_ENDX": "partial_contig_end",
    "HMM": "hmm_only", "INTERNAL_STOP": "internal_stop",
}


def main():
    mb = load_microbigge()
    ann = read_tsv(OUT / "symbol_annotation.tsv")
    clusters = read_tsv(ALL_ISOLATES_TSV)[["target_acc", "PDS_acc"]]

    mb = mb.rename(columns={"class": "class_as_called", "subclass": "subclass_as_called"})
    mb = mb.merge(ann[["symbol", "element_type", "gene_family", "class", "subclass",
                       "subclass_ambiguous", "refseq_protein_accession", "match_level"]]
                  .rename(columns={"refseq_protein_accession": "reference_protein_accession"}),
                  on=["symbol", "element_type"], how="left", validate="many_to_one")
    mb = mb.merge(clusters, on="target_acc", how="left", validate="many_to_one")
    mb["hit_quality"] = mb["method"].map(HIT_QUALITY).fillna("other")
    for col in ("start", "stop"):
        mb[col] = pd.to_numeric(mb[col])

    # AMR genes sharing a contig (plasmid / integron linkage)
    amr = mb[mb["element_type"] == "AMR"].sort_values(["target_acc", "contig", "start"])
    per_contig = (amr.groupby(["target_acc", "contig"])["symbol"]
                  .agg(n_amr_on_contig="size", amr_symbols_on_contig=";".join).reset_index())
    mb = mb.merge(per_contig, on=["target_acc", "contig"], how="left")
    mb["n_amr_on_contig"] = mb["n_amr_on_contig"].fillna(0).astype(int)
    mb["amr_symbols_on_contig"] = mb["amr_symbols_on_contig"].fillna("")

    mb = mb.sort_values(["target_acc", "contig", "start"]).reset_index(drop=True)
    mb = mb[["target_acc", "biosample_acc", "scientific_name", "PDS_acc", "contig", "start", "stop",
             "strand", "symbol", "element_name", "element_type", "subtype", "scope", "gene_family",
             "class", "subclass", "subclass_ambiguous", "class_as_called", "subclass_as_called",
             "method", "hit_quality", "pct_coverage", "pct_identity", "protein",
             "reference_protein_accession", "match_level", "n_amr_on_contig",
             "amr_symbols_on_contig"]]
    write_tsv(mb, "microbigge_hits.tsv")

    # gene/protein-level view. GenBank protein accessions are unique per assembly, so identical
    # proteins only share an accession when it is a RefSeq WP_ (non-redundant) one. Variants are
    # therefore summarised by their match signature against the reference protein.
    hits = mb.copy()
    hits["species"] = hits["scientific_name"].str.split().str[:2].str.join(" ")
    hits["refseq_protein"] = hits["protein"].where(hits["protein"].str.startswith("WP_"), "")
    distinct = lambda s: sorted(set(s) - {""})
    variants = (hits.groupby(["element_type", "symbol", "hit_quality", "pct_identity", "pct_coverage"])
                .agg(gene_family=("gene_family", "first"), class_=("class", "first"),
                     subclass=("subclass", "first"), reference_protein_accession=
                     ("reference_protein_accession", "first"),
                     n_hits=("target_acc", "size"), n_isolates=("target_acc", "nunique"),
                     n_clusters=("PDS_acc", lambda s: len(distinct(s))),
                     n_species=("species", "nunique"),
                     n_refseq_proteins=("refseq_protein", lambda s: len(distinct(s))),
                     refseq_proteins=("refseq_protein", lambda s: ";".join(distinct(s)[:10])),
                     example_protein=("protein", lambda s: next(iter(distinct(s)), "")))
                .reset_index().rename(columns={"class_": "class"}))
    for col in ("pct_identity", "pct_coverage"):
        variants[col] = pd.to_numeric(variants[col])
    variants = variants.sort_values(["element_type", "symbol", "n_isolates"],
                                    ascending=[True, True, False])
    write_tsv(variants, "gene_variant_summary.tsv")

    # QC 1: same element set per isolate as in the Pathogen Detection metadata.
    # The metadata lists each distinct symbol[=qualifier] once; MicroBIGG-E has one row per
    # copy, so sets must agree while MicroBIGG-E can have more rows (extra gene copies).
    md_long = explode_genotypes(load_metadata())
    md_long = md_long[md_long["target_acc"].isin(set(mb["target_acc"]))]
    for etype in ("AMR", "STRESS", "VIRULENCE"):
        a = md_long[md_long["element_type"] == etype].groupby("target_acc")["symbol"].agg(set)
        b = mb[mb["element_type"] == etype].groupby("target_acc")["symbol"].agg(set)
        a, b = a.align(b)
        same = sum(x == y for x, y in zip(a, b))
        extra = int((mb["element_type"] == etype).sum() - (md_long["element_type"] == etype).sum())
        print(f"QC {etype}: isolates with identical element sets in metadata and MicroBIGG-E: "
              f"{same:,} / {len(a):,}  (MicroBIGG-E has {extra:,} additional copies)")

    # QC 2: harmonised class vs the class MicroBIGG-E reported when the isolate was analysed
    amr = mb[mb["element_type"] == "AMR"]
    agree = (amr["class"] == amr["class_as_called"]).mean()
    agree_sub = (amr["subclass"] == amr["subclass_as_called"]).mean()
    print(f"QC AMR hits: class == class_as_called {agree:.2%}; subclass agreement {agree_sub:.2%}")
    diff = (amr[amr["subclass"] != amr["subclass_as_called"]]
            .groupby(["symbol", "class", "subclass", "class_as_called", "subclass_as_called"])
            .size().rename("n_hits").reset_index().sort_values("n_hits", ascending=False))
    write_tsv(diff, "qc_class_differences.tsv")
    print("\nhit quality (AMR):")
    print(amr["hit_quality"].value_counts().to_string())


if __name__ == "__main__":
    main()
