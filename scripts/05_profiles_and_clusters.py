"""Per-isolate resistance profile (all isolates) and per-SNP-cluster AMR summary."""
import pandas as pd

from amr_common import (ALL_ISOLATES_TSV, KEY_SUBCLASSES, MARKER_COLUMN, MARKER_SUBCLASS, OUT,
                        SPOT_CHECKS, load_metadata, read_tsv, split_classes, write_tsv)

META_COLS = ["target_acc", "biosample_acc", "asm_acc", "scientific_name", "collection_date",
             "geo_loc_name", "host", "isolation_source", "epi_type"]
INCOMPLETE = {"PARTIAL", "PARTIAL_END_OF_CONTIG", "MISTRANSLATION"}


def colname(prefix, value):
    return f"{prefix}_{value.lower().replace(' ', '_').replace('-', '_')}"


def join_sorted(values):
    return ";".join(sorted(set(values)))


def presence_matrix(amr, col, prefix, keep=None):
    """0/1 isolate x class matrix; multi-valued entries like 'A/B' count for both A and B."""
    pairs = amr[["target_acc", col]].drop_duplicates()
    pairs = pairs.assign(value=pairs[col].map(split_classes)).explode("value").dropna(subset="value")
    if keep is not None:
        pairs = pairs[pairs["value"].isin(keep)]
    mat = pd.crosstab(pairs["target_acc"], pairs["value"]).clip(upper=1)
    return mat.rename(columns=lambda v: colname(prefix, v))


def main():
    md = load_metadata()
    long = read_tsv(OUT / "isolate_gene_long.tsv")
    hits = read_tsv(OUT / "microbigge_hits.tsv")
    clusters = read_tsv(ALL_ISOLATES_TSV)[["target_acc", "PDS_acc"]]
    ann = read_tsv(OUT / "symbol_annotation.tsv").set_index("symbol")

    checked = [s for s in SPOT_CHECKS if s in ann.index]          # genes absent from a group are skipped
    for symbol in checked:
        got = list(ann.loc[symbol, ["class", "subclass"]])
        assert got == list(SPOT_CHECKS[symbol]), f"spot check failed for {symbol}: {got} != {SPOT_CHECKS[symbol]}"
    print(f"QC spot checks passed: {', '.join(checked)}")

    amr = long[long["element_type"] == "AMR"]
    amr_hits = hits[hits["element_type"] == "AMR"]

    # marker genes (MARKER_SUBCLASS, e.g. carbapenemases for Enterobacterales): harmonised subclass,
    # plus the subclass MicroBIGG-E reported for the actual hit (more specific than the family-level
    # value for e.g. blaOXA variants)
    # MARKER_SUBCLASS may be null (groups with no high-concern acquired determinant): no marker column
    if MARKER_SUBCLASS:
        is_carba = lambda s: s.str.split("/").map(lambda v: MARKER_SUBCLASS in v)
        carba = pd.concat([amr.loc[is_carba(amr["subclass"]), ["target_acc", "symbol"]],
                           amr_hits.loc[is_carba(amr_hits["subclass_as_called"]), ["target_acc", "symbol"]]])
    else:
        carba = pd.DataFrame(columns=["target_acc", "symbol"], dtype=str)

    profile = md[META_COLS].merge(clusters, on="target_acc", how="left", validate="one_to_one")
    profile["in_microbigge"] = profile["target_acc"].isin(set(hits["target_acc"])).astype(int)
    per_isolate = pd.DataFrame({
        "n_amr_genes": amr.groupby("target_acc").size(),
        "n_amr_incomplete": amr[amr["call_qualifier"].isin(INCOMPLETE)].groupby("target_acc").size(),
        "n_amr_point_mutations": amr[amr["subtype"] == "POINT"].groupby("target_acc").size(),
        "amr_genes": amr.groupby("target_acc")["symbol"].agg(join_sorted),
        **({MARKER_COLUMN: carba.groupby("target_acc")["symbol"].agg(join_sorted)}
           if MARKER_SUBCLASS else {}),
        # copy-level information only exists for isolates in MicroBIGG-E
        "n_amr_hits_microbigge": amr_hits.groupby("target_acc").size(),
        "n_contigs_with_amr": amr_hits.groupby("target_acc")["contig"].nunique(),
        "max_amr_on_one_contig": pd.to_numeric(amr_hits["n_amr_on_contig"])
                                   .groupby(amr_hits["target_acc"]).max(),
    })
    classes = presence_matrix(amr, "class", "class")
    subclasses = presence_matrix(amr, "subclass", "subclass", keep=KEY_SUBCLASSES)
    if MARKER_SUBCLASS:
        marker_col = colname("subclass", MARKER_SUBCLASS)
        subclasses.loc[subclasses.index.intersection(carba["target_acc"]), marker_col] = 1
        for t in set(carba["target_acc"]) - set(subclasses.index):
            subclasses.loc[t] = 0
            subclasses.loc[t, marker_col] = 1

    profile = (profile.set_index("target_acc").join(per_isolate).join(classes).join(subclasses)
               .reset_index())
    class_cols = list(classes.columns) + list(subclasses.columns)
    count_cols = ["n_amr_genes", "n_amr_incomplete", "n_amr_point_mutations"]
    profile[count_cols + class_cols] = profile[count_cols + class_cols].fillna(0).astype(int)
    for col in ("amr_genes", *( (MARKER_COLUMN,) if MARKER_SUBCLASS else () )):
        profile[col] = profile[col].fillna("")
    mb_only = ["n_amr_hits_microbigge", "n_contigs_with_amr", "max_amr_on_one_contig"]
    profile[mb_only] = profile[mb_only].astype("Int64")
    profile.loc[profile["in_microbigge"] == 1, mb_only] = \
        profile.loc[profile["in_microbigge"] == 1, mb_only].fillna(0)
    profile["n_amr_classes"] = profile[list(classes.columns)].sum(axis=1)
    class_lists = (amr["class"].map(split_classes).groupby(amr["target_acc"])
                   .agg(lambda s: join_sorted(c for v in s for c in v)))
    profile["amr_classes"] = profile["target_acc"].map(class_lists).fillna("")
    write_tsv(profile, "isolate_resistance_profile.tsv")
    assert len(profile) == len(md) and profile["target_acc"].is_unique

    # ---- per SNP cluster -------------------------------------------------------------------
    members = profile[profile["PDS_acc"] != ""].copy()
    members["year"] = pd.to_numeric(members["collection_date"].str[:4], errors="coerce")
    members["country"] = members["geo_loc_name"].str.split(":").str[0].str.strip()
    members["species"] = members["scientific_name"].str.split().str[:2].str.join(" ")
    nonempty = lambda s: join_sorted(v for v in s if v)
    summary = members.groupby("PDS_acc").agg(
        n_isolates=("target_acc", "size"), species=("species", nonempty),
        year_min=("year", "min"), year_max=("year", "max"), countries=("country", nonempty),
        median_amr_genes=("n_amr_genes", "median"))
    summary[["year_min", "year_max"]] = summary[["year_min", "year_max"]].astype("Int64")
    pct = (members.groupby("PDS_acc")[class_cols].mean() * 100).round(1).add_prefix("pct_")

    gene_counts = (amr.drop_duplicates(["target_acc", "symbol"])
                   .merge(members[["target_acc", "PDS_acc"]], on="target_acc")
                   .groupby(["PDS_acc", "symbol"]).size().rename("n").reset_index()
                   .merge(summary["n_isolates"].reset_index(), on="PDS_acc"))
    gene_counts["label"] = gene_counts["symbol"] + "(" + gene_counts["n"].astype(str) + ")"
    common = (gene_counts[gene_counts["n"] >= 0.5 * gene_counts["n_isolates"]]
              .groupby("PDS_acc")["label"].agg(";".join).rename("genes_in_half_or_more"))
    summary = summary.join(common)
    if MARKER_SUBCLASS:
        carba_counts = (carba.drop_duplicates().merge(members[["target_acc", "PDS_acc"]], on="target_acc")
                        .groupby(["PDS_acc", "symbol"]).size().rename("n").reset_index())
        carba_counts["label"] = carba_counts["symbol"] + "(" + carba_counts["n"].astype(str) + ")"
        summary = summary.join(carba_counts.groupby("PDS_acc")["label"].agg(";".join).rename(MARKER_COLUMN))
    summary = (summary.join(pct).fillna({
        "genes_in_half_or_more": "", **({MARKER_COLUMN: ""} if MARKER_SUBCLASS else {})})
        .sort_values("n_isolates", ascending=False).reset_index())
    write_tsv(summary, "cluster_amr_summary.tsv")

    # ---- overview ---------------------------------------------------------------------------
    print(f"\nisolates: {len(profile):,}; in an SNP cluster: {len(members):,} "
          f"({summary.shape[0]:,} clusters); in MicroBIGG-E: {int(profile['in_microbigge'].sum()):,}")
    print("\n% of isolates carrying >=1 gene per class:")
    print((profile[class_cols].mean() * 100).round(1).sort_values(ascending=False).to_string())
    if MARKER_SUBCLASS:
        top = carba.drop_duplicates().groupby("symbol").size().sort_values(ascending=False).head(10)
        print(f"\ntop {MARKER_SUBCLASS.lower()} genes -> {MARKER_COLUMN} (isolates):")
        print(top.to_string())


if __name__ == "__main__":
    main()
