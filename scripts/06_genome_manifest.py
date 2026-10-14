"""Decide which assembly accession to download for every isolate with a public genome.

MicroBIGG-E reports contigs / proteins either on the RefSeq copy of an assembly (NZ_ contigs,
WP_ proteins -> GCF_ accession) or on the GenBank copy (GCA_). The download accession follows
MicroBIGG-E so that its contig coordinates and protein accessions can be looked up directly.

Genomes/manifest/assembly_summary.jsonl is fetched with the NCBI datasets tool when missing.
"""
import json
import subprocess

import pandas as pd

from amr_common import GENOMES, OUT, TOOLS, read_tsv

MANIFEST = GENOMES / "manifest"


def main():
    profile = read_tsv(OUT / "isolate_resistance_profile.tsv")
    profile = profile.loc[profile["asm_acc"] != "",
                          ["target_acc", "biosample_acc", "scientific_name", "PDS_acc", "asm_acc"]]

    hits = read_tsv(OUT / "microbigge_hits.tsv", usecols=["target_acc", "contig", "protein"])
    hits["refseq"] = hits["contig"].str.match(r"^N[CZ]_")
    per_isolate = hits.groupby("target_acc").agg(
        microbigge_source=("refseq", lambda s: "RefSeq" if s.all() else "GenBank" if not s.any() else "mixed"),
        n_hits=("contig", "size"), n_hits_with_protein=("protein", lambda s: int((s != "").sum())))

    MANIFEST.mkdir(parents=True, exist_ok=True)
    summary_file = MANIFEST / "assembly_summary.jsonl"
    if not summary_file.exists():
        gca_list = MANIFEST / "gca_accessions.txt"
        gca_list.write_text("\n".join(profile["asm_acc"]) + "\n")
        with open(summary_file, "w") as out:
            subprocess.run([TOOLS / "datasets", "summary", "genome", "accession", "--inputfile",
                            gca_list, "--as-json-lines"], stdout=out, check=True)

    records = []
    with open(summary_file) as fh:
        for line in fh:
            r = json.loads(line)
            info, paired = r["assembly_info"], r["assembly_info"].get("paired_assembly", {})
            records.append({
                "asm_acc": r["accession"], "assembly_name": info.get("assembly_name", ""),
                "assembly_level": info.get("assembly_level", ""),
                "assembly_status": info.get("assembly_status", ""),
                "genbank_annotated": int("annotation_info" in r),
                "paired_refseq_acc": paired.get("accession", ""),
                "paired_refseq_status": paired.get("status", ""),
                "total_length": r.get("assembly_stats", {}).get("total_sequence_length", ""),
                "n_contigs": r.get("assembly_stats", {}).get("number_of_contigs", ""),
            })
    summary = pd.DataFrame(records)

    man = (profile.merge(per_isolate, on="target_acc", how="left")
           .merge(summary, on="asm_acc", how="left", validate="one_to_one"))
    use_refseq = (man["microbigge_source"] == "RefSeq") & (man["paired_refseq_status"] == "current")
    man["download_acc"] = man["asm_acc"].where(~use_refseq, man["paired_refseq_acc"])
    man["expect_proteins"] = ((use_refseq) | (man["genbank_annotated"] == 1)).astype(int)
    man.to_csv(MANIFEST / "genome_manifest.tsv", sep="\t", index=False)
    man["download_acc"].to_csv(MANIFEST / "download_accessions.txt", index=False, header=False)

    print(f"isolates with a public assembly: {len(man):,}")
    print(f"not found in NCBI summary: {int(man['assembly_name'].isna().sum())}")
    print("\nMicroBIGG-E source x download accession type:")
    print(pd.crosstab(man["microbigge_source"], man["download_acc"].str[:3]).to_string())
    print("\nassembly status:", man["assembly_status"].value_counts().to_dict())
    print("assembly level:", man["assembly_level"].value_counts().to_dict())
    print(f"expect a protein FASTA: {int(man['expect_proteins'].sum()):,}")
    odd = man[(man["microbigge_source"] == "RefSeq") & ~use_refseq]
    print(f"RefSeq in MicroBIGG-E but no current GCF: {len(odd)}")
    print(f"total genome size: {pd.to_numeric(man['total_length']).sum() / 1e9:.1f} Gbp")


if __name__ == "__main__":
    main()
