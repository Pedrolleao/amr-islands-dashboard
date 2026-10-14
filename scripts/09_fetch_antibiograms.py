"""Fetch the full antibiogram tables (MIC / disk values, method, standard) from NCBI BioSample.

Two searches, so nothing is missed:
  1. every BioSample of the dataset, by accession, restricted to records with an antibiogram
  2. the organism query of the group (groups.json), which also finds records outside the dataset

Output (<Organism_group>/Antibiograms/):
  biosample_antibiograms.xml        raw BioSample records as downloaded
  antibiogram_full.tsv              one row per isolate x antibiotic, isolates of this dataset
  antibiogram_isolates.tsv          one row per isolate with an antibiogram (genome? cluster? n tested)
  antibiogram_outside_dataset.tsv   same table for BioSamples found by the organism query but not in the dataset
"""
import re
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

import pandas as pd

from amr_common import ANTIBIOGRAM_QUERY, DATA, OUT, read_tsv

AST_OUT = DATA / "Antibiograms"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
PHENOTYPE_CODE = {"resistant": "R", "susceptible": "S", "intermediate": "I", "not defined": "ND",
                  "susceptible-dose dependent": "SDD", "nonsusceptible": "NS"}


def eutils(tool, **params):
    data = urllib.parse.urlencode(params).encode()
    for attempt in range(4):
        try:
            with urllib.request.urlopen(EUTILS + tool, data=data, timeout=180) as r:
                text = r.read().decode("utf-8")
            time.sleep(0.4)                       # NCBI allows 3 requests / second without a key
            return text
        except Exception as err:                  # transient NCBI errors: wait and retry
            print(f"  retry {tool} ({err})", flush=True)
            time.sleep(5 * (attempt + 1))
    raise SystemExit(f"NCBI {tool} kept failing")


def search_ids(term):
    xml = eutils("esearch.fcgi", db="biosample", term=term, retmax=10000)
    return re.findall(r"<Id>(\d+)</Id>", xml)


def snake(text):
    return re.sub(r"[^a-z0-9]+", "_", text.lower()).strip("_")


def parse_records(root):
    rows = []
    for sample in root.iter("BioSample"):
        table = next((t for t in sample.iter("Table") if t.get("class", "").startswith("Antibiogram")), None)
        if table is None:
            continue
        header = [snake(c.text or "") for c in table.find("Header")]
        organism = sample.find("Description/Organism")
        for row in table.find("Body"):
            cells = [(c.text or "").strip() for c in row]
            rows.append({"biosample_acc": sample.get("accession"),
                         "biosample_organism": organism.get("taxonomy_name", "") if organism is not None else "",
                         "antibiogram_version": table.get("class"), **dict(zip(header, cells))})
    return pd.DataFrame(rows)


def main():
    AST_OUT.mkdir(parents=True, exist_ok=True)
    profile = read_tsv(OUT / "isolate_resistance_profile.tsv")
    accessions = sorted(set(profile["biosample_acc"]) - {""})

    ids = set(search_ids(ANTIBIOGRAM_QUERY))
    print(f"BioSamples matching '{ANTIBIOGRAM_QUERY}': {len(ids)}")
    for i in range(0, len(accessions), 200):
        batch = " OR ".join(f"{a}[accn]" for a in accessions[i:i + 200])
        ids |= set(search_ids(f"antibiogram[filter] AND ({batch})"))
    print(f"after searching the {len(accessions):,} BioSamples of the dataset by accession: {len(ids)}")

    ids = sorted(ids, key=int)
    chunks = [eutils("efetch.fcgi", db="biosample", id=",".join(ids[i:i + 100]))
              for i in range(0, len(ids), 100)]
    root = ET.Element("BioSampleSet")
    for chunk in chunks:
        root.extend(ET.fromstring(chunk))
    ET.ElementTree(root).write(AST_OUT / "biosample_antibiograms.xml", encoding="utf-8",
                               xml_declaration=True)

    table = parse_records(root)
    if table.empty:
        print("no antibiogram tables found for this group")
        return
    table["result"] = table["resistance_phenotype"].str.lower().map(PHENOTYPE_CODE).fillna("")
    keys = profile[["biosample_acc", "target_acc", "asm_acc", "scientific_name", "PDS_acc",
                    "in_microbigge"]]
    merged = table.merge(keys, on="biosample_acc", how="left")
    first = ["target_acc", "biosample_acc", "asm_acc", "scientific_name", "PDS_acc", "antibiotic",
             "result", "resistance_phenotype", "measurement_sign", "measurement", "measurement_units"]
    rest = [c for c in merged.columns if c not in first + ["in_microbigge", "biosample_organism"]]

    ours = merged[merged["target_acc"].notna()].sort_values(["target_acc", "antibiotic"])
    ours[first + rest].to_csv(AST_OUT / "antibiogram_full.tsv", sep="\t", index=False)
    other = merged[merged["target_acc"].isna()].sort_values(["biosample_acc", "antibiotic"])
    other[["biosample_acc", "biosample_organism"] + first[5:] + rest].to_csv(
        AST_OUT / "antibiogram_outside_dataset.tsv", sep="\t", index=False)

    count = lambda code: (lambda s: int((s == code).sum()))
    isolates = (ours.groupby(["target_acc", "biosample_acc", "asm_acc", "scientific_name", "PDS_acc"])
                .agg(n_antibiotics=("antibiotic", "nunique"), n_resistant=("result", count("R")),
                     n_intermediate=("result", count("I")), n_susceptible=("result", count("S")),
                     has_measurements=("measurement", lambda s: int((s != "").any())),
                     methods=("laboratory_typing_method", lambda s: ";".join(sorted(set(s) - {""}))),
                     standards=("testing_standard", lambda s: ";".join(sorted(set(s) - {""}))))
                .reset_index())
    isolates.insert(3, "has_genome", (isolates["asm_acc"] != "").astype(int))
    isolates.to_csv(AST_OUT / "antibiogram_isolates.tsv", sep="\t", index=False)

    # QC against the S/I/R summary that came with the Pathogen Detection metadata
    summary = read_tsv(OUT / "ast_long.tsv")
    a = set(map(tuple, summary[["target_acc", "antibiotic", "result"]].values))
    b = set(map(tuple, ours[["target_acc", "antibiotic", "result"]].values))
    print(f"\nisolates of the dataset with an antibiogram: {isolates['target_acc'].nunique()} "
          f"({int(isolates['has_genome'].sum())} with a downloaded genome); "
          f"rows: {len(ours):,}; with a measured value: {int((ours['measurement'] != '').sum()):,}")
    print(f"QC metadata AST_phenotypes isolates covered: "
          f"{len(set(summary['target_acc']) & set(ours['target_acc']))} / {summary['target_acc'].nunique()}")
    print(f"QC isolate x antibiotic x result in both: {len(a & b):,}; only in metadata: {len(a - b)}; "
          f"only in BioSample: {len(b - a)}")
    print(f"BioSamples with an antibiogram outside the dataset: "
          f"{other['biosample_acc'].nunique()}")
    print("\ntyping methods:", ours["laboratory_typing_method"].value_counts().to_dict())
    print("testing standards:", ours["testing_standard"].value_counts().to_dict())


if __name__ == "__main__":
    main()
