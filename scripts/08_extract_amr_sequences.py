"""Pull the AMR / stress / virulence sequences of every isolate out of its downloaded genome.

For each MicroBIGG-E hit:
  protein     looked up by accession in the isolate's protein.faa
  nucleotide  cut from the genome by contig + start/stop/strand (also covers hits that have
              no protein accession, e.g. genes with an internal stop or found by BLASTX only)

Output (<Organism_group>/AMR_sequences/):
  proteins/<TYPE>/<gene_family>.faa         every copy, header = isolate|protein|symbol|quality|location
  proteins_unique/<TYPE>/<gene_family>.faa  identical sequences collapsed, with isolate counts
  nucleotide/<TYPE>/<gene_family>.fna       gene nucleotide sequences
  sequence_index.tsv                        one row per hit: found?, lengths, sequence md5
  protein_variants.tsv                      one row per distinct protein sequence of each symbol
Usage: 08_extract_amr_sequences.py [assemblies_dir]   (default Genomes/assemblies)
"""
import hashlib
import re
import shutil
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

from amr_common import DATA, GENOMES, OUT, read_tsv

SEQ_OUT = DATA / "AMR_sequences"
COMPLEMENT = str.maketrans("ACGTRYKMSWBDHVNacgtrykmswbdhvn", "TGCAYRMKSWVHDBNtgcayrmkswvhdbn")
BASES = "TCAG"
CODON_TABLE = dict(zip((a + b + c for a in BASES for b in BASES for c in BASES),
                       "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"))


def read_fasta(path, wanted=None):
    """{id: sequence}; only ids in `wanted` when given."""
    seqs, name, keep = {}, None, False
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].split(None, 1)[0]
                keep = wanted is None or name in wanted
                if keep:
                    seqs[name] = []
            elif keep:
                seqs[name].append(line.strip())
    return {k: "".join(v) for k, v in seqs.items()}


def translate(nt):
    nt = nt.upper()
    return "".join(CODON_TABLE.get(nt[i:i + 3], "X") for i in range(0, len(nt) - 2, 3))


def md5(seq):
    return hashlib.md5(seq.encode()).hexdigest()


def safe_name(text):
    return re.sub(r"[^A-Za-z0-9._-]+", "_", text).strip("_") or "unnamed"


def write_fastas(records, folder, ext):
    """records: {(type, family): [(header, seq), ...]}"""
    if folder.exists():
        shutil.rmtree(folder)
    used = {}                     # file names must stay distinct on case-insensitive disks too
    for (etype, family), items in sorted(records.items()):
        name = safe_name(family)
        n = 1
        while used.get((etype, name.lower()), family) != family:
            n += 1
            name = f"{safe_name(family)}_{n}"
        used[(etype, name.lower())] = family
        path = folder / etype / f"{name}.{ext}"
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, "w") as fh:
            for header, seq in items:
                fh.write(f">{header}\n")
                fh.writelines(seq[i:i + 80] + "\n" for i in range(0, len(seq), 80))


def main():
    assemblies = Path(sys.argv[1]) if len(sys.argv) > 1 else GENOMES / "assemblies"
    man = read_tsv(GENOMES / "manifest" / "genome_manifest.tsv")
    hits = read_tsv(OUT / "microbigge_hits.tsv")
    hits["family"] = hits["gene_family"].where(hits["gene_family"] != "", hits["symbol"])
    hits[["start", "stop"]] = hits[["start", "stop"]].astype(int)
    acc_of = dict(zip(man["target_acc"], man["download_acc"]))

    proteins, nucleotides = defaultdict(list), defaultdict(list)
    index = []
    n_done = 0
    for target, grp in hits.groupby("target_acc", sort=False):
        folder = assemblies / acc_of[target]
        fna = next(folder.glob("*_genomic.fna"), None) if folder.is_dir() else None
        if fna is None:
            continue
        n_done += 1
        if n_done % 500 == 0:
            print(f"  {n_done:,} genomes processed", flush=True)
        contigs = read_fasta(fna, set(grp["contig"]))
        faa = folder / "protein.faa"
        prots = read_fasta(faa, set(grp["protein"]) - {""}) if faa.exists() else {}

        for h in grp.itertuples(index=False):
            key = (h.element_type, h.family)
            location = f"{h.contig}:{h.start}-{h.stop}({h.strand})"
            nt = contigs.get(h.contig, "")[h.start - 1:h.stop]
            if nt and h.strand == "-":
                nt = nt.translate(COMPLEMENT)[::-1]
            prot = prots.get(h.protein, "")
            base = f"{target}|{h.protein or 'no_protein_acc'}|{h.symbol}|{h.hit_quality}|{location}"
            if nt:
                nucleotides[key].append((base, nt))
            if prot:
                proteins[key].append((base, prot))
            # does the annotated protein agree with the translated coordinates?
            agrees = ""
            if nt and prot:
                agrees = int(translate(nt).rstrip("*")[1:] == prot[1:])
            index.append((target, acc_of[target], h.contig, h.start, h.stop, h.strand, h.element_type,
                          h.symbol, h.family, h.hit_quality, h.protein, len(nt), len(prot),
                          md5(prot) if prot else "", md5(nt.upper()) if nt else "", agrees))

    idx = pd.DataFrame(index, columns=[
        "target_acc", "assembly_acc", "contig", "start", "stop", "strand", "element_type", "symbol",
        "gene_family", "hit_quality", "protein", "nt_length", "protein_length", "protein_md5",
        "nt_md5", "translation_matches_protein"])
    SEQ_OUT.mkdir(parents=True, exist_ok=True)
    idx.to_csv(SEQ_OUT / "sequence_index.tsv", sep="\t", index=False)
    write_fastas(proteins, SEQ_OUT / "proteins", "faa")
    write_fastas(nucleotides, SEQ_OUT / "nucleotide", "fna")

    # distinct protein sequences: identical proteins share an md5 whatever their accession
    meta = hits[["target_acc", "PDS_acc", "scientific_name"]].drop_duplicates("target_acc")
    with_prot = idx[idx["protein_md5"] != ""].merge(meta, on="target_acc")
    with_prot["species"] = with_prot["scientific_name"].str.split().str[:2].str.join(" ")
    refseq_first = lambda s: next((p for p in sorted(set(s)) if p.startswith("WP_")), sorted(set(s))[0])
    variants = (with_prot.groupby(["element_type", "gene_family", "symbol", "protein_md5"])
                .agg(protein_length=("protein_length", "first"), n_copies=("target_acc", "size"),
                     n_isolates=("target_acc", "nunique"),
                     n_clusters=("PDS_acc", lambda s: len(set(s) - {""})),
                     n_species=("species", "nunique"),
                     hit_quality=("hit_quality", lambda s: ";".join(sorted(set(s)))),
                     representative_protein=("protein", refseq_first),
                     example_isolate=("target_acc", "first"))
                .reset_index()
                .sort_values(["element_type", "gene_family", "n_isolates"], ascending=[True, True, False]))
    variants["variant_id"] = (variants["gene_family"].map(safe_name) + "_"
                              + variants["protein_md5"].str[:8])
    variants.to_csv(SEQ_OUT / "protein_variants.tsv", sep="\t", index=False)

    seq_of = {}
    for items in proteins.values():
        for _, seq in items:
            seq_of.setdefault(md5(seq), seq)
    unique = defaultdict(list)
    for v in variants.drop_duplicates(["element_type", "gene_family", "protein_md5"]).itertuples(index=False):
        same = variants[(variants["gene_family"] == v.gene_family) & (variants["protein_md5"] == v.protein_md5)]
        header = (f"{v.variant_id}|{';'.join(sorted(set(same['symbol'])))}|len={v.protein_length}"
                  f"|n_isolates={int(same['n_isolates'].sum())}|{v.representative_protein}")
        unique[(v.element_type, v.gene_family)].append((header, seq_of[v.protein_md5]))
    write_fastas(unique, SEQ_OUT / "proteins_unique", "faa")

    # ---- QC -----------------------------------------------------------------------------------
    expected = hits[hits["target_acc"].isin(set(idx["target_acc"]))]
    print(f"genomes processed: {n_done:,} / {hits['target_acc'].nunique():,}")
    print(f"QC nucleotide sequences cut: {(idx['nt_length'] > 0).sum():,} / {len(expected):,} hits")
    span_ok = (idx["nt_length"] == idx["stop"] - idx["start"] + 1).mean()
    print(f"QC nucleotide length == stop-start+1: {span_ok:.2%}")
    has_acc = expected["protein"] != ""
    print(f"QC proteins found by accession: {(idx['protein_length'] > 0).sum():,} / {int(has_acc.sum()):,}")
    checked = idx[idx["translation_matches_protein"] != ""]
    print(f"QC translated coordinates == annotated protein: "
          f"{(checked['translation_matches_protein'] == 1).mean():.2%} of {len(checked):,}")
    print(f"distinct protein sequences: {variants['protein_md5'].nunique():,} "
          f"(from {len(with_prot):,} protein copies); FASTA files: "
          f"{sum(1 for _ in (SEQ_OUT / 'proteins').rglob('*.faa')):,} gene families")


if __name__ == "__main__":
    main()
