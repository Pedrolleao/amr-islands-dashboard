"""11 - gene index: gene ORDER and annotation of every CDS in every downloaded genome.

Nothing else in this pipeline reads `genomic.gff`. Steps 01-10 answer "which isolate carries which
resistance gene"; this file is what makes "where does that gene sit, and what surrounds it"
answerable, and it feeds the protein clustering (13) and the islands (17).

Per genome, in ONE pass over its two annotation files:
  genomic.gff  -> one row per CDS: contig, gene order along the contig, coordinates, strand,
                  locus_tag, protein_id, product, pseudogene / partial flags
  protein.faa  -> amino-acid length, md5 and sequence of each protein

THREE THINGS THAT LOOK LIKE DETAILS AND ARE NOT

1. PSEUDOGENE CDS ARE KEPT. ~2% of CDS carry `pseudo=true`, have no `protein_id` and appear in no
   FASTA - and they are disproportionately TRANSPOSASES (frameshifted Tn3, IS110, IS5 ...). Dropping
   them would shift every downstream gene ordinal AND systematically understate how much of an
   island is mobile-element machinery. They get a row, an ordinal, and `protein_md5 = ""`.

2. MULTI-SEGMENT CDS ARE COLLAPSED BY THEIR GFF `ID`. Programmed frameshifts (~0.09% of CDS) are
   written as two GFF lines for one gene. Counting lines instead of genes silently drifts the
   ordinals, which is the one number everything downstream depends on.

3. `aa_len` COMES FROM THE FASTA, NOT FROM THE COORDINATE SPAN. The span disagrees for ~0.8% of CDS
   (partial genes, frameshifts). `aa_len_gff` keeps the span-derived value so the two can be compared.

OUTPUT FORMAT. `gene_index.tsv.gz` is written as ONE GZIP MEMBER PER GENOME, concatenated. That is a
valid .gz for any reader, and because gzip members are independent, `gene_index_offsets.tsv` lets a
consumer decompress a single genome's block (~2 ms) instead of scanning 30M rows. There is no
columnar engine on this machine (no duckdb, no pyarrow), and making this pipeline depend on one when
its other ten scripts are pandas-only would be a regression.

DEDUPLICATION. Across 6,063 Citrobacter genomes there are ~30.8M CDS but only ~1.5M DISTINCT protein
sequences (measured: 150 genomes give 743,229 proteins / 188,156 distinct, ratio still falling).
`proteins_nr.faa` holds one record per distinct sequence, named by its md5 and written in md5-sorted
order so the clustering in 13 is reproducible. Clustering 1.5M sequences instead of 30.8M is the
difference between minutes and days, and loses nothing: 15 maps every copy back to its cluster.

Usage: python3 scripts/11_gene_index.py [--group G] [--workers N] [--limit N]
"""
import gzip
import hashlib
import shutil
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from urllib.parse import unquote

import pandas as pd

from amr_common import GENOMES, PANGENOME, read_tsv, write_tsv

WORKERS = 12          # I/O saturates around 8 workers on this drive; the work is light on CPU
CHUNK = 128           # genomes per shard - a static split, so a re-run is deterministic

COLUMNS = ["gene_id", "target_acc", "download_acc", "contig", "contig_ordinal", "ordinal",
           "start", "end", "strand", "locus_tag", "protein_id", "gene_symbol", "product",
           "aa_len", "aa_len_gff", "protein_md5", "is_pseudo", "is_partial", "n_segments",
           "contig_n_genes", "contig_len", "contig_topology", "contig_molecule"]
STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def clean(text):
    """GFF attribute value -> plain text (URL-encoded, and may hold tabs/newlines)."""
    return unquote(text).replace("\t", " ").replace("\n", " ").replace("\r", " ").strip()


def attrs(field):
    out = {}
    for item in field.rstrip().split(";"):
        if "=" in item:
            k, _, v = item.partition("=")
            out[k] = v
    return out


def iter_fasta(path):
    """(accession, sequence) per record; accession is the first token of the header."""
    acc, chunks = None, []
    try:
        fh = open(path)
    except FileNotFoundError:
        return
    with fh:
        for line in fh:
            if line.startswith(">"):
                if acc is not None:
                    yield acc, "".join(chunks)
                acc = line[1:].split(None, 1)[0]
                chunks = []
            else:
                chunks.append(line.strip())
        if acc is not None:
            yield acc, "".join(chunks)


def parse_gff(path):
    """-> (genes in first-seen order, contigs). One entry per CDS gene, segments collapsed by ID."""
    genes, order, contigs = {}, [], {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            seqid, ftype, start, end, strand, attr = f[0], f[2], f[3], f[4], f[6], f[8]
            if ftype == "region":
                a = attrs(attr)
                if seqid not in contigs:
                    contigs[seqid] = {
                        "contig_len": end, "n_genes": 0,
                        "topology": "circular" if a.get("Is_circular") == "true" else "linear",
                        "molecule": a.get("genome", ""),
                        "plasmid_name": clean(a.get("plasmid-name", "")),
                        "seen": len(contigs) + 1}
                continue
            if ftype != "CDS":
                continue
            a = attrs(attr)
            gid = a.get("ID", "")
            if gid in genes:                       # second segment of a frameshifted CDS
                g = genes[gid]
                g["start"] = min(g["start"], int(start))
                g["end"] = max(g["end"], int(end))
                g["n_segments"] += 1
                continue
            if seqid not in contigs:               # CDS before its region line: still index it
                contigs[seqid] = {"contig_len": "", "n_genes": 0, "topology": "", "molecule": "",
                                  "plasmid_name": "", "seen": len(contigs) + 1}
            genes[gid] = {
                "contig": seqid, "start": int(start), "end": int(end), "strand": strand,
                "locus_tag": a.get("locus_tag", ""), "protein_id": a.get("protein_id", ""),
                "gene_symbol": clean(a.get("gene", "")), "product": clean(a.get("product", "")),
                "is_pseudo": 1 if a.get("pseudo") == "true" else 0,
                "is_partial": 1 if (a.get("partial") == "true" or "start_range" in a
                                    or "end_range" in a) else 0,
                "n_segments": 1}
            order.append(gid)
    return [genes[g] for g in order], contigs


def process_genome(target_acc, acc):
    """Parse one genome -> (row strings, per-contig rows, per-md5 info, counters)."""
    folder = GENOMES / "assemblies" / acc
    gff = folder / "genomic.gff"
    if not gff.exists():
        return None

    genes, contigs = parse_gff(gff)
    seq_of_acc, md5_of_acc = {}, {}
    for pacc, seq in iter_fasta(folder / "protein.faa"):
        md5 = hashlib.md5(seq.encode()).hexdigest()
        md5_of_acc[pacc] = (len(seq), md5)
        if md5 not in seq_of_acc:
            seq_of_acc[md5] = (seq, pacc)

    by_contig = {}
    for g in genes:
        by_contig.setdefault(g["contig"], []).append(g)
    for contig, lst in by_contig.items():
        lst.sort(key=lambda g: (g["start"], g["end"]))      # gene order along the contig
        for i, g in enumerate(lst, 1):
            g["ordinal"] = i
        contigs[contig]["n_genes"] = len(lst)

    rows, copies = [], {}
    st = {"n_cds": 0, "n_pseudo": 0, "n_with_protein": 0, "n_seg2": 0,
          "n_len_ok": 0, "n_len_cmp": 0, "n_contigs": len(by_contig)}
    for contig in sorted(by_contig, key=lambda c: contigs[c]["seen"]):
        cinfo = contigs[contig]
        for g in by_contig[contig]:
            aa_len, md5 = md5_of_acc.get(g["protein_id"], (0, "")) if g["protein_id"] else (0, "")
            aa_len_gff = max((g["end"] - g["start"] + 1) // 3 - 1, 0)
            if md5:
                st["n_with_protein"] += 1
                copies[md5] = copies.get(md5, 0) + 1
                st["n_len_cmp"] += 1
                st["n_len_ok"] += int(aa_len == aa_len_gff)
            st["n_pseudo"] += g["is_pseudo"]
            st["n_seg2"] += int(g["n_segments"] > 1)
            rows.append("\t".join(map(str, [
                f'{acc}|{g["locus_tag"]}', target_acc, acc, contig, cinfo["seen"], g["ordinal"],
                g["start"], g["end"], g["strand"], g["locus_tag"], g["protein_id"],
                g["gene_symbol"], g["product"], aa_len, aa_len_gff, md5, g["is_pseudo"],
                g["is_partial"], g["n_segments"], cinfo["n_genes"], cinfo["contig_len"],
                cinfo["topology"], cinfo["molecule"]])))
    st["n_cds"] = len(rows)
    crows = [(acc, target_acc, c, contigs[c]["seen"], contigs[c]["contig_len"],
              contigs[c]["n_genes"], contigs[c]["topology"], contigs[c]["molecule"],
              contigs[c]["plasmid_name"]) for c in by_contig]
    return rows, crows, copies, seq_of_acc, st


def run_chunk(args):
    """One worker, one slice of genomes: its own gzip shard, nr-sequence shard and count shard.

    Genomes are partitioned across workers, so summing each worker's "distinct genomes holding this
    md5" over workers gives the exact global genome count - no global set is ever needed.
    """
    idx, jobs, shard_dir = args
    shard_dir = Path(shard_dir)
    gz_path = shard_dir / f"genes_{idx:04d}.tsv.gz"
    seq_path = shard_dir / f"seqs_{idx:04d}.tsv"
    cnt_path = shard_dir / f"counts_{idx:04d}.tsv"
    offsets, crows, missing, stats = [], [], [], []
    copies, genomes = {}, {}
    written = set()
    with open(gz_path, "wb") as gz_out, open(seq_path, "w") as seq_out:
        for target_acc, acc in jobs:
            res = process_genome(target_acc, acc)
            if res is None:
                missing.append((acc, target_acc, "genomic.gff"))
                continue
            rows, c_rows, g_copies, seq_of_md5, st = res
            start_off = gz_out.tell()
            comp = gzip.compress(("\n".join(rows) + "\n").encode(), 6)
            gz_out.write(comp)                       # one independent gzip member per genome
            offsets.append((acc, start_off, len(comp), st["n_cds"]))
            crows.extend(c_rows)
            for md5, n in g_copies.items():
                copies[md5] = copies.get(md5, 0) + n
                genomes[md5] = genomes.get(md5, 0) + 1
                if md5 not in written:
                    written.add(md5)
                    seq, pacc = seq_of_md5[md5]
                    seq_out.write(f"{md5}\t{pacc}\t{acc}\t{seq}\n")
            st.update(acc=acc, target_acc=target_acc)
            stats.append(st)
    with open(cnt_path, "w") as fh:
        for md5, n in copies.items():
            fh.write(f"{md5}\t{n}\t{genomes[md5]}\n")
    return {"idx": idx, "gz": str(gz_path), "seqs": str(seq_path), "counts": str(cnt_path),
            "offsets": offsets, "crows": crows, "missing": missing, "stats": stats}


def main():
    argv = sys.argv[1:]
    workers = int(argv[argv.index("--workers") + 1]) if "--workers" in argv else WORKERS
    limit = int(argv[argv.index("--limit") + 1]) if "--limit" in argv else None

    man = read_tsv(GENOMES / "manifest" / "genome_manifest.tsv")
    jobs = [(r.target_acc, r.download_acc) for r in man.itertuples()
            if r.download_acc and (GENOMES / "assemblies" / r.download_acc).is_dir()]
    jobs.sort(key=lambda j: j[1])
    if limit:
        jobs = jobs[:limit]
    log(f"{len(jobs):,} genomes, {workers} workers")

    PANGENOME.mkdir(parents=True, exist_ok=True)
    shard_dir = PANGENOME / ".shards"
    if shard_dir.exists():
        shutil.rmtree(shard_dir)
    shard_dir.mkdir()

    chunks = [(i, jobs[s:s + CHUNK], str(shard_dir))
              for i, s in enumerate(range(0, len(jobs), CHUNK))]
    t0 = time.time()
    results = []
    with ProcessPoolExecutor(max_workers=workers) as ex:
        for n, res in enumerate(ex.map(run_chunk, chunks), 1):
            results.append(res)
            if n % 8 == 0 or n == len(chunks):
                done = sum(len(r["stats"]) for r in results)
                log(f"  {done:,}/{len(jobs):,} genomes ({time.time()-t0:.0f}s)")
    results.sort(key=lambda r: r["idx"])

    # ---- merge gzip shards byte-wise (gzip members concatenate) -------------------------------
    log("merging shards")
    index_path = PANGENOME / "gene_index.tsv.gz"
    offsets = []
    with open(index_path, "wb") as out:
        out.write(gzip.compress(("\t".join(COLUMNS) + "\n").encode(), 6))
        for res in results:
            base = out.tell()
            with open(res["gz"], "rb") as fh:
                shutil.copyfileobj(fh, out, 1 << 20)
            offsets.extend((acc, base + off, length, n) for acc, off, length, n in res["offsets"])
    write_tsv(pd.DataFrame(offsets, columns=["download_acc", "byte_offset", "byte_length", "n_cds"]),
              "gene_index_offsets.tsv", folder=PANGENOME)

    # ---- per-genome / per-contig summaries -----------------------------------------------------
    src_of = dict(zip(man["download_acc"], man["microbigge_source"]))
    lvl_of = dict(zip(man["download_acc"], man["assembly_level"]))
    grows, crows, missing = [], [], []
    tot = {k: 0 for k in ("n_cds", "n_pseudo", "n_with_protein", "n_seg2", "n_len_ok", "n_len_cmp")}
    for res in results:
        missing.extend(res["missing"])
        crows.extend(res["crows"])
        for st in res["stats"]:
            grows.append((st["acc"], st["target_acc"], src_of.get(st["acc"], ""),
                          lvl_of.get(st["acc"], ""), st["n_contigs"], st["n_cds"], st["n_pseudo"],
                          st["n_with_protein"]))
            for k in tot:
                tot[k] += st[k]
    write_tsv(pd.DataFrame(grows, columns=["download_acc", "target_acc", "source", "assembly_level",
                                           "n_contigs", "n_cds", "n_pseudo", "n_cds_with_protein"]),
              "genome_index.tsv", folder=PANGENOME)
    write_tsv(pd.DataFrame(crows, columns=["download_acc", "target_acc", "contig", "contig_ordinal",
                                           "contig_len", "n_genes", "topology", "molecule",
                                           "plasmid_name"]),
              "contig_index.tsv", folder=PANGENOME)

    # ---- non-redundant proteome ----------------------------------------------------------------
    log("building the non-redundant proteome")
    copies, genomes = {}, {}
    for res in results:
        with open(res["counts"]) as fh:
            for line in fh:
                md5, n, g = line.rstrip("\n").split("\t")
                copies[md5] = copies.get(md5, 0) + int(n)
                genomes[md5] = genomes.get(md5, 0) + int(g)
    seqs, example = {}, {}
    for res in results:
        with open(res["seqs"]) as fh:
            for line in fh:
                md5, pacc, acc, seq = line.rstrip("\n").split("\t")
                if md5 not in seqs:
                    seqs[md5] = seq
                    example[md5] = (pacc, acc)

    nr_path = PANGENOME / "proteins_nr.faa"
    odd = 0
    with open(nr_path, "w") as out:
        for md5 in sorted(seqs):                      # md5 order -> reproducible clustering input
            seq = seqs[md5]
            if set(seq) - STANDARD_AA:
                odd += 1
            out.write(f">{md5}\n")
            out.writelines(seq[i:i + 80] + "\n" for i in range(0, len(seq), 80))
    write_tsv(pd.DataFrame([(m, len(seqs[m]), copies[m], genomes[m], example[m][0], example[m][1])
                            for m in sorted(seqs)],
                           columns=["protein_md5", "aa_len", "n_copies", "n_genomes",
                                    "example_protein_id", "example_genome"]),
              "protein_nr.tsv", folder=PANGENOME)

    # ---- QC -------------------------------------------------------------------------------------
    n_cds, n_prot = tot["n_cds"], tot["n_with_protein"]
    print("\n===== QC 11_gene_index =====")
    print(f"genomes indexed: {len(grows):,}   annotation missing: {len(missing)}")
    for acc, target, what in missing[:10]:
        print(f"  MISSING {what}: {acc} ({target})")
    print(f"CDS rows: {n_cds:,}   with a protein: {n_prot:,}   "
          f"pseudogenes: {tot['n_pseudo']:,} ({100*tot['n_pseudo']/max(n_cds,1):.2f}%)")
    print(f"multi-segment CDS (programmed frameshift): {tot['n_seg2']:,} "
          f"({100*tot['n_seg2']/max(n_cds,1):.3f}%)  [expect ~0.09%; 0 means ID collapse is broken]")
    print(f"aa_len == span-derived length: {100*tot['n_len_ok']/max(tot['n_len_cmp'],1):.2f}%  "
          f"[expect ~99.2%]")
    print(f"distinct protein sequences: {len(seqs):,}  "
          f"(dedup ratio {len(seqs)/max(n_prot,1):.4f} of {n_prot:,} copies)")
    print(f"sum(n_copies) == CDS with a protein: {sum(copies.values()) == n_prot}")
    print(f"proteins with a non-standard residue: {odd}")
    print(f"contigs: {len(crows):,}")
    print(f"gene_index.tsv.gz {index_path.stat().st_size/1e9:.2f} GB   "
          f"proteins_nr.faa {nr_path.stat().st_size/1e9:.2f} GB")
    shutil.rmtree(shard_dir)
    log(f"done in {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
