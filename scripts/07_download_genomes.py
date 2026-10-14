"""Download genome FASTA + protein FASTA + GFF3 for every accession in the genome manifest.

Uses NCBI datasets "dehydrated" packages (the recommended route for thousands of genomes) in
batches. Resumable: accessions already complete in Genomes/assemblies/ are skipped, so the
script can simply be re-run after an interruption.

Result:  Genomes/assemblies/<accession>/{<accession>_<name>_genomic.fna, protein.faa, genomic.gff}
"""
import shutil
import subprocess
import sys
import zipfile

from amr_common import GENOMES, TOOLS, read_tsv

BATCH_SIZE = 500
WORKERS = "10"
ASSEMBLIES = GENOMES / "assemblies"
STAGING = GENOMES / "_staging"
MANIFEST = GENOMES / "manifest"


def is_complete(acc, expect_proteins):
    d = ASSEMBLIES / acc
    if not any(d.glob("*_genomic.fna")):
        return False
    return not expect_proteins or ((d / "protein.faa").exists() and (d / "genomic.gff").exists())


def run(cmd):
    print("$", " ".join(str(c) for c in cmd), flush=True)
    return subprocess.run(cmd, stdout=sys.stdout, stderr=subprocess.STDOUT).returncode


def download_batch(accessions, label):
    if STAGING.exists():
        shutil.rmtree(STAGING)
    STAGING.mkdir(parents=True)
    acc_file = STAGING / "accessions.txt"
    acc_file.write_text("\n".join(accessions) + "\n")
    package = STAGING / "package.zip"
    if run([TOOLS / "datasets", "download", "genome", "accession", "--inputfile", acc_file,
            "--include", "genome,protein,gff3", "--dehydrated", "--filename", package,
            "--no-progressbar"]):
        return
    with zipfile.ZipFile(package) as z:
        z.extractall(STAGING)
    run([TOOLS / "datasets", "rehydrate", "--directory", STAGING, "--max-workers", WORKERS,
         "--no-progressbar"])

    data = STAGING / "ncbi_dataset" / "data"
    report = data / "assembly_data_report.jsonl"
    if report.exists():
        with open(MANIFEST / "assembly_data_report.jsonl", "a") as out:
            out.write(report.read_text())
    moved = 0
    for d in data.iterdir():
        if d.is_dir() and any(d.glob("*_genomic.fna")):
            target = ASSEMBLIES / d.name
            if target.exists():
                shutil.rmtree(target)
            shutil.move(str(d), str(target))
            moved += 1
    print(f"[{label}] stored {moved} / {len(accessions)} assemblies", flush=True)
    shutil.rmtree(STAGING)


def main():
    man = read_tsv(MANIFEST / "genome_manifest.tsv")
    ASSEMBLIES.mkdir(parents=True, exist_ok=True)
    expect = dict(zip(man["download_acc"], man["expect_proteins"] == "1"))

    for attempt in (1, 2):                      # second pass retries whatever failed
        todo = [a for a in man["download_acc"] if not is_complete(a, expect[a])]
        print(f"pass {attempt}: {len(man) - len(todo):,} complete, {len(todo):,} to download", flush=True)
        if not todo:
            break
        batches = [todo[i:i + BATCH_SIZE] for i in range(0, len(todo), BATCH_SIZE)]
        for i, batch in enumerate(batches, 1):
            download_batch(batch, f"pass {attempt} batch {i}/{len(batches)}")

    missing = [a for a in man["download_acc"] if not is_complete(a, expect[a])]
    (MANIFEST / "download_missing.txt").write_text("\n".join(missing) + ("\n" if missing else ""))
    print(f"done: {len(man) - len(missing):,} complete, {len(missing):,} missing "
          f"(listed in manifest/download_missing.txt)", flush=True)


if __name__ == "__main__":
    main()
