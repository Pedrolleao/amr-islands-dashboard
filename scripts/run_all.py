"""Run the whole AMR annotation pipeline and save every script's QC output to qc_report.txt.

Usage: run_all.py [--group <Organism_group>] [--sequences]
"""
import subprocess
import sys
from pathlib import Path

from amr_common import GROUP, OUT

STEPS = ["01_fetch_reference.py", "02_symbol_annotation.py", "03_isolate_gene_long.py",
         "04_microbigge_hits.py", "05_profiles_and_clusters.py"]
# genomes + sequences (large download, resumable): run_all.py --sequences
SEQUENCE_STEPS = ["06_genome_manifest.py", "07_download_genomes.py", "08_extract_amr_sequences.py"]
# written last so that it sees every table produced in this run
FINAL_STEPS = ["10_group_readme.py"]


def main():
    here = Path(__file__).resolve().parent
    print(f"organism group: {GROUP}")
    OUT.mkdir(parents=True, exist_ok=True)
    with open(OUT / "qc_report.txt", "w") as report:
        for step in STEPS + (SEQUENCE_STEPS if "--sequences" in sys.argv else []) + FINAL_STEPS:
            header = f"\n===== {step} =====\n"
            print(header, end="")
            report.write(header)
            run = subprocess.run([sys.executable, str(here / step)], capture_output=True, text=True)
            print(run.stdout, end="")
            report.write(run.stdout)
            if run.returncode:
                print(run.stderr, file=sys.stderr)
                report.write(run.stderr)
                raise SystemExit(f"{step} failed")
    print(f"\nQC report: {OUT / 'qc_report.txt'}")


if __name__ == "__main__":
    main()
