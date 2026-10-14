"""Run the resistance-island pipeline (steps 11-22) and the per-cluster variation profiles (23) and save every QC block to one report.

Usage: run_islands.py [--group <Organism_group>] [--from N] [--skip-clustering]

Step 13 (MMseqs2) is the only step that needs a large scratch filesystem; --skip-clustering resumes
from an existing clusters_raw.tsv.
"""
import subprocess
import sys
from pathlib import Path

from amr_common import GROUP, ISLANDS

STEPS = ["11_gene_index.py", "12_link_hits_to_genes.py", "13_protein_clusters.py",
         "14_cluster_registry.py", "15_cluster_annotation.py", "16_amr_anchors.py",
         "17_define_islands.py", "18_island_gene_table.py", "19_island_types.py",
         "20_island_figures.py", "21_island_summary_figures.py", "22_island_qc.py",
         "23_cluster_variation.py"]


def main():
    here = Path(__file__).resolve().parent
    start = int(sys.argv[sys.argv.index("--from") + 1]) if "--from" in sys.argv else 11
    steps = [s for s in STEPS if int(s[:2]) >= start]
    if "--skip-clustering" in sys.argv:
        steps = [s for s in steps if not s.startswith("13")]
    print(f"organism group: {GROUP}\nsteps: {', '.join(s[:2] for s in steps)}")
    ISLANDS.mkdir(parents=True, exist_ok=True)
    with open(ISLANDS / "qc_islands_report.txt", "w") as report:
        for step in steps:
            header = f"\n===== {step} =====\n"
            print(header, end="", flush=True)
            report.write(header)
            run = subprocess.run([sys.executable, str(here / step), "--group", GROUP],
                                 capture_output=True, text=True)
            print(run.stdout, end="")
            report.write(run.stdout)
            # 22 exits non-zero when a check fails; that is a result, not a crash
            if run.returncode and not step.startswith("22"):
                print(run.stderr, file=sys.stderr)
                report.write(run.stderr)
                raise SystemExit(f"{step} failed")
    print(f"\nQC report: {ISLANDS / 'qc_islands_report.txt'}")


if __name__ == "__main__":
    main()
