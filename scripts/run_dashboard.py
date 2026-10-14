"""Build the dashboard: payloads (31-34, 37), the three pages (35) and the round-trip check (36).

Usage: run_dashboard.py [--group <Organism_group>] [--pages-only]

Every input is an output of an earlier step, so the first thing done is to check that they exist and
say which step to rerun if one is missing. --pages-only skips 31-34 and 36 and just rebuilds the
HTML from the existing payloads, which is the loop when editing a template.
"""
import subprocess
import sys
from pathlib import Path

from amr_common import DASHBOARD, DOMAINS, GROUP, ISLANDS, PANGENOME, STRUCTURES, VARIATION

# (file, the step that writes it). Stress and virulence sets are optional: 31 and 34 skip a set
# whose summary is absent, so only the resistance set is required here.
REQUIRED = [
    (PANGENOME / "cluster_annotation.tsv", "15_cluster_annotation.py"),
    (ISLANDS / "islands.tsv", "17_define_islands.py"),
    (ISLANDS / "island_gene_table.tsv", "18_island_gene_table.py"),
    (ISLANDS / "gene_island.tsv", "18_island_gene_table.py"),
    (ISLANDS / "island_types.tsv", "19_island_types.py"),
    (ISLANDS / "island_type_membership.tsv", "19_island_types.py"),
    (ISLANDS / "amr_family_cooccurrence.tsv", "21_island_summary_figures.py"),
    (VARIATION / "position_variation.tsv", "23_cluster_variation.py"),
    (STRUCTURES / "structure_summary.tsv", "24_structure_prediction.py"),
    (DOMAINS / "cluster_domains.tsv", "27_pfam_domains.py"),
    (DOMAINS.with_name("Domains_all") / "cluster_domains.tsv", "27_pfam_domains.py --set all"),
]
PAYLOAD_STEPS = ["31_cluster_payload.py", "32_island_payload.py",
                 "33_island_genemap_payload.py", "34_structure_payload.py",
                 "37_antibiogram_payload.py"]      # 37 writes an empty payload if step 09 never ran
PAGE_STEP = "35_build_pages.py"
CHECK_STEP = "36_verify_payloads.py"


def main():
    here = Path(__file__).resolve().parent
    pages_only = "--pages-only" in sys.argv
    print(f"organism group: {GROUP}")
    missing = [(f, step) for f, step in REQUIRED if not f.exists()]
    if missing and not pages_only:
        for f, step in missing:
            print(f"  missing {f.relative_to(f.parents[2])}  <- rerun scripts/{step}")
        raise SystemExit(f"{len(missing)} upstream table(s) missing; nothing built")
    steps = ([PAGE_STEP] if pages_only else PAYLOAD_STEPS + [PAGE_STEP, CHECK_STEP])
    DASHBOARD.mkdir(parents=True, exist_ok=True)
    with open(DASHBOARD / "qc_dashboard_report.txt", "w" if not pages_only else "a") as report:
        for step in steps:
            header = f"\n===== {step} =====\n"
            print(header, end="", flush=True)
            report.write(header)
            run = subprocess.run([sys.executable, str(here / step), "--group", GROUP],
                                 capture_output=True, text=True)
            print(run.stdout, end="")
            report.write(run.stdout)
            if run.returncode:
                print(run.stderr, file=sys.stderr)
                report.write(run.stderr)
                raise SystemExit(f"{step} failed")
    print(f"\npages: {DASHBOARD}/{{islands,clusters,antibiograms}}.html"
          f"\nQC report: {DASHBOARD / 'qc_dashboard_report.txt'}")


if __name__ == "__main__":
    main()
