"""Download and pin the AMRFinderPlus Reference Gene Catalog + gene hierarchy (fam.tsv)."""
import shutil
import urllib.request

from amr_common import REF_BASE_URL, REF_FILES, REF_ROOT


def fetch(name, dest):
    with urllib.request.urlopen(REF_BASE_URL + name, timeout=300) as r, open(dest, "wb") as fh:
        shutil.copyfileobj(r, fh)


def main():
    with urllib.request.urlopen(REF_BASE_URL + "version.txt", timeout=60) as r:
        version = r.read().decode().strip()
    dest = REF_ROOT / f"AMRFinderPlus_{version}"
    dest.mkdir(parents=True, exist_ok=True)
    for name in REF_FILES:
        path = dest / name
        if path.exists() and path.stat().st_size > 0:
            print(f"have  {path.name}")
            continue
        fetch(name, path)
        print(f"fetched {path.name}  ({path.stat().st_size:,} bytes)")
    print(f"reference pinned at {dest}  (db version {version})")


if __name__ == "__main__":
    main()
