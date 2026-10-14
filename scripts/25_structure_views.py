"""25 - py3Dmol views of every predicted structure, for the per-cluster dashboard.

For each cluster of Structures/structure_summary.tsv:
  views/html/<cluster>__<gene>.html   an interactive 3Dmol page, cartoon style, built with py3Dmol and
                                      then wired with a colour-mode switch. The coordinates are embedded,
                                      the library is loaded from views/lib/3Dmol-min.js, so the page works
                                      offline and can be dropped into a dashboard as an <iframe>.
  views/png/<cluster>__<gene>__<mode>.png   a still of the same scene per colour mode, rendered headless
                                      with Playwright through the viewer's own pngURI().

Colour modes (both read the B-factor column, which is why step 24 writes two PDB variants):
  plddt      AlphaFold confidence colours, from pdb/
  variation  share of lineages differing from the reference at that position, from pdb_variation/,
             a light-to-dark blue ramp matching the blue bars of the step-23 figures.

Usage: python3 scripts/25_structure_views.py [--group G] [--only CFR-000010,...] [--no-png] [--redo]
Needs the viz environment: tools/viz_env/bin/python (py3Dmol + playwright; Chromium already installed).
"""
import base64
import json
import re
import shutil
import sys
import time
from pathlib import Path

import pandas as pd
import py3Dmol

from amr_common import DATA, SET, STRUCTURES, read_tsv, set_dir, write_tsv

STRUCTURES = set_dir(STRUCTURES)

VIEWS = STRUCTURES / "views"
W, H = 900, 700

# pLDDT bands, the AlphaFold convention
PLDDT_BANDS = [(90, 100, "#0053d6"), (70, 90, "#65cbf3"), (50, 70, "#ffdb13"), (0, 50, "#ff7d45")]
# variation: one hue, light to dark (the dataviz skill's blue sequential ramp)
# grey for "no variation at all", then the blue ramp: a conserved position must not read as a faint hit
VAR_RAMP = ["#e2e1dc", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]
VAR_STOPS = [0, 1, 5, 10, 25, 50, 100]          # percent of lineages differing


def colour_js():
    """Two B-factor colour functions, shared by the page and the still renderer."""
    return f"""
const PLDDT_BANDS = {json.dumps(PLDDT_BANDS)};
const VAR_RAMP = {json.dumps(VAR_RAMP)};
const VAR_STOPS = {json.dumps(VAR_STOPS)};
function colourPlddt(atom) {{
  for (const [lo, hi, c] of PLDDT_BANDS) if (atom.b >= lo && atom.b <= hi) return c;
  return "#cccccc";
}}
function colourVariation(atom) {{
  let c = VAR_RAMP[0];
  for (let i = 0; i < VAR_STOPS.length; i++) if (atom.b >= VAR_STOPS[i]) c = VAR_RAMP[i];
  if (atom.b <= 0) c = VAR_RAMP[0];
  return c;
}}
"""


def page(cluster, title, subtitle, pdb_plddt, pdb_var, has_variation, mode="plddt"):
    """py3Dmol builds the viewer; the page then adds the colour switch and the still-render hook."""
    view = py3Dmol.view(width=W, height=H)
    view.addModel(pdb_plddt, "pdb")
    view.setStyle({"cartoon": {"color": "spectrum"}})
    view.zoomTo()
    base = view._make_html()
    # py3Dmol fetches 3Dmol from a CDN and keeps its viewer in a numbered global; use the local copy and
    # remember that name, so the page works offline and the colour switch can reach the viewer
    base = re.sub(r"https://[^'\"]*3Dmol-min\.js", "../lib/3Dmol-min.js", base)
    uuid = re.search(r"viewer_(\d+)", base).group(1)
    body = base[base.index("<div"):]

    legend_plddt = "".join(
        f'<span class="sw" style="background:{c}"></span>{lo}-{hi}' for lo, hi, c in PLDDT_BANDS)
    legend_var = "".join(
        f'<span class="sw" style="background:{VAR_RAMP[i]}"></span>'
        f'{"none" if VAR_STOPS[i] == 0 else str(VAR_STOPS[i]) + "%"}'
        for i in range(len(VAR_STOPS)))
    return f"""<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8"><title>{title}</title>
<script src="../lib/3Dmol-min.js"></script>
<style>
 :root {{ --surface:#fcfcfb; --ink:#1f1f1e; --ink2:#52514e; --line:#e6e5e0; }}
 body {{ margin:0; background:var(--surface); color:var(--ink);
        font:14px/1.5 -apple-system,BlinkMacSystemFont,"Segoe UI",Helvetica,Arial,sans-serif; }}
 header {{ padding:12px 16px 8px; }}
 h1 {{ font-size:15px; margin:0 0 2px; }}
 p.sub {{ margin:0; font-size:12px; color:var(--ink2); }}
 #wrap {{ position:relative; width:{W}px; height:{H}px; max-width:100%; margin:0 16px; }}
 .bar {{ display:flex; gap:12px; align-items:center; flex-wrap:wrap; padding:8px 16px; font-size:12px;
        color:var(--ink2); }}
 button {{ font:inherit; color:var(--ink); background:#fff; border:1px solid var(--line);
          border-radius:6px; padding:4px 10px; cursor:pointer; }}
 button[aria-pressed="true"] {{ background:var(--ink); color:#fff; border-color:var(--ink); }}
 .sw {{ display:inline-block; width:11px; height:11px; border-radius:2px; margin:0 4px 0 10px;
       vertical-align:-1px; }}
 .legend:first-child .sw:first-child {{ margin-left:0; }}
</style></head><body>
<header><h1>{title}</h1><p class="sub">{subtitle}</p></header>
<div class="bar">
  <button id="b-plddt" aria-pressed="true" onclick="setMode('plddt')">confidence (pLDDT)</button>
  <button id="b-variation" aria-pressed="false" {"" if has_variation else "disabled title='all copies identical'"}
          onclick="setMode('variation')">variation across lineages</button>
  <button onclick="toggleSpin()">spin</button>
  <span id="legend" class="legend"></span>
</div>
<div id="wrap">{body}</div>
<script>
{colour_js()}
const PDB = {{ plddt: {json.dumps(pdb_plddt)}, variation: {json.dumps(pdb_var)} }};
const LEGEND = {{ plddt: {json.dumps(legend_plddt)}, variation: {json.dumps(legend_var)} }};
const HAS_VARIATION = {json.dumps(bool(has_variation))};
const VIEWER_VAR = "viewer_{uuid}";
let viewer = null, spinning = false;
function setMode(mode) {{
  if (mode === "variation" && !HAS_VARIATION) return;
  viewer.clear();
  viewer.addModel(PDB[mode], "pdb");
  viewer.setStyle({{}}, {{cartoon: {{colorfunc: mode === "plddt" ? colourPlddt : colourVariation}}}});
  viewer.zoomTo();
  viewer.render();
  document.getElementById("legend").innerHTML =
    (mode === "plddt" ? "pLDDT " : "% of lineages differing ") + LEGEND[mode];
  for (const m of ["plddt", "variation"])
    document.getElementById("b-" + m).setAttribute("aria-pressed", String(m === mode));
  window.__mode = mode;
}}
function toggleSpin() {{ spinning = !spinning; viewer.spin(spinning ? "y" : false); }}
// py3Dmol builds the viewer inside a promise; wait for the global it assigns
const waitForViewer = setInterval(() => {{
  if (!window[VIEWER_VAR]) return;
  clearInterval(waitForViewer);
  viewer = window[VIEWER_VAR];
  setMode(new URLSearchParams(location.search).get("mode") || {json.dumps(mode)});
  window.__ready = true;
}}, 50);
</script></body></html>
"""


def render_stills(jobs, out_dir):
    """One headless browser for every still; the viewer renders its own PNG (no page chrome)."""
    from playwright.sync_api import sync_playwright
    done = []
    # this machine has an older cached Chromium than the playwright package expects; use what is there
    cached = sorted(Path.home().glob(".cache/ms-playwright/chromium-*/chrome-linux*/chrome"))
    with sync_playwright() as pw:
        browser = pw.chromium.launch(executable_path=str(cached[-1]) if cached else None,
                                     args=["--use-gl=angle", "--use-angle=swiftshader",
                                           "--enable-unsafe-swiftshader"])
        pg = browser.new_page(viewport={"width": W + 40, "height": H + 160})
        for html, name, mode in jobs:
            pg.goto(f"file://{html}?mode={mode}")
            pg.wait_for_function("window.__ready === true && window.__mode", timeout=60000)
            pg.wait_for_timeout(250)
            uri = pg.evaluate("window[VIEWER_VAR].pngURI()")
            (out_dir / f"{name}__{mode}.png").write_bytes(base64.b64decode(uri.split(",", 1)[1]))
            done.append((name, mode))
    return done


def main():
    t0 = time.time()
    only = set(sys.argv[sys.argv.index("--only") + 1].split(",")) if "--only" in sys.argv else None
    redo, no_png = "--redo" in sys.argv, "--no-png" in sys.argv

    s = read_tsv(STRUCTURES / "structure_summary.tsv")
    if only:
        s = s[s.cluster_id.isin(only)]
    (VIEWS / "html").mkdir(parents=True, exist_ok=True)
    (VIEWS / "png").mkdir(parents=True, exist_ok=True)
    lib = VIEWS / "lib" / "3Dmol-min.js"
    if not lib.exists():
        # every set keeps its own copy so a folder can be moved on its own; reuse one already downloaded
        have = sorted(DATA.glob("Structures*/views/lib/3Dmol-min.js"))
        if not have:
            raise SystemExit("views/lib/3Dmol-min.js is missing - download it from https://3Dmol.org/build/")
        lib.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy(have[0], lib)
        print(f"copied 3Dmol-min.js from {have[0].parent.parent.parent.name}")

    # clusters deferred to ColabFold have a summary row but no structure yet: nothing to show for them
    missing = [c.name for c in s.itertuples() if not (STRUCTURES / "pdb" / f"{c.name}.pdb").exists()]
    if missing:
        s = s[~s.name.isin(missing)]
        print(f"no structure yet for {len(missing)} cluster(s), skipped: {', '.join(sorted(missing))}")

    rows, jobs = [], []
    for c in s.itertuples():
        html_path = VIEWS / "html" / f"{c.name}.html"
        pdb_plddt = (STRUCTURES / "pdb" / f"{c.name}.pdb").read_text()
        pdb_var = (STRUCTURES / "pdb_variation" / f"{c.name}.pdb").read_text()
        has_var = float(c.max_variation_pct) > 0
        gene = c.family or c.product_top
        title = f"{c.cluster_id} · {gene} · {c.amr_classes or 'point mutation'}"
        ptm = "" if str(getattr(c, "ptm", "")) in ("", "nan", "None") else f", pTM {c.ptm}"
        sub = (f"{c.aa_len} aa · ESMFold, mean pLDDT {c.mean_plddt}{ptm}"
               f"{' · FRAGMENT of a longer gene' if c.is_fragment == 'True' else ''}"
               f" · {int(float(c.n_members)):,} protein copies in {c.n_distinct_seqs} distinct sequences"
               + (f" · most variable position differs in {c.max_variation_pct}% of lineages"
                  if has_var else " · every copy identical"))
        if redo or not html_path.exists():
            html_path.write_text(page(c.cluster_id, title, sub, pdb_plddt, pdb_var, has_var))
        modes = ["plddt"] + (["variation"] if has_var else [])
        for m in modes:
            png = VIEWS / "png" / f"{c.name}__{m}.png"
            if redo or not png.exists():
                jobs.append((str(html_path), c.name, m))
        rows.append({"cluster_id": c.cluster_id, "name": c.name, "gene": gene,
                     "html": f"views/html/{c.name}.html",
                     "png_plddt": f"views/png/{c.name}__plddt.png",
                     "png_variation": f"views/png/{c.name}__variation.png" if has_var else "",
                     "aa_len": c.aa_len, "mean_plddt": c.mean_plddt, "ptm": c.ptm,
                     "is_fragment": c.is_fragment, "has_variation": has_var,
                     "anchor_tier": c.anchor_tier, "amr_classes": c.amr_classes})
    print(f"pages: {len(rows)}; stills to render: {len(jobs)} ({time.time() - t0:.0f} s)", flush=True)

    if jobs and not no_png:
        batch = 50
        for i in range(0, len(jobs), batch):
            render_stills(jobs[i:i + batch], VIEWS / "png")
            print(f"  {min(i + batch, len(jobs))}/{len(jobs)} stills ({time.time() - t0:.0f} s)", flush=True)

    write_tsv(pd.DataFrame(rows), "view_manifest.tsv", VIEWS)
    n_html = len(list((VIEWS / "html").glob("*.html")))
    n_png = len(list((VIEWS / "png").glob("*.png")))
    print("\nchecks")
    print(f"  pages written: {n_html}/{len(s)}")
    print(f"  stills: {n_png} (expected {sum(1 for r in rows) + sum(1 for r in rows if r['has_variation'])})")
    empty = [p.name for p in (VIEWS / "png").glob("*.png") if p.stat().st_size < 5000]
    print(f"  stills under 5 kB (suspect blank): {len(empty)} {empty[:5]}")
    print(f"\ndone in {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
