"""Build the browser version of SeqGuide (a single static web page).

    python web/build.py            -> site/index.html (+ .nojekyll)

The page runs SeqGuide with stlite (Streamlit on Pyodide/WebAssembly) entirely in the visitor's
browser: no server, no installation, data never leaves the computer. External tools are not
available there, so SeqGuide uses its built-in engine automatically.
"""

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
STLITE = "0.85.1"   # pinned: stlite's API can change between versions

sys.path.insert(0, str(ROOT))
from seqguide import __version__  # noqa: E402

files = {"app/streamlit_app.py": (ROOT / "app" / "streamlit_app.py").read_text()}
for p in sorted((ROOT / "seqguide").glob("*.py")):
    files[f"seqguide/{p.name}"] = p.read_text()

html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>SeqGuide</title>
<link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/@stlite/browser@{STLITE}/build/stlite.css">
<style>
  body {{ margin: 0; font-family: -apple-system, Segoe UI, Roboto, Helvetica, Arial, sans-serif; }}
  .splash {{ max-width: 560px; margin: 18vh auto; padding: 0 24px; color: #333; line-height: 1.6; }}
  .splash h1 {{ font-size: 26px; margin-bottom: 4px; }}
  .bar {{ height: 4px; background: #e8e8e8; border-radius: 2px; overflow: hidden; margin-top: 18px; }}
  .bar div {{ width: 30%; height: 100%; background: #2a78d6; animation: slide 1.4s infinite ease-in-out; }}
  @keyframes slide {{ 0% {{ margin-left: -30%; }} 100% {{ margin-left: 100%; }} }}
</style>
</head>
<body>
<div id="root">
  <div class="splash">
    <h1>🧬 SeqGuide</h1>
    <div>Guided analysis of plasmid and amplicon sequencing data · version {__version__}</div>
    <p>Loading… The first visit takes about 30–60 seconds while your browser downloads the analysis engine;
       later visits are faster. Everything runs on <b>your</b> computer: your files are never uploaded.</p>
    <p style="font-size:14px;color:#666">Works in Chrome, Edge and Safari.</p>
    <div class="bar"><div></div></div>
  </div>
</div>
<script type="module">
import {{ mount }} from "https://cdn.jsdelivr.net/npm/@stlite/browser@{STLITE}/build/stlite.js";
const files = {json.dumps(files)};
mount({{
  requirements: ["numpy", "matplotlib", "jinja2"],
  entrypoint: "app/streamlit_app.py",
  files: files,
  streamlitConfig: {{
    "client.toolbarMode": "viewer",
    "server.maxUploadSize": 2000,
    "theme.base": "light",
  }},
}}, document.getElementById("root"));
</script>
</body>
</html>
"""

out = ROOT / "site"
out.mkdir(exist_ok=True)
(out / "index.html").write_text(html)
(out / ".nojekyll").write_text("")
print(f"Wrote {out / 'index.html'} ({len(html) / 1e3:.0f} kB, {len(files)} Python files, stlite {STLITE})")
