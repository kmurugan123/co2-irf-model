"""Build ../../index.html from page.html and its parts.  Run from this folder: python build_page.py"""
from pathlib import Path
here = Path(__file__).parent
s = (here / "page.html").read_text()
for key, name in [("MODEL", "model.js"), ("MCJS", "mc.js"), ("DATA", "page_data.json"), ("MC", "mc_band.json"),
                  ("MCF", "mc_future.json"), ("CAL", "calib.json")]:
    s = s.replace(f"/*__{key}__*/", (here / name).read_text())
head = ('<!doctype html><html lang="en"><head><meta charset="utf-8">'
        '<meta name="viewport" content="width=device-width,initial-scale=1,viewport-fit=cover">'
        '<style>:root{color-scheme:light}body{margin:0}img{max-width:100%}[hidden]{display:none!important}</style></head><body>')
(here.parent.parent / "index.html").write_text(head + s + "</body></html>")
print("wrote index.html")
