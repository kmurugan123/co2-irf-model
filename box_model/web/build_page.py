"""Build ../../index.html from page.html + model.js + page_data.json + mc_band.json + mc_future.json.  Run from this folder: python build_page.py"""
from pathlib import Path
here = Path(__file__).parent
s = (here / "page.html").read_text()
s = (s.replace("/*__MODEL__*/", (here / "model.js").read_text())
      .replace("/*__DATA__*/", (here / "page_data.json").read_text())
      .replace("/*__MC__*/", (here / "mc_band.json").read_text())
      .replace("/*__MCF__*/", (here / "mc_future.json").read_text()))
head = ('<!doctype html><html lang="en"><head><meta charset="utf-8">'
        '<meta name="viewport" content="width=device-width,initial-scale=1,viewport-fit=cover">'
        '<style>:root{color-scheme:light}body{margin:0}img{max-width:100%}[hidden]{display:none!important}</style></head><body>')
(here.parent.parent / "index.html").write_text(head + s + "</body></html>")
print("wrote index.html")
