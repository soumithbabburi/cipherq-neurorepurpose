"""Build a styled HTML + PDF from VALIDATION_NARRATIVE.md.

  python build_validation_doc.py

Outputs (repo root):
  VALIDATION_NARRATIVE.html  — styled, print-ready (open in browser, Ctrl+P → Save as PDF)
  VALIDATION_NARRATIVE.pdf   — generated via xhtml2pdf
"""
import re
from pathlib import Path

import markdown
from xhtml2pdf import pisa

ROOT = Path(__file__).resolve().parent
MD = ROOT / "VALIDATION_NARRATIVE.md"
HTML = ROOT / "VALIDATION_NARRATIVE.html"
PDF = ROOT / "VALIDATION_NARRATIVE.pdf"

CSS = """
@page { size: A4; margin: 1.8cm 1.7cm; }
body { font-family: Helvetica, Arial, sans-serif; font-size: 10.5pt; color: #1f2a26; line-height: 1.5; }
h1 { font-size: 20pt; color: #2f6f5b; margin: 0 0 2pt; }
h2 { font-size: 14pt; color: #2f6f5b; border-bottom: 1.5px solid #d9e2dd; padding-bottom: 3px; margin: 18px 0 8px; }
h3 { font-size: 11.5pt; color: #2a3b34; margin: 12px 0 4px; }
p, li { font-size: 10.5pt; }
strong { color: #1a241f; }
code { font-family: Courier, monospace; font-size: 9.5pt; background: #f3f6f4; padding: 1px 3px; }
pre { background: #f3f6f4; border: 1px solid #e0e8e3; padding: 8px 10px; font-size: 9pt; }
table { width: 100%; border-collapse: collapse; margin: 8px 0 14px; }
th { background: #eef3f0; color: #2a3b34; text-align: left; font-size: 9.5pt; padding: 5px 7px; border: 1px solid #d2ddd7; }
td { font-size: 9.5pt; padding: 5px 7px; border: 1px solid #e2e9e5; vertical-align: top; }
hr { border: none; border-top: 1px solid #d9e2dd; margin: 16px 0; }
a { color: #2f6f5b; text-decoration: none; }
.cover { border-left: 4px solid #2f6f5b; padding: 4px 0 4px 14px; margin: 0 0 18px; }
.cover .sub { color: #6b7a73; font-size: 10pt; }
"""

HTML_TMPL = """<!DOCTYPE html><html><head><meta charset="utf-8">
<title>RepurposeIQ — Data Validation</title><style>{css}</style></head>
<body>{body}</body></html>"""

# unicode → PDF-safe (xhtml2pdf's core fonts lack some glyphs)
PDF_SAFE = {
    "→": "-&gt;", "←": "&lt;-", "≤": "&le;", "≥": "&ge;",
    "×": "x", "•": "&bull;", "–": "-", "—": "&mdash;",
    "‘": "'", "’": "'", "“": '"', "”": '"', "…": "...",
}


def render(md_text, for_pdf=False):
    body = markdown.markdown(md_text, extensions=["tables", "fenced_code", "sane_lists"])
    if for_pdf:
        for u, rep in PDF_SAFE.items():
            body = body.replace(u, rep)
    return HTML_TMPL.format(css=CSS, body=body)


def main():
    md_text = MD.read_text(encoding="utf-8")
    # styled HTML (keeps full unicode — best for browser Print-to-PDF)
    HTML.write_text(render(md_text, for_pdf=False), encoding="utf-8")
    print(f"wrote {HTML.name}")
    # PDF via xhtml2pdf
    with open(PDF, "wb") as f:
        result = pisa.CreatePDF(render(md_text, for_pdf=True), dest=f, encoding="utf-8")
    print("PDF errors:" , result.err)
    print(f"wrote {PDF.name}  ({PDF.stat().st_size//1024} KB)")


if __name__ == "__main__":
    main()
