#!/usr/bin/env python3
"""
Zip the files arXiv needs to compile the manuscript.

The .tex, .bib, figures, and class/bst files sit at the archive root
(not under paper/). The compiled PDF, .bbl, and aux files are omitted; arXiv
regenerates the .bbl from the .bib with bibtex.

    python scripts/arxiv_zip.py              # write engine-lrds-arxiv.zip
    python scripts/arxiv_zip.py --out PATH   # choose the zip path
    python scripts/arxiv_zip.py --dry-run    # list files, do not write
"""

import argparse
import sys
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PAPER = ROOT / "paper"
DEFAULT_ZIP = ROOT / "engine-lrds-arxiv.zip"

# Names as they appear in paper/ and at the zip root.
FILES = [
    "lrd_engine.tex",
    "lrd.bib",
    "aasjournal.bst",
    "twindow.pdf",
    "fplane.pdf",
]


def collect():
    missing = [name for name in FILES if not (PAPER / name).is_file()]
    if missing:
        print("Missing from paper/:", file=sys.stderr)
        for name in missing:
            print(f"  {name}", file=sys.stderr)
        sys.exit(1)
    return [PAPER / name for name in FILES]


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--out", type=Path, default=DEFAULT_ZIP, help="zip path")
    ap.add_argument("--dry-run", action="store_true", help="list files only")
    args = ap.parse_args()

    files = collect()
    if args.dry_run:
        for f in files:
            print(f.name)
        return
    with zipfile.ZipFile(args.out, "w", zipfile.ZIP_DEFLATED) as zf:
        for f in files:
            zf.write(f, f.name)
    print(f"wrote {args.out} ({len(files)} files)")


if __name__ == "__main__":
    main()
