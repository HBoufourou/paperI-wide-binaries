#!/usr/bin/env python3
"""Build the corrected preprint PDF from ``paperI.tex`` using latexmk."""

from __future__ import annotations

import shutil
import subprocess

from data_source import REPO_ROOT


SOURCE = REPO_ROOT / "paperI.tex"
OUTPUT = REPO_ROOT / "Boufourou_2026_PaperI_corrected_v3.pdf"


def main():
    if shutil.which("latexmk") is None:
        raise SystemExit("latexmk is required to build the manuscript PDF")
    subprocess.run(
        [
            "latexmk",
            "-pdf",
            "-g",
            "-interaction=nonstopmode",
            "-halt-on-error",
            SOURCE.name,
        ],
        cwd=REPO_ROOT,
        check=True,
    )
    built = SOURCE.with_suffix(".pdf")
    built.replace(OUTPUT)
    print(f"wrote {OUTPUT}")


if __name__ == "__main__":
    main()
