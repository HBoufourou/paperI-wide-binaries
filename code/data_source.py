#!/usr/bin/env python3
"""Authoritative data loader for the corrected DR3 analysis.

The v1/v2 repository accidentally analysed Chae's virtual Newtonian catalogue
(`Newton_dr3_MSMS_d200pc_5.csv`).  The corrected analysis downloads the public
Gaia catalogue with repaired RUWE values from Zenodo record 10986733, verifies
its checksum, and refuses files that look like the 81,088-row mock catalogue.
"""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path
from urllib.request import urlopen

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DATA_DIR = REPO_ROOT / "data" / "external"
CATALOG_NAME = "gaia_dr3_MSMS_d200pc_ruwe.csv"
CATALOG_URL = (
    "https://zenodo.org/api/records/10986733/files/"
    "gaia_dr3_MSMS_d200pc_ruwe.csv/content"
)
CATALOG_MD5 = "1b6c5063163a4e6c07043d13aeb70f55"
CATALOG_ROWS = 81_880
MOCK_ROWS = 81_088

REQUIRED_COLUMNS = {
    "source_id1",
    "source_id2",
    "R_chance",
    "s[kau]",
    "d1[pc]",
    "d1_err[pc]",
    "d2[pc]",
    "d2_err[pc]",
    "M1[Msun]",
    "M2[Msun]",
    "mu1ra[mas/yr]",
    "mu1ra_err[mas/yr]",
    "mu1dec[mas/yr]",
    "mu1dec_err[mas/yr]",
    "mu2ra[mas/yr]",
    "mu2ra_err[mas/yr]",
    "mu2dec[mas/yr]",
    "mu2dec_err[mas/yr]",
    "RV1[km/s]",
    "RV1_err[km/s]",
    "RV2[km/s]",
    "RV2_err[km/s]",
    "ruwe1",
    "ruwe2",
    "e",
    "MagG1",
    "MagG2",
    "bp_rp1",
    "bp_rp2",
    "RA1[deg]",
    "DEC1[deg]",
    "RA2[deg]",
    "DEC2[deg]",
    "A_G1[mag]",
    "A_G2[mag]",
}


def md5sum(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def download_catalog(destination: Path | None = None) -> Path:
    destination = Path(destination or (DEFAULT_DATA_DIR / CATALOG_NAME))
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists() and md5sum(destination) == CATALOG_MD5:
        return destination

    temporary = destination.with_suffix(destination.suffix + ".part")
    with urlopen(CATALOG_URL) as response, temporary.open("wb") as output:
        while block := response.read(1024 * 1024):
            output.write(block)
    checksum = md5sum(temporary)
    if checksum != CATALOG_MD5:
        temporary.unlink(missing_ok=True)
        raise ValueError(
            f"Checksum mismatch for {CATALOG_NAME}: {checksum} != {CATALOG_MD5}"
        )
    temporary.replace(destination)
    return destination


def load_catalog(path: str | Path | None = None, *, download: bool = True) -> pd.DataFrame:
    catalog_path = Path(path) if path else DEFAULT_DATA_DIR / CATALOG_NAME
    if catalog_path.name.startswith("Newton_"):
        raise ValueError(
            "Refusing a Newton_* virtual catalogue for the observational DR3 analysis."
        )
    if not catalog_path.exists():
        if path or not download:
            raise FileNotFoundError(catalog_path)
        catalog_path = download_catalog(catalog_path)

    checksum = md5sum(catalog_path)
    if checksum != CATALOG_MD5:
        raise ValueError(
            f"Checksum mismatch for the observational catalogue: "
            f"{checksum} != {CATALOG_MD5}"
        )

    frame = pd.read_csv(
        catalog_path,
        dtype={"source_id1": "string", "source_id2": "string"},
    )
    missing = sorted(REQUIRED_COLUMNS.difference(frame.columns))
    if missing:
        raise ValueError(f"Catalogue is missing required columns: {missing}")
    if len(frame) == MOCK_ROWS:
        raise ValueError(
            "Refusing an 81,088-row catalogue: this is the size of Chae's virtual "
            "Newtonian sample, not the 81,880-row Gaia catalogue."
        )
    if len(frame) != CATALOG_ROWS:
        raise ValueError(f"Unexpected catalogue length: {len(frame)} != {CATALOG_ROWS}")
    return frame


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--path", type=Path, help="Optional existing catalogue path")
    parser.add_argument(
        "--no-download", action="store_true", help="Fail rather than download missing data"
    )
    args = parser.parse_args()
    frame = load_catalog(args.path, download=not args.no_download)
    selected = args.path or (DEFAULT_DATA_DIR / CATALOG_NAME)
    print(f"verified catalogue: {selected}")
    print(f"rows: {len(frame):,}; columns: {len(frame.columns)}")
    if Path(selected).name == CATALOG_NAME:
        print(f"md5: {md5sum(Path(selected))}")


if __name__ == "__main__":
    main()
