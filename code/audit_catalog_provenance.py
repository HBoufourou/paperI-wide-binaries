#!/usr/bin/env python3
"""Reproduce the catalogue-provenance finding that triggered the correction."""

from __future__ import annotations

import json
from pathlib import Path
from urllib.request import urlopen

import numpy as np
import pandas as pd

from data_source import DEFAULT_DATA_DIR, REPO_ROOT, load_catalog, md5sum

MOCK_NAME = "Newton_dr3_MSMS_d200pc_5.csv"
MOCK_URL = f"https://zenodo.org/api/records/10652994/files/{MOCK_NAME}/content"
MOCK_MD5 = "997069b5635200853896d005a238dde9"
LEGACY = REPO_ROOT / "data" / "legacy_newton_mock_v5_reduced.csv"


def download(url: str, destination: Path, checksum: str) -> Path:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists() and md5sum(destination) == checksum:
        return destination
    temporary = destination.with_suffix(destination.suffix + ".part")
    with urlopen(url) as response, temporary.open("wb") as output:
        while block := response.read(1024 * 1024):
            output.write(block)
    actual = md5sum(temporary)
    if actual != checksum:
        temporary.unlink(missing_ok=True)
        raise ValueError(f"Checksum mismatch: {actual} != {checksum}")
    temporary.replace(destination)
    return destination


def close_fraction(left, right):
    return float(np.isclose(left, right, rtol=1e-9, atol=0, equal_nan=True).mean())


def main():
    mock_path = download(MOCK_URL, DEFAULT_DATA_DIR / MOCK_NAME, MOCK_MD5)
    real = load_catalog()
    legacy = pd.read_csv(LEGACY, dtype={"source_id1": "string", "source_id2": "string"})
    mock = pd.read_csv(mock_path, dtype={"source_id1": "string", "source_id2": "string"})

    keys = ["source_id1", "source_id2"]
    mock_match = legacy.merge(mock, on=keys, suffixes=("_legacy", "_mock"), validate="1:1")
    real_match = legacy.merge(real, on=keys, suffixes=("_legacy", "_real"), validate="1:1")
    column_match = {}
    for column in legacy.columns:
        if column in keys:
            continue
        column_match[column] = close_fraction(
            mock_match[f"{column}_legacy"].to_numpy(),
            mock_match[f"{column}_mock"].to_numpy(),
        )

    w1 = 1 / real_match["d1_err[pc]_real"] ** 2
    w2 = 1 / real_match["d2_err[pc]_real"] ** 2
    mean_distance = (
        real_match["d1[pc]_real"] * w1 + real_match["d2[pc]_real"] * w2
    ) / (w1 + w2)
    allowed_offset = real_match["s[kau]_legacy"] * 1000 / 206265

    relative_mock = np.hypot(
        real_match["mu1ra[mas/yr]_legacy"] - real_match["mu2ra[mas/yr]_legacy"],
        real_match["mu1dec[mas/yr]_legacy"] - real_match["mu2dec[mas/yr]_legacy"],
    )
    relative_real = np.hypot(
        real_match["mu1ra[mas/yr]_real"] - real_match["mu2ra[mas/yr]_real"],
        real_match["mu1dec[mas/yr]_real"] - real_match["mu2dec[mas/yr]_real"],
    )

    report = {
        "finding": "legacy observational input is Chae's virtual Newtonian v5 catalogue",
        "legacy_rows": int(len(legacy)),
        "official_mock_rows": int(len(mock)),
        "corrected_real_rows": int(len(real)),
        "matched_to_mock_by_source_id": int(len(mock_match)),
        "matched_to_real_by_source_id": int(len(real_match)),
        "all_reduced_columns_match_mock_at_rtol_1e_9": bool(
            all(value == 1.0 for value in column_match.values())
        ),
        "column_match_fraction_to_mock": column_match,
        "distance_fingerprint": {
            "legacy_median_abs_d1_minus_real_weighted_mean_pc": float(
                np.median(np.abs(real_match["d1[pc]_legacy"] - mean_distance))
            ),
            "real_median_abs_d1_minus_real_weighted_mean_pc": float(
                np.median(np.abs(real_match["d1[pc]_real"] - mean_distance))
            ),
            "legacy_fraction_within_projected_orbital_offset": float(
                np.mean(np.abs(real_match["d1[pc]_legacy"] - mean_distance) <= allowed_offset * 1.0001)
            ),
            "real_fraction_within_projected_orbital_offset": float(
                np.mean(np.abs(real_match["d1[pc]_real"] - mean_distance) <= allowed_offset * 1.0001)
            ),
        },
        "relative_proper_motion_correlation_with_real": float(
            np.corrcoef(relative_mock, relative_real)[0, 1]
        ),
        "sources": {
            "mock_zenodo_record": 10652994,
            "mock_filename": MOCK_NAME,
            "mock_md5": MOCK_MD5,
            "real_zenodo_record": 10986733,
        },
    }
    results = REPO_ROOT / "results"
    results.mkdir(exist_ok=True)
    destination = results / "catalog_provenance_audit.json"
    destination.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    print(f"\nwrote {destination}")


if __name__ == "__main__":
    main()
