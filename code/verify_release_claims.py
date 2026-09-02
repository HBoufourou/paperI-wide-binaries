#!/usr/bin/env python3
"""Fail if the committed v3 result files no longer support the stated claims."""

import json

import pandas as pd

from data_source import REPO_ROOT


audit = json.loads((REPO_ROOT / "results/catalog_provenance_audit.json").read_text())
dr3 = json.loads((REPO_ROOT / "results/corrected_dr3_summary.json").read_text())
phase_g = pd.read_csv(REPO_ROOT / "results/corrected_phase_G_results.csv")

assert audit["all_reduced_columns_match_mock_at_rtol_1e_9"] is True
assert audit["legacy_rows"] == 81_088
assert audit["corrected_real_rows"] == 81_880
assert dr3["catalogue"]["rows"] == 81_880
assert dr3["catalogue"]["md5"] == "1b6c5063163a4e6c07043d13aeb70f55"
assert dr3["status"] == "model-dependent / non-conclusive"
assert dr3["validation_control"]["empirical"]["passes_2_percent"] is False
assert dr3["validation_control"]["thermal"]["passes_2_percent"] is True
assert dr3["validation_control"]["superthermal"]["passes_2_percent"] is False

newton = phase_g[(phase_g.gamma == 1.0) & (phase_g.f_trip >= 0.2)]
injected = phase_g[phase_g.gamma == 1.4]
assert 1.061 < newton.E1.min() < 1.063
assert 1.127 < newton.E1.max() < 1.129
assert 1.001 < newton.E2.min() < 1.003
assert 1.043 < newton.E2.max() < 1.045
assert 1.305 < injected.E2.min() < 1.307
assert 1.345 < injected.E2.max() < 1.347

print("release claims verified against the machine-readable results")
