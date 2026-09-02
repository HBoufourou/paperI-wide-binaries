#!/usr/bin/env python3
"""Print the corrected, model-dependent DR3 verdict from the audited run."""

import json

from data_source import REPO_ROOT


path = REPO_ROOT / "results" / "corrected_dr3_summary.json"
if not path.exists():
    raise SystemExit("Run code/corrected_dr3_analysis.py first.")

summary = json.loads(path.read_text())
print("CORRECTED DR3 STATUS:", summary["status"])
print("validation control:", json.dumps(summary["validation_control"], indent=2))
print("E2 test:", json.dumps(summary["e2"]["test"], indent=2))
print("E2 deep:", json.dumps(summary["e2"]["deep"], indent=2))
print("E3:", json.dumps(summary["e3"], indent=2))
print("\nVerdict: model-dependent / non-conclusive; v1/v2 observational claims withdrawn.")
