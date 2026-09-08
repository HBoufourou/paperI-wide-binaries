# Reproducing the final statistical branch

Ce dossier conserve les calculs conditionnels, leurs paramètres figés et les vérifications. Reproduire ces nombres ne prouve pas que les hypothèses décrivent entièrement Gaia. Les anciens résultats restent des travaux de développement, distincts de cette expérience prospective.

## Environment and inputs

Run commands from the package root, keeping the `work/final_paper` and `work/dr4_chain` layout. The statistical engine, summaries, tests and independent statistical verifier need Python and NumPy only. This execution used Python3.12.14 and NumPy2.3.5. Reproducing the exact random indices requires the same NumPy sampling implementation; matching library versions is recommended. Physical population/observation regeneration has additional dependencies documented by the root package, including pandas for inherited empirical covariates.

For reproducibility and moderate memory use, set `OPENBLAS_NUM_THREADS=1`, `OMP_NUM_THREADS=1`, `MKL_NUM_THREADS=1` and `NUMEXPR_NUM_THREADS=1` before starting Python. In PowerShell use `$env:OPENBLAS_NUM_THREADS='1'` and the analogous assignments.

The statistical calculation reads the44 population/observation records in `work/final_paper/bank_manifest.json`, their hashes, `work/final_paper/PROTOCOL.md`, and the frozen `stats/STATISTICAL_DECLARATION.json`. Paths inside the input manifest resolve against its parent folder. Some point to `../dr4_chain/` for the historical fitting and first calibration banks. No network access is required when these bank files are included. The model sources and declaration are hash-checked before every stage.

The independent bank proof is `work/final_paper/review/physical_banks_verification.json`. The original source/case declaration followed a13-test unit suite, the full software fixture, the independent fixture verifier and the recorded seed audit. `stats/SEED_AUDIT.json` and `stats/audit_history/` preserve that provenance; they are not a requirement to rerun the statistical arithmetic against the already frozen declaration. The audit itself additionally reads the historical `work/dr4_chain/decision/SEED_AUDIT.json`, `INPUT_INSPECTION.json`, and the completed software fixture proof.

## Verify the recorded run without regenerating it

These commands leave the recorded scientific profiles and decisions unchanged. The second writes a new independent report, whose path must be distinct from any verification proof one wishes to preserve.

```text
python -m unittest discover -s work/final_paper/stats -p test_engine.py -v
python work/final_paper/review/verify_statistics.py --root work/final_paper --out work/final_paper/review/statistics_verification_replay.json
```

The independent verifier imports no production engine. It reads all44 raw banks, the cached codes/templates, sampled indices, all calibration/test profiles and saved summary tables. It reconstructs measured categories, weights, mixtures, direct conditional deviance before profiling, integer rank probabilities, within-library combination, envelope, decisions, counts and paired tables. It does not regenerate the physical forces or observing process. A numerical PASS and physical/survey adequacy are different claims.

`python work/final_paper/stats/report_results.py --out work/final_paper/stats/results` rebuilds only descriptive CSV/Markdown/JSON summaries from hash-verified saved results. It does not alter scientific tables, profiles, p-values or thresholds. A later report-script edit changes its recorded source hash and should remain visible.

## Repeat all statistical sampling in a separate output folder

The following regenerates the same prescribed414 independent case streams, not a new scientific experiment. It reuses the recorded raw banks and frozen sources/declaration. Each output stage refuses to overwrite an existing stage, so choose an initially absent folder and retain old results. Budget approximately5GB of additional space for complete new profiles before any archive or duplicate banks; measure available storage first. Time depends on hardware and compression.

```text
python work/final_paper/stats/engine.py declare
python work/final_paper/stats/engine.py freeze --manifest work/final_paper/bank_manifest.json --out work/final_paper/stats/replay_results
python work/final_paper/stats/engine.py calibrate --manifest work/final_paper/bank_manifest.json --out work/final_paper/stats/replay_results
python work/final_paper/stats/engine.py test --manifest work/final_paper/bank_manifest.json --out work/final_paper/stats/replay_results
python work/final_paper/stats/report_results.py --out work/final_paper/stats/replay_results
python work/final_paper/review/verify_statistics.py --root work/final_paper --results work/final_paper/stats/replay_results --out work/final_paper/review/statistics_verification_resampling.json
```

`declare` accepts the existing declaration only when its content exactly matches current sources and protocol. It does not replace a changed declaration. Do not copy a newly generated manifest hash into an old recorded execution. If regenerating physical banks changes container hashes or execution metadata, perform a new input freeze into another directory, retain both manifests, and label it a reproduction run. A source/protocol change instead requires a separately declared experiment, not this replay command.

The small software fixture has44 wholly artificial banks and repeats2. Its original proof and complete copied source are under `stats/software_fixture`. Its very coarse rank resolution intentionally gives only indeterminate decisions and supplies no physical evidence. `software_smoke.py` refuses to replace that directory; run it in a disposable copied source layout if a fresh full fixture is needed. The unit command above is normally enough to repeat the light checks without duplicating inputs.

## Interpretation and accounting

Calibration contains234 cases, three paired cadences and1000 repeats:702000 catalogue views. Test contains180 cases, three paired cadences and1000 repeats:540000 views,540 distinct physical/cadence cells. Each view feeds all five information variants and both calibration rules, producing5400 metric rows. These are paired representations, not independent new physical samples. `paired.csv` contains7020 contrasts, and `envelope_effects.csv` preserves2700 envelope comparisons including true-model retention and adverse transitions.

The39-point nuisance grid contains beta1/1.3/1.6, f0/.075/.15/.225/.3 and mixtures of inner40/80 priors. f=.15 and eta=.5 are calibration/profiling points without dedicated independent test cells. Test priors60/160 and beta1.15/1.45 exercise specified transfer only. No continuous nuisance domain is certified by that finite list.

One index array per case is shared by modes and cadences. Count arrays are reconstructed and canonically hashed instead of being saved repeatedly. Full nuisance surfaces remain float64. `OUTPUT_CONTRACT.md` specifies all axes, formulas, count hashing, output labels and the measured-input replay API.

Correct/wrong count criteria (95%/1%) are separate from true-model-retention99%. Wilson intervals are marginal cellwise intervals and do not become simultaneous guarantees when worst cells are selected. The calibration-library envelope covers the listed empirical laws; transfer to other finite libraries or real Gaia must be assessed empirically. The conditional score removes a marginal likelihood term but its calibration is not conditioned on exactly fixed observed stratum counts. All reported conclusions remain conditional on the physical parent, companion prescription, covariates and observing model; `sky_ready` stays false.
