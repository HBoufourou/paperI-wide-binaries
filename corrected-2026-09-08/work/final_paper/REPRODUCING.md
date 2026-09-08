# Reproducing the conditional experiment

The archive root contains `work/`. Run commands from that root, preserving the relative directory layout. The released files are the reference record. Perform new writes and regeneration in a separate extraction or working copy; retain the original archive and its checksums.

This is a synthetic, conditional Newton/QUMOND experiment. Reproducing its calculations does not validate its population assumptions on Gaia, establish a general error guarantee, or rebuild every possible fitting library. The actual physical and statistical scopes are in `work/final_paper/PROTOCOL.md` and the manuscript.

## 1. Environment and small checks

The numerical execution used **Python3.12.14, NumPy2.3.5 and pandas3.0.1**, on Windows64-bit. The inference engine itself needs Python and NumPy. Pandas is an import dependency of the inherited empirical-covariate provider used for physical regeneration. Scientific execution does not require the PDF-rendering environment, a GPU, SciPy or an internet download. Preserve the included empirical arrays, model, validation records and historical inputs.

The examples below use PowerShell. Set `$Python` to the matching interpreter, either available on PATH or an explicit executable path. Limit numerical libraries before launching Python; run heavy stages sequentially.

```powershell
$Python = 'python'
$env:OMP_NUM_THREADS = '1'
$env:OPENBLAS_NUM_THREADS = '1'
$env:MKL_NUM_THREADS = '1'
& $Python -c "import sys,numpy,pandas; print(sys.version); print(numpy.__version__,pandas.__version__)"
& $Python -m unittest discover -s work/final_paper/stats -p test_engine.py -v
& $Python work/final_paper/stats/engine.py declare
```

The unit suite contains13 tests: measured-field restrictions, units and source relabelling, invalid covariance, histogram projection, conditional-mixture order, envelope order, tail ties, streams, shared indices, classification versus retention, and the catalogue API. The command above was checked against the released source. `declare` verifies the existing declaration against the current code, protocol and complete case schedule; it refuses a changed declaration. In an archive with the declaration already present, it does not create a new scientific run.

`stats/SEED_AUDIT.json` and `stats/audit_history/` retain the seed review. In a working copy, the complete audit can be repeated with:

```powershell
& $Python work/final_paper/stats/audit_and_declare.py
```

This checks414 scientific catalogue streams,64 new physical/observer streams, isolated fixture streams and the inherited audit. It requires the included software-fixture proof and historical audit records. It rewrites the current audit JSON and preserves a content-addressed audit history; it does not draw physical populations or scientific catalogues. A different path, runtime or failed check must not be concealed by editing a frozen scientific declaration.

## 2. Inspect the saved experiment

Start with `work/final_paper/stats/results/REPORT.md`. The authoritative cell-level tables are:

| File in `stats/results/` | Contents |
| --- | --- |
| `metrics.csv` | All5400 method/rule/cell rows: counts, correct/wrong/indeterminate rates, true-model retention, Wilson intervals, support and nuisance summaries. |
| `paired.csv` | All7020 paired comparisons; every difference is first minus second, with discordant counts. |
| `method_summary.csv`, `family_summary.csv`, `realization_summary.csv` | Descriptive summaries; the full cell tables remain authoritative. |
| `count_failures.csv`, `confidence_failures.csv`, `worst_cells.csv` | Failed targets and exact parameters behind the extremes. |
| `envelope_effects.csv`, `envelope_summary.csv` | Paired changes under envelope versus baseline calibration. |
| `information_summary.csv` | Signed information-block comparisons, with positive/zero/negative changes and discordant totals. |
| `library_support.csv`, `fitting_component_support.csv` | Weight concentration and component support; these are not guarantees of physical support in every histogram bin. |
| `report_artifacts.json` | Summary provenance and recorded execution times. |

The540 physical-case/scenario cells share their catalogue indices across five variants and two rules. They do not represent5400 independent physical experiments. Each cell contains1000 catalogues of1000 parents. Inspect both test realizations and both causes of indeterminacy; do not replace a failed cell with a pooled mean. Wilson intervals are marginal and conditional on the realized libraries.

`stats/OUTPUT_CONTRACT.md` defines all array axes, integer codes, mixture order, canonical count hashes and the interpretation of every output. Each calibration/test `caseNNN.npz` saves a single parent-index array, full nuisance surfaces and profiled scores. Test files also save probabilities, decisions and reasons. Count arrays are reconstructed from saved indices and measured codes rather than duplicated on disk. The calibration reference contains both statistics for all39 nuisance points, two gravities and three calibration realizations.

## 3. Repeat the verification

The released reports have distinct scopes:

- `review/physical_banks_verification.json`:44 complete banks,880000 accepted initial rows and132 observation files;6143 checks passed in the released run. It independently reconstructs initial invariants, masses, selection and weights and checks recorded observation identities. It does not replay every rejected proposal or regenerate every trajectory/noise/GLS fit.
- `review/statistics_verification.json`: independently reconstructs all measured codes, templates, indices/counts, score surfaces, rank values, decisions and descriptive tables, without importing the production inference engine. Inspect its actual `complete`, `passed`, scope and failure fields.
- `review/delivery_api_verification.json`: bounded replay of the catalogue interface using four saved catalogues from test case3, three scenarios, five modes and two rules. Its120 calls are not a replacement for the complete independent statistical reconstruction.

In a working copy, these commands reproduce the checks. The physical checker writes its two reports beside itself; the other output paths are explicitly selectable.

```powershell
& $Python work/final_paper/review/verify_fresh_banks.py
& $Python work/final_paper/review/verify_statistics.py --out work/final_paper/review/statistics_recheck.json
& $Python work/final_paper/verify_delivery_api.py --out work/final_paper/review/delivery_api_recheck.json
```

The full statistical check is computational work, not merely a checksum scan. The independent checker requires the final test manifest and will not report a completed scientific pass from an incomplete run. Its `--wait` option can wait for that manifest. No physics is regenerated by this statistical checker. The API replay uses the released `stats/results` path; its `--out` changes only its verification report.

## 4. Recreate statistics using the unchanged physical inputs

Use a separate working copy with the released model and all44 referenced banks intact. An alternative statistical destination is supported. The following example must start with a **new** `stats/reproduced` directory:

```powershell
& $Python work/final_paper/stats/engine.py declare
& $Python work/final_paper/stats/engine.py freeze --manifest work/final_paper/bank_manifest.json --out work/final_paper/stats/reproduced
& $Python work/final_paper/stats/engine.py calibrate --manifest work/final_paper/bank_manifest.json --out work/final_paper/stats/reproduced
& $Python work/final_paper/stats/engine.py test --manifest work/final_paper/bank_manifest.json --out work/final_paper/stats/reproduced
& $Python work/final_paper/review/verify_statistics.py --manifest work/final_paper/bank_manifest.json --results work/final_paper/stats/reproduced --out work/final_paper/review/statistics_reproduced.json
& $Python work/final_paper/stats/report_results.py --out work/final_paper/stats/reproduced
```

`freeze` hashes the inputs, reconstructs the codes/templates and writes a new execution record. `calibrate` uses the frozen234 cases and `test` uses180 cases; both use the declared seeds and1000 repetitions. No alternative binning, nuisance grid or threshold is introduced. This reproduces the same experiment; it is not an additional independent power test.

The source/protocol declaration remains unchanged. The new execution JSON records its actual manifest location, so some provenance hashes legitimately differ after relocation even when numerical arrays agree. Runtime and elapsed-time metadata also differ. Do not compare entire JSON or archive bytes as a substitute for comparing the scientific arrays and decisions. Preserve fresh hashes and investigate any failed numerical check; binary-identical output is not promised across different NumPy, linear-algebra or compression builds.

**Restart behavior:** completed frozen inputs can be followed by the next unstarted stage. A stage directory that already exists is refused: there is no mid-calibration or mid-test case-resume feature. After interruption, preserve that incomplete output directory and begin again in another new `--out` destination. Do not delete partial cases and pretend the stage resumed. `--fixture-repeats` is allowed only with explicit software-fixture mode and is not a scientific acceleration option.

## 5. Regenerate the32 fresh physical banks

This optional heavier reproduction retains the potential, empirical reservoir and12 historical fitting/calibration inputs. It does **not** rebuild the six original fitting populations or solve the potential grid again. The producer imports byte-checked upstream sources. Keep the entire included dependency layout, the upstream manifests and the scalar-model validation files intact.

The producer has no `--out`: it writes relative to its own `work/final_paper/` directory. Use a second clean extraction of the archive. In that copy only, preserve the released `work/final_paper/banks/` directory under another name and preserve `work/final_paper/bank_manifest.json` under another name. Leave `GENERATION_PROTOCOL.json`, all scientific sources and all upstream inputs unchanged. The original archive or first extraction remains the authoritative release; renaming the local copies does not retarget their historical manifest paths.

With `banks/` and `bank_manifest.json` absent **only in that reproduction copy**, run:

```powershell
& $Python work/final_paper/generate_fresh_banks.py --freeze-only
& $Python work/final_paper/generate_fresh_banks.py
& $Python work/final_paper/review/verify_fresh_banks.py
```

The first command validates the frozen generation specification and dependencies. The second generates32 fresh banks using the original64 streams,20,000 accepted parents each and the declared effective-size rule permitting40,000. It creates a manifest containing those32 plus12 fixed historical references and authorizes it only on completion. The third recomputes the physical-bank checks. Read its completed result before running the statistical commands in section4 against the newly generated manifest.

Every completed bank has its population, three observed views, hashes and metadata. Rerunning the producer checks and skips completed banks already recorded in the top-level manifest. If it finds NPZ files in an unrecorded/incomplete bank directory, it stops for inspection; it does not silently reuse partial output. Preserve that reproduction attempt and start a clean one, or inspect it explicitly without rewriting the released evidence.

The elapsed times in regenerated bank metadata change the top-level manifest hash. Regenerated statistical outputs must therefore receive their own execution record; the old statistical execution cannot simply be pointed at a different manifest. Reusing the same frozen seeds tests reproducibility, not a third independently seeded calibration/test-library realization.

## 6. Measured-catalogue API

The released `stats/engine.py` includes `predict_catalogue`; its source SHA-256 is `87bc9208e25acc6b97e7a64f444c3e47f43e5cb8fe3dc2aac8d41d037ca8a363`. This API was added before the final statistical source freeze. It uses the same scores, templates and calibration as the saved experiment.

Supply the matching verified template dictionary and `reference.npz["sorted"]`, and1000 rows of measured inputs. Parent fields are `distance_pc` in pc, `Mphot1` and `Mphot2` in solar masses, and `rproj_kau` with shape(N,2). `features.OBS_FIELDS[mode]` gives the exact observation whitelist, including proper motions in mas/yr, covariance in(mas/yr)², and the defined diagnostic scores. The reduced modes need only their own listed fields. Latent gravity, energy, angular momentum, true mass and companion labels are not inputs to this prediction.

For caller-supplied measured dictionaries and already verified matching artifacts, the invocation is:

```python
answer = engine.predict_catalogue(
    measured_parent, measured_observations, "uniform_white",
    template, reference, mode="full66", rule="envelope",
    declared_scope_confirmed=True,
)
```

Before using artifacts, check the declaration and the execution/template/reference hashes as in `verify_delivery_api.py`; the low-level prediction function expects the caller to provide matching validated artifacts. Its scope flag is an explicit caller assertion, not a numerical test proving that an external population satisfies the simulator.

A supported response has `status="conditional"`, with `decision` equal to `Newton`, `QUMOND` or `indeterminate`; reasons distinguish `one_retained`, `both_retained` and `both_rejected`. Probabilities follow hypothesis order Newton,QUMOND. Detected unsupported input returns `status="unsupported"`, `decision="indeterminate"` and the reason. Every response keeps `sky_ready=false`. Neither successful invocation nor an exclusive simulated decision establishes that a real Gaia catalogue is calibrated for this method.
