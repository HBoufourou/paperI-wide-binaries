# Coupled-library decision engine

L’expérience conditionnelle est exécutée : 22 bibliothèques nouvelles et vérifiées, 156 000 vues de calibration et 216 000 vues de test. Les 216 cellules atteignent les cibles observées ≥95 % de décisions correctes et ≤1 % de décisions fausses ; 208 atteignent aussi les bornes Wilson plus strictes. Les huit exceptions concernent Newton avec beta=1,3, hors grille ajustée. Ce résultat sur des populations simulées ne constitue pas un verdict sur Gaia ni une garantie pour DR4. La reconstruction statistique indépendante a passé 11 558 contrôles ; sa preuve et sa portée sont présentées séparément dans ../review/.

The parent protocol and exact bank schema are in PROTOCOL.md and BANK_CONTRACT.md. Runtime dependencies are Python3.12 and NumPy2.3.5; CSV/JSON output uses the standard library. No orbit generator, optimizer package, internet, plotting or prior experiment result is imported by this engine.

```powershell
$env:OPENBLAS_NUM_THREADS='1'
$env:OMP_NUM_THREADS='1'
python -m unittest discover -s work/dr4_chain/decision -p test_engine.py -v
```

Eight executed tests cover units and signed projections; source relabelling and ignoring hidden fields; invalid shapes/covariance; category and weighted-probability construction; profiled scores against a scalar deviance; p-values against direct counting with ties; distinct seeds/count-versus-confidence criteria; and weighted index pairing across cadences.

software_smoke.py additionally exercises the whole manifest→templates→calibration→test path on22 artificial48-row banks. It uses eight calibration/catalogue repeats, so every output must be indeterminate at the fixed1% rejection threshold. It is not an orbit simulation, independent scientific validation or demonstration of physical separation. Its corrected execution is in software_fixture/, with software_only flags throughout. A first plumbing run accidentally inherited unused future scientific seed numbers; its complete archive and source are retained under software_fixture_initial_seed_audit/. SOFTWARE_SEED_AUDIT.md records the pre-physical correction. All corrected fixture sampling uses279M streams and future scientific sampling uses270/271M calibration,274M Newton tests and277M QUMOND tests. The subsequent seed-only amendment is recorded in AMENDMENT_01_TEST_SEEDS.md; byte-identical sources for the corrected executed fixture are preserved in software_fixture/executed_source/, and its existing verification is unchanged.

The software fixture script refuses to overwrite its existing output. Re-execute it only in a fresh source copy if desired; do not remove the archived fixture history to rerun. Tests may be rerun without generating physical data. A report of the supplied fixture can be rendered without scoring again:

```powershell
python work/dr4_chain/decision/report_results.py --out work/dr4_chain/decision/software_fixture/outputs
```

The physical domain, trajectories and all 22 banks passed their declared and independent checks before the following scientific sequence was executed. These commands now refuse the occupied execution directory; reproduce in a separate clean source/output copy or inspect the supplied saved results:

```powershell
python work/dr4_chain/decision/engine.py freeze --manifest work/dr4_chain/bank_manifest.json --out work/dr4_chain/decision/results
python work/dr4_chain/decision/engine.py calibrate --manifest work/dr4_chain/bank_manifest.json --out work/dr4_chain/decision/results
python work/dr4_chain/decision/engine.py test --manifest work/dr4_chain/bank_manifest.json --out work/dr4_chain/decision/results
python work/dr4_chain/decision/report_results.py --out work/dr4_chain/decision/results
```

Scientific defaults are fixed:1000 parents/catalogue,1000 catalogues/cell,26 discrete nuisance points/hypothesis,156000 calibration views and216000 test views. Repeats may be reduced only with an explicit software-fixture purpose at freeze; the override cannot be used for production. Scenarios intentionally share sampled parent index arrays. Files retain complete counts, profile surfaces, best and near-best nuisance indices, p-values, decisions, reasons, hashes and sampling lineage. Stages refuse overwriting frozen output. Independent reconstruction of physical and statistical outputs belongs to the separate ../review/ task; its actual status must be read from its verification files.

For an already-loaded model, engine.infer_observed accepts measured population/observation dictionaries, one declared scenario, templates and sorted calibration reference arrays. It requires exactly1000 parents and explicit confirmation that the declared scope applies; otherwise it returns indeterminate/unsupported. The physical provenance assertion is not a proof that the actual sky follows the model. The engine cannot certify unmodelled selection, continuous nuisance interpolation, photometric calibration or a different observing process.

The saved scientific results are results/REPORT.md, metrics.csv (all216 cells), paired_cadences.csv (144 comparisons), group_summary.csv (18 groups), worst_cells.csv (54 descriptive extrema), full test files and calibration references. All192 beta-endpoint cells, including actual60/160 AU prior tests, meet the confidence targets. The outside-beta-grid Newton case reaches a worst observed950/1000 correct decisions with Wilson interval[93.469%,96.187%]; the maximum wrong rate is10/1000 with interval[0.544%,1.831%]. These are different cells and the latter cannot establish a≤1% underlying wrong rate.

Every indeterminate result in this executed test is a double rejection. Within the fitted family, minimum true-model retention is98.3% with marginal Wilson interval[97.294%,98.936%]. Thus good classification targets do not establish nominal99% model retention after transfer between independently generated finite libraries. The reported intervals are marginal; the repeated catalogues reuse finite banks.

In test records eta=0 is an ignored interface placeholder when prior is supplied. Prior40 corresponds to fitted eta=0, prior80 to eta=1, and actual60/160 priors have no true mixture eta. At f=0 the prior is immaterial. Read the physical prior column rather than treating that placeholder as a truth label.

Full reproduction instructions and the exact eight external dependencies are in ../REPRODUCING.md and ../PACKAGE_DEPENDENCIES.json. The pre-execution draft remains preserved as ../REPRODUCING_DRAFT.md. Engine and PROTOCOL hashes remain those fixed in results/EXECUTION_PROTOCOL.json.

The independent verifier reconstructs all66 observation-view category arrays, all124 million sampled indices, templates,156,000 calibration profiles and216,000 test profiles/p-values/decisions,216 metric cells and144 cadence pairs. Its recorded numerical probabilities and decisions agree exactly; maximum score difference is3.82e−11. This verifies implementation consistency, not a new simulation trial or Gaia validation. Reproduce with `python work/dr4_chain/review/verify_decisions.py --root work/dr4_chain --out work/dr4_chain/review/decision_verification_reproduced.json`.

The final seed audit has59 passing checks, including exact declaration matching for that independent replay. Its first58-check pre-freeze outcome is preserved from the execution transcript in audit_history/PREFREEZE_RECORDED_EVIDENCE.json; the exact checker source is archived alongside it. The original mutable audit JSON had been refreshed, so that evidence record does not claim to recover its original complete bytes. The first post-review scanner signal is separately preserved and explained in audit_history/README.md.
