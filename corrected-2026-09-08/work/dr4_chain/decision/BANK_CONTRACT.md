# Input manifest contract

All paths are relative to the directory containing the manifest. Absolute input paths are accepted only as supplied metadata; a portable package should use relative paths. Hashes are required before files are read. The engine owns neither the physical generator nor the observer. Exactly22 bank records are required: both gravities crossed with fitting3/calibration3/test5 components, each record supplying the3 paired scenario views.

```json
{
  "schema_version": 1,
  "purpose": "physics_production",
  "production_authorized": true,
  "physical_validation": {
    "status": "passed_for_declared_domain",
    "path": "physical_validation.json",
    "sha256": "provided by upstream"
  },
  "physical_scope": "The explicitly validated potential and selection domain, with limitations",
  "scenarios": ["uniform_white", "clustered_correlated", "sparse_floor"],
  "records": [
    {
      "gravity": "Newton",
      "role": "fitting",
      "component": "pure",
      "prior_au": null,
      "population_path": "banks/newton_fitting_pure/population.npz",
      "observations": {
        "uniform_white": "banks/newton_fitting_pure/uniform_white.npz",
        "clustered_correlated": "banks/newton_fitting_pure/clustered_correlated.npz",
        "sparse_floor": "banks/newton_fitting_pure/sparse_floor.npz"
      },
      "sha256": {
        "population": "SHA256",
        "observations": {
          "uniform_white": "SHA256",
          "clustered_correlated": "SHA256",
          "sparse_floor": "SHA256"
        }
      },
      "metadata": {"seeds": [], "admissibility": "upstream common selection", "source_lineage": "upstream"}
    }
  ]
}
```

Gravity strings are exactly Newton/QUMOND. Roles are fitting/calibration/test. Components are pure,comp40,comp80,comp60,comp160; the last two exist only for test. Each fitting/calibration role therefore has6 records; test has10. All records need weight_beta1 and weight_beta1p6. Test records also need weight_beta1p3. Weight arrays need not be normalized on disk, but must be finite, nonnegative, have positive sum and match the row count. They are normalized separately per admissible component. All observed arrays correspond to the same row order in that record, and all three scenarios share the identical population.

Population and observation required array names/shapes are listed in PROTOCOL.md and checked in engine.measured_features. Additional true positions, energies, companion parameters, original velocities, gravity labels and bias arrays may exist for upstream audit but are never read by the measured feature function. Component/gravity metadata and beta weights belong only to the simulator distribution construction, supervised templates and evaluation labels.

The switch production_authorized is an upstream permission/status assertion, not proof in itself. Before freezing a production run the caller must reserve the decision seeds and review the referenced numerical validation for the declared physical domain. The engine is not authorized to turn an incomplete physical calculation into an accepted model. Software fixtures must instead set purpose=software_fixture and be invoked explicitly with --software-fixture; their output must not be reported as physical validation.

Commands, once banks are authorized and all hashes are final:

```powershell
python work/dr4_chain/decision/engine.py freeze --manifest work/dr4_chain/bank_manifest.json --out work/dr4_chain/decision/results
python work/dr4_chain/decision/engine.py calibrate --manifest work/dr4_chain/bank_manifest.json --out work/dr4_chain/decision/results
python work/dr4_chain/decision/engine.py test --manifest work/dr4_chain/bank_manifest.json --out work/dr4_chain/decision/results
```

Every source/protocol/input hash is fixed at freeze. The engine refuses existing output stages. For a relocated reproduction use an intact relative-path manifest or a separately recorded new manifest and a fresh execution directory; do not replace hashes inside archived results. Computational stages use only Python and NumPy (CSV output uses the standard library).
