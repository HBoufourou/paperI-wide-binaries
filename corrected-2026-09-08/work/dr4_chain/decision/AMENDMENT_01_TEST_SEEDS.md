# Before-production amendment: reserve the executed analytic fixture

8 September 2026. This amendment was approved before any physical decision freeze, calibration or test. Physical-grid refinement was still in progress. No physical decision outcome was available.

The independent analytic action-integral fixture in `../review/check_energy.py` had already used seed `275000001`. The original prospective QUMOND test namespace overlapped that seed. Preserve the fixture, its source and its executed evidence. Reserve `275000001` for that historical analytic fixture; it is not a scientific population or decision catalogue.

Only future QUMOND decision-test seeds change: the test rule is now `274000001 + 3000000*gravity_index + 10000*beta_index + 1000*prior_index + 100*f_case_index`. Thus Newton tests remain at 274M and QUMOND tests use 277M. Calibration remains at 270/271M. All indices and declared cases retain the definitions in PROTOCOL.md. No score, bin, smoothing, nuisance grid, catalogue size, repeat count or decision threshold changes.

The earlier software fixture remains intact. Its exact executed source and protocol are additionally preserved under `software_fixture/executed_source/`; the archived initial plumbing run remains under `software_fixture_initial_seed_audit/`. These archives used artificial inputs and are not new physical evidence. In particular, the initial plumbing run's historical overlap with prospective scientific seed numbers is disclosed in SOFTWARE_SEED_AUDIT.md, rather than erased or described as independent physical sampling.

`seed_audit.py` independently enumerates all 124 future decision streams, the 44 upstream population/observer streams and the declared upstream validation/analytic streams. It compares the enumeration with the current decision case factories, checks upstream source expressions and scans current root/review/physics Python sources for additional large integer seed candidates. The existing analytic fixture must be present in the reservation audit and absent from future scientific streams. Source hashes make this an audit of a specific source snapshot; rerun the audit if upstream sources change before production.

The separate upstream amendment `../AMENDMENT_01_PHOTOMETRIC_ANCHOR.md` is also part of the current input conditions. It anchors total photometric mass before true-mass orbital scaling and does not change any decision rule. Neither amendment establishes a guarantee for real Gaia data.
