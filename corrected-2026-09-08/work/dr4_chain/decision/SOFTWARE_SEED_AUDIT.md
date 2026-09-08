# Software-only execution namespace correction before any physical outcomes

The first full plumbing fixture used artificial48-row input banks and eight calibration repeats. It correctly exercised the interface and returned indeterminate throughout, but the fixture execution specification inherited the not-yet-used scientific catalogue seed numbers270–275M. This was found immediately after that software-only run, before any physical library generation, physical calibration or physical test result.

The complete initial fixture, initial engine source and initial fixture driver are preserved under software_fixture_initial_seed_audit/. It is an implementation-history artifact, not physical evidence and not a new scientific library. No old scientific result was changed or relabelled.

Fixture executions now use279100001+1000*calibration_case_index and279300001+1000*test_case_index. Artificial input rows use279000010; isolated unit fixtures use279000001–3. All are distinct from the124 prescribed scientific streams270–275M and upstream population/observer/trajectory streams260–269M. The corrected software fixture is rerun in software_fixture/ and its declaration records these separate seeds. Scientific seed rules, scores, thresholds, nuisance grid and observation inputs are unchanged.
