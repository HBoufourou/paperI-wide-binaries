# Statistical output contract

The controlling scientific specification is `../PROTOCOL.md`. This new implementation imports no old decision module. `features.py` reads only declared measured fields. `engine.py` legitimately uses population energy and gravity/companion labels only for simulation weights, supervised templates and known evaluation truth. Python and NumPy suffice.

All paths in bank manifests are relative to the manifest directory. A bank key is `(gravity, role, realization, component)`. Exactly44 records are required: six historical fitting records at realization0; six historical calibration records at0; twelve new calibration records at1/2; twenty new test records at1/2. The new physical manifest must be authorized and all referenced population/observation hashes must match. The independent physical verification and root approval precede statistical input freeze.

`STATISTICAL_DECLARATION.json` is written once by `engine.py declare`, before scientific draws. It fixes root-protocol SHA, engine/features SHA, cases, seeds and dimensions. Scientific repeat count is1000, catalogue size1000. A separate input freeze creates `results/EXECUTION_PROTOCOL.json`; replacing an existing declaration/stage is refused. If implementation work changes code before production, do not silently overwrite an existing frozen declaration.

The mode axis order is `full66, no_diagnostics66, full34, covariates_only66, conditional66`. Rule order is `baseline, envelope`. Scenario order is `uniform_white, clustered_correlated, sparse_floor`. Hypothesis order is `Newton, QUMOND`. Theta order is beta{1,1.3,1.6}, then ascending f{0,.075,.15,.225,.3}, then eta{0,.5,1}, retaining eta0 only at f0:39 points.

Calibration realization order is0,1,2. Calibration cases iterate realization, hypothesis, theta (234 cases). Test cases iterate realization1,2, hypothesis, beta{1,1.3,1.6,1.15,1.45}; each fitted beta first has the pure case, then f{.075,.225,.3} crossed with prior{40,80,60,160}; each outside beta first has pure then f.3 prior60/160. There are45 cases per hypothesis/realization,180 total. Pure is counted once, with `prior=null`. Test eta is `null`, avoiding a false physical interpretation of an unused placeholder.

The f=.15 points and eta=.5 mixtures are profiled and calibrated, but have no dedicated test-library transfer cells. They must not be described as independently tested merely because the calibration grid contains them.

Scientific calibration seed is `320000001+1000000*realization+100000*h+1000*beta_index+10*mixture_index`. Scientific test seed is `325000001+1000000*(realization-1)+100000*h+1000*local_case_index`. These414 seeds do not depend on mode, cadence or calibration rule. Software input uses310000010, isolated unit fixtures310000001–4, and full software-run cases use310100001+1000*global_case_index across calibration followed by test. The software override changes only repeats and case seeds in an explicitly software-only execution; its declaration remains distinct from scientific data.

## Cached measured codes and templates

`results/codes/bankXX.npz` stores `full66` and `full34` codes[scenario,Nparent],uint16, and weights for beta1/1.3/1.6/1.15/1.45. Historical beta1/1.3/1.6 arrays are checked against weight1×(-E)^(beta−1); test-only weights are calculated by that formula. Code cache rows and hashes are fixed in EXECUTION_PROTOCOL. IndexXX follows sorted bank keys; metadata records the exact source key and raw-file lineage.

Full codes encode `((((sep*8+parallel)*8+perpendicular)*4+diagnostic)*2+precision)`. The no-diagnostic projection is `(code//8)*2+code%2` (512 cells). The covariate projection is `(code//512)*2+code%2` (8 cells). Full34 has its own2048 codes built with its own measured PM/covariance and only source residual/acceleration scores34. It never requests a66/delta/cross field. The other reduced modes likewise accept only their minimum field sets; they do not first call full66 and discard prohibited inputs.

`templates.npz` holds probability arrays for full66/full34/no_diagnostics66/covariates_only66, each[scenario,hypothesis,theta,bin]. Marginal templates are exact projections of the already mixed, smoothed full66 probabilities. `raw_MODE` arrays[scenario,bin] report the fitting support union before smoothing. Conditional66 uses the full66 probability array and the covariate marginal rather than a separate componentwise conditional mixture. Template hash and component support are in template_manifest.json.

## Shared sampled indices and scores

Each calibration/test case has one NPZ `caseNNN.npz`, with one `parent_index`[repeat,1000],uint32. For that case, concatenate the participating components in `lineage.components` order; normalize weights separately within each component before applying the mixture fractions. `component_offsets` identify each parent block. The same indices address all modes and three paired observations. No count array is saved. Each case's JSON record contains canonical count hashes for full66,full34,no_diagnostics66,covariates_only66. Reconstruct by bincount of the saved indices and codes. Conditional66 shares full66 counts.

Hash a count array as SHA256 of sorted-key JSON containing dtype and shape, followed by contiguous little-endian array bytes. Counts have dtypeuint16 and shape[scenario,repeat,bin]. This is an array hash, not a ZIP-container hash.

Each case NPZ contains `surface`[mode,scenario,repeat,hypothesis,theta], `D`/`R`[mode,scenario,repeat,hypothesis], and `best`/`near_best_count` with the same shape. Surfaces are float64. D is the minimum per-hypothesis surface and R_h=D_h−D_other. Near-best counts use the descriptive2-deviance window and are not confidence intervals.

For conditional66, compute the full-mixture surface minus its marginal-z surface separately at every theta, then minimize. The conditional mixture is not built by mixing component conditionals, and the conditional minimum is not obtained by subtracting independently minimized full and marginal deviances. This score removes the direct marginal likelihood term while calibration still mixes varying n_z. It is not an exactly n_z-conditioned test.

`calibration/reference.npz` contains `sorted`[mode,calibration_library,scenario,hypothesis,theta,repeat,statistic], statistic orderD,R, sorted on repeat. Its manifest has234 cases,702000 catalogue views and3510000 score views. References are independently drawn from oldcal0 and freshcal1/2 under new seeds.

Test files additionally contain `p_marginal`[mode,library,scenario,repeat,hypothesis,statistic] after maximum over theta; `p_library`[mode,library,scenario,repeat,hypothesis] after Bonferroni combination; `p_joint`/`rejected`[rule,mode,scenario,repeat,hypothesis]; and `decision`/`reason`[rule,mode,scenario,repeat]. Baseline selects library0. Envelope maximizes the already combined p_library over libraries. Maximizing separately over libraries within each statistic is prohibited and is covered by a counterexample unit test.

Decision is−1/0/1 for indeterminate/Newton/QUMOND. Reason is0 both_retained,1 both_rejected,2 one_retained. Increasing p-values under the envelope reduces hypothesis rejection, but need not monotonically reduce wrong classifications: a previous double rejection can become a single wrong retention. Paired outcomes must be measured, not assumed.

## Tables

`metrics.csv` has5400 rows:540 physical case/scenario cells×5 modes×2 rules. It includes all counts, true-model retention, Wilson intervals, separate classification and retention gates, combined gates, raw-support diagnostics and nuisance summaries. No repeated prior labels for pure systems are created.

`paired.csv` has7020 rows:4320 information comparisons (four contrasts, two rules,three cadences,180cases) plus2700 envelope-versus-baseline comparisons (five modes,three cadences,180cases). All differences are first−second, and each correct/wrong/indeterminate outcome has first-only and second-only discordant counts. These are paired comparisons conditional on finite libraries; separate cadences are not treated as independent replicates.

Stage manifests record file hashes, count hashes, lineage and dimensions. Final scientific conclusions require the independent reconstruction and a report preserving every failure. No old result is overwritten or relabelled as fresh evidence.

## Conditional measured-input interface

`engine.predict_catalogue(population, observations, scenario, template, reference, mode='full66', rule='envelope', declared_scope_confirmed=False)` accepts measured arrays for exactly1000 parents. The caller supplies the matching validated template dictionary and sorted calibration reference. The selected feature whitelist is defined in `features.py`; latent gravity/energy/companion labels are not requested by this interface. Only the requested mode/scenario is returned.

The default scope assertion is false. Missing, nonfinite, malformed or unsupported measured inputs return `status='unsupported', decision='indeterminate'` with a reason. A valid, scope-confirmed call returns `status='conditional'`, the decision and reason, two hypothesis p-values, best nuisance indices, the fraction in raw-empty fitting cells, and `sky_ready=False`. Scope confirmation is a caller assertion, not a test of physical truth, real-Gaia calibration or correct survey selection. This API is a conditional replay mechanism, not an independently certified astronomical decision service.
