# Conditional information-block and calibration-library experiment

The completed execution used 702,000 paired calibration catalogue views and 540,000 test views, each evaluated in five fixed information variants. Each test cell contains 1,000 catalogue draws of 1,000 parents. Baseline and envelope rules act on those same scores and catalogues. The10 method/rule rows below each summarize540 cells; they are not10 independent physical experiments.

All previous releases remain unchanged. This protocol was fixed after inspection of the earlier development experiment. It expands the beta grid prospectively, retains the original fitting banks, adds two calibration-library realizations and reserves two fresh test-library realizations. It is not a blind human study, a preregistration or a validation of the actual Gaia sky.

## Separate classification and retention criteria

| Mode | Rule | Cells | Classification_count | Retention_count | Classification_CI | Retention_CI | Both_CI | Worst_correct | Worst_wrong | Worst_true_retained |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| full66 | baseline | 540 | 540 | 532 | 532 | 488 | 488 | 97.5% | 0.9% | 97.5% |
| no_diagnostics66 | baseline | 540 | 540 | 534 | 540 | 469 | 469 | 98.5% | 0.1% | 98.5% |
| full34 | baseline | 540 | 536 | 533 | 489 | 449 | 439 | 92.6% | 0.9% | 98.6% |
| covariates_only66 | baseline | 540 | 0 | 506 | 0 | 375 | 0 | 0.0% | 2.2% | 97.7% |
| conditional66 | baseline | 540 | 528 | 535 | 502 | 501 | 480 | 90.1% | 2.5% | 97.3% |
| full66 | envelope | 540 | 540 | 537 | 531 | 501 | 501 | 98.1% | 1.0% | 98.1% |
| no_diagnostics66 | envelope | 540 | 540 | 538 | 540 | 493 | 493 | 98.7% | 0.1% | 98.7% |
| full34 | envelope | 540 | 535 | 538 | 495 | 485 | 468 | 92.7% | 0.7% | 98.8% |
| covariates_only66 | envelope | 540 | 0 | 528 | 0 | 436 | 0 | 0.0% | 1.4% | 98.4% |
| conditional66 | envelope | 540 | 528 | 537 | 508 | 514 | 492 | 90.1% | 2.1% | 97.9% |

Classification count-pass requires observed correct≥95% and wrong≤1%. Retention count-pass separately requires true-model retention≥99%. The stronger CI gates require Wilson95% lower(correct)≥95%, upper(wrong)≤1%, and lower(true-retained)≥99%, with classification and retention kept separate. Both-CI requires both gates. These are marginal descriptive criteria; they are not simultaneous confidence statements over all cells or selected extrema.

There are 1,145 method/rule/cell rows failing at least one count criterion and 1,570 failing at least one stronger confidence criterion. Every failure is retained in count_failures.csv or confidence_failures.csv. metrics.csv contains all5,400 rows; family_summary.csv separates fitted-family, outside-inner-prior and outside-beta-grid tests; realization_summary.csv preserves the two test realizations. worst_cells.csv supplies the exact parameters and intervals behind each reported extreme. A good average never replaces a failed cell.

## Paired effect of the calibration-library envelope

| Mode | Correct_increased | Correct_decreased | Wrong_increased | Wrong_decreased | Retention_increased | Correct_delta_range | Wrong_delta_range |
| --- | --- | --- | --- | --- | --- | --- | --- |
| full66 | 61 | 0 | 15 | 0 | 61 | [0.00, 1.00] pp | [0.00, 0.10] pp |
| no_diagnostics66 | 85 | 0 | 0 | 0 | 85 | [0.00, 0.90] pp | [0.00, 0.00] pp |
| full34 | 17 | 212 | 3 | 66 | 113 | [-3.80, 0.40] pp | [-0.50, 0.10] pp |
| covariates_only66 | 1 | 520 | 13 | 268 | 315 | [-22.20, 0.10] pp | [-1.60, 0.20] pp |
| conditional66 | 39 | 35 | 0 | 33 | 62 | [-1.30, 0.60] pp | [-0.40, 0.00] pp |

An increased wrong rate is adverse, even though its signed change is positive. All rates above compare exactly paired catalogues; envelope_effects.csv retains all2,700 per-cell/mode contrasts, recovered/lost correct and wrong decisions, true-model retention changes and transitions from double rejection. paired.csv also retains all prescribed information-block comparisons under both calibration rules. The tables do not treat the three cadences as independent replicates.

The envelope is the maximum across libraries of each already Bonferroni-combined hypothesis p-value. It never lowers an individual p-value, so lost true-model retention must be zero; this identity was checked in the descriptive reconstruction. That does not imply wrong classifications can never increase: an earlier double rejection may become a single wrong retention. The observed paired transitions, rather than an assumed superiority, determine the interpretation.

## What the information comparisons identify

| Rule | Comparison | Correct_increased | Correct_decreased | Correct_delta_range | Wrong_delta_range |
| --- | --- | --- | --- | --- | --- |
| baseline | full66 minus no_diagnostics66 | 146 | 98 | [-1.90, 1.20] pp | [-0.10, 0.90] pp |
| baseline | full66 minus full34 | 402 | 11 | [-0.50, 7.10] pp | [-0.90, 0.80] pp |
| baseline | full66 minus covariates_only66 | 540 | 0 | [40.70, 100.00] pp | [-2.20, 0.70] pp |
| baseline | full66 minus conditional66 | 299 | 10 | [-0.20, 9.60] pp | [-1.60, 0.20] pp |
| envelope | full66 minus no_diagnostics66 | 122 | 103 | [-1.50, 1.00] pp | [-0.10, 1.00] pp |
| envelope | full66 minus full34 | 424 | 9 | [-0.40, 7.10] pp | [-0.70, 0.90] pp |
| envelope | full66 minus covariates_only66 | 540 | 0 | [42.80, 100.00] pp | [-1.40, 0.90] pp |
| envelope | full66 minus conditional66 | 318 | 11 | [-0.30, 9.60] pp | [-1.10, 0.30] pp |

Each contrast retains all540 paired cells; information_summary.csv includes positive, zero and negative differences plus first-only/second-only discordant totals for correct, wrong and indeterminate outcomes. A positive wrong-rate difference is adverse. These summaries do not replace the exact-cell paired.csv or transform reused cadences into independent observations.

full66 versus no_diagnostics66 removes exactly the diagnostic block while retaining the same projected velocity and uncertainty/separation cells. full34 instead uses its own measured PM/covariance and only34-month source residual/acceleration diagnostics; its interface never reads66-month or cross-release information. That contrast concerns the whole later information package, not observation duration alone, and full34 is not a recreation of the public Gaia DR3 products.

covariates_only66 retains normalized separation and formal66-month precision. Separation is a physical outcome and selection variable. The control therefore measures information supplied by the stipulated parent and selection; successful classification by that control is not itself label leakage. conditional66 subtracts the marginal-z deviance for each nuisance point before profiling, using the full mixed distribution. It removes the direct marginal likelihood term while retaining differences in available stratum counts and conditional populations. Its calibration is not exactly conditional on each possible observed n_z, and it does not prove that parent assumptions have been eliminated.

The fitted grid is beta1/1.3/1.6 and the comp40/comp80 mixture family. The new1.15/1.45 cases are outside that discrete grid, and actual60/160 AU inner priors remain physical transfer tests. Pure systems are counted once per beta/hypothesis/realization and have no meaningful inner-prior or mixture-eta truth. Eta is null in test records; templates use eta only to mix40/80.

The f=0.15 points and eta=0.5 mixtures enter fitting/profiling and calibration, but have no dedicated test-library cells. Calibration coverage of those nuisance points is not an independent transfer test of them.

## Finite-library and physical limits

Rank calibration applies under its calibrating empirical law. The library envelope covers a finite union of those laws, not every fresh empirical library or continuous nuisance value. The two new test realizations measure transfer and provide a limited sensitivity check. They do not measure fitting-library variability: the original six fitting banks stay fixed. Repeated catalogues reuse finite banks; counts, effective row numbers and smoothing do not certify a larger sample size or a future-survey guarantee.

Across the44 libraries and five stored beta weight sets, the minimum effective row count is 16320.7 and the largest normalized parent weight is 0.000193509. These are global weight-concentration diagnostics, not binwise support guarantees or counts of independent repeated catalogues. library_support.csv and fitting_component_support.csv retain the source records; every metric row separately records the mean fraction entering cells empty in the raw fitting union. Near-best nuisance counts use a descriptive2-deviance window and are not parameter confidence intervals.

The same scalar two-mass Newton/QUMOND model, fixed synthetic external field, restricted stationary parent, empirical covariates, photocentre prescription and observation operator are inherited. True total mass is tied to a fixed photometric anchor before orbital scaling. The inherited empirical pool already reflects chance-alignment, RUWE, precision, RV and transverse-velocity selection. Source-level error/photometry relations, distance/parallax uncertainties and real-survey selection remain conditional approximations. The unresolved companion enters the outer interaction through a monopole approximation; this is not a complete QUMOND three-body solution.

All five variants use measured fields only. Gravity/companion labels and energy weights define supervised simulated laws and evaluation truth, never the individual measured feature vector. The numerical results require an independent implementation check; that review is a separate artifact and is not external human peer review. The public-sky readiness status stays false.

## Reproducibility

STATISTICAL_DECLARATION.json fixes sources, protocol,414 catalogue streams and all cases before execution. EXECUTION_PROTOCOL.json adds immutable input hashes and cached measured codes. Every calibration/test case saves one parent-index array shared by all methods/cadences, fullfloat64 nuisance surfaces, scores and decisions. Counts are reconstructible from those indices and codes and have canonical hashes; they are not duplicated on disk. The independent checker can reconstruct every case.

Recorded stage durations (seconds, descriptive single-thread measurements): {"calibration": 527.4323837000411, "freeze": 36.675012699794024, "test": 524.2472886999603}. Code hashes, data hashes and any reproduction-run timing differences must be retained rather than edited into an old manifest.
