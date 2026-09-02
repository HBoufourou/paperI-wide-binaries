# Data-provenance correction and scope assessment

Date: 2026-09-02  
Affected public versions: arXiv:2608.24556v1; Zenodo v1/v2; repository commit
`af8ddf07` and earlier observational outputs.

## Finding

The v1/v2 scripts loaded a file named `Newton_dr3_MSMS_d200pc_5.csv` and the
repository redistributed an 81,088-row reduction of that file as though it were
the Chae (2024) observational sample. Chae's public record identifies
`Newton_dr3_MSMS_d200pc_5.csv` as a virtual Newtonian catalogue.

The independent audit in `code/audit_catalog_provenance.py` establishes that:

- all 81,088 identifiers match the official mock-v5 catalogue;
- every one of the 19 redistributed analysis columns matches that mock at
  relative tolerance `1e-9`;
- the median `abs(d1 - d_M)` is 0.000768 pc in the legacy file, versus 0.137569
  pc in the real Gaia catalogue;
- 94.19% of legacy distances lie within the projected orbital offset, versus
  8.42% for real Gaia;
- relative proper motions correlate with the real catalogue at only 0.626;
- the authoritative corrected-RUWE real catalogue contains 81,880 pairs.

The report is saved as `results/catalog_provenance_audit.json`.

## Consequence

The v1/v2 Section 5 computation did not measure gravity from the sky. It
measured the response of the estimators to a simulated Newtonian universe with
hidden companions. Therefore the following observational statements are
withdrawn, not merely updated:

- `gamma_test = 1.045 [1.025, 1.068]` on the claimed Chae observational sample;
- the joint-fit `gamma = 1.05 [1.04, 1.07]`, `f_trip = 0.15` interpretation;
- `Delta ln L = 129` and the derived approximately `16 sigma` statement;
- the claimed deep-bin exclusion of `gamma = 1.35--1.40`;
- the description of Gaia DR3 as decisively Newtonian in this work.

The numerical identity `sqrt(2 Delta ln L)` was not the problem; the likelihood
was evaluated on the wrong class of catalogue, so it cannot support an
observational significance claim.

## Corrected real-catalogue analysis

The new pipeline uses the complete 81,880-row
`gaia_dr3_MSMS_d200pc_ruwe.csv` from Zenodo record 10986733 and verifies MD5
`1b6c5063163a4e6c07043d13aeb70f55`. With the fiducial cuts, 38,472 pairs
remain; 28,221 are in the 0.2--2 kau validation zone, 7,900 in the 2--30 kau
test zone, and 403 in the deep bin.

At the pre-registered `T = sqrt(2)` truncation:

| eccentricity model | validation offset | E2, 2--30 kau (68%) | E2, deep bin (68%) | E3 best fit; f_trip |
|---|---:|---:|---:|---:|
| empirical per zone | +10.99% — FAIL | 1.084 [1.063, 1.105] | 1.308 [1.232, 1.418] | 1.13 [1.13, 1.14]; 0.14 |
| thermal | +1.77% — PASS | 1.007 [0.987, 1.023] | 1.221 [1.132, 1.316] | 1.01 [1.00, 1.02]; 0.16 |
| superthermal | +4.06% — FAIL | 1.007 [0.988, 1.022] | 1.222 [1.146, 1.314] | 1.04 [1.02, 1.05]; 0.15 |

The empirical branch, used for the original primary claim, fails the mandatory
2% validation control. The thermal branch passes and is consistent with Newton
in the global test zone, but the deep-bin estimates are higher and have broad
intervals. Results also vary materially with the eccentricity family. The
corrected outcome is therefore **model-dependent and non-conclusive** under the
pre-registered decision logic. No high-significance gravity verdict is claimed.

## Synthetic experiment

The known-truth injection–recovery experiment is logically separate and was
rerun with noise resampled from the corrected real catalogue. At residual
triple fractions 0.2--0.3, the full-median estimator returns 1.062--1.128 on
Newtonian truth; the truncated estimator returns 1.002--1.044 and the mixture
fit 0.98--1.00. For injected `gamma = 1.4`, the truncated and mixture
estimators return 1.306--1.346 and 1.36--1.40. This supports the narrower
methodological conclusion that tail treatment can bias a full-median
estimator. It does not determine the gravity law from Gaia DR3.

## RUWE and additional code corrections

The corrected release is used as a whole, rather than transplanting only its
RUWE columns. The new pipeline also propagates radial-velocity uncertainty in
the perspective correction and applies each trial `gamma` before noise and
selection to both the binary orbit and the low-acceleration outer orbit of the
triple template. Compact inner photocentre motion remains Newtonian.

## Reproducibility

Run `code/audit_catalog_provenance.py`, `code/phase_G.py`, and
`code/corrected_dr3_analysis.py`. Legacy inputs and scripts are explicitly
named and quarantined. The verified loader refuses both a `Newton_*` filename
and the 81,088-row fingerprint.
