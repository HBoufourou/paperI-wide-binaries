# Conditional gravity decisions for wide binaries

**Corrected scientific version — 8 September 2026**  
Hicham Boufourou

The current manuscript is **Conditional gravity decisions for wide binaries: orbital populations, shared astrometry and calibration-library transfer**.

**Correction:** the 81,088-pair catalogue analysed as Gaia observations in the earlier paper is numerically identified as Chae's official Newtonian simulation number 5. The affected observational interpretation and gravity-exclusion claims are withdrawn. This replacement studies conditional decisions in simulations; it reports no Newton–MOND verdict from the real sky. See the [correction and provenance evidence](CORRECTION.md).

## Read the corrected work

- [Corrected manuscript — PDF](corrected-2026-09-08/ARTICLE_FINAL.pdf)
- [Manuscript source — Markdown](corrected-2026-09-08/ARTICLE_FINAL.md)
- [Scientific files and reproduction instructions](corrected-2026-09-08/REPRODUCIBILITY.md)
- [Current publication and data availability](PUBLICATION_STATUS.md)

The corrected directory contains the manuscript, frozen protocols, inference and population code, detailed result tables, and verification records. It is the active scientific version in this repository. Large physical banks and saved calibration/test arrays require the accompanying data package; the new Zenodo deposit is still pending. The older Zenodo DOIs do **not** identify this corrected experiment.

## What the experiment establishes

The pipeline links declared Newtonian and QUMOND orbital populations, simulated shared astrometry, nuisance profiling, and separate fitting, calibration and test libraries. Five matched variants and two calibration rules are compared on the same catalogue draws. There are 702,000 calibration catalogue views and 540,000 test views; these reuse finite physical libraries and are not independent new stellar populations.

For Full66 with the calibration envelope:

| Criterion across 540 test conditions | Conditions passing |
| --- | ---: |
| Observed classification targets | 540 / 540 |
| Classification targets using Wilson bounds | 531 / 540 |
| True-hypothesis retention using Wilson bounds | 501 / 540 |
| Both strengthened criteria together | 501 / 540 |

The minimum observed correct-decision rate is 98.1% and the maximum wrong-decision rate is 1.0% in this synthetic experiment. These are recovery rates for a simulated truth, not probabilities that Newton or MOND is true. The minimum true-hypothesis retention is 98.1%, below the 99% target. No tested variant and rule passes every strengthened target.

The later simulated astrometric information package improves separation in the matched comparison. Diagnostics do not show a uniform advantage, and the calibration envelope can turn previously indeterminate cases into wrong exclusive decisions. These failures are retained in the manuscript and tables. The method is a reproducible preparation for future tests; calibration of real contamination, selection and Gaia products remains necessary before a sky verdict or a claim that DR4 will be decisive.

## Earlier versions

The previous manuscript, scripts, reductions and figures are preserved under [legacy/af8ddf07](legacy/af8ddf07), with their original bytes and Git history. They are historical evidence, not the active scientific results. The original [DR4 protocol](protocole/PREREGISTRATION_DR4.md) is retained unchanged; read its [dated status note](protocole/STATUS_2026-09-08.md) alongside it. The earlier draft [pull request #1](https://github.com/HBoufourou/paperI-wide-binaries/pull/1) predates the present expanded experiment and is not the active manuscript.

The original arXiv entry is [arXiv:2608.24556](https://arxiv.org/abs/2608.24556); the corrected manuscript here is a separate version and does not imply that the arXiv record has been updated or that MNRAS has accepted the replacement.

## Attribution and citation

See [CITATION.cff](CITATION.cff) for this dated repository version and [licensing and attribution](LICENSES_AND_ATTRIBUTION.md). AI assistance was used in development and checking; numerical reconstructions are not represented as external human peer review. Responsibility for the scientific content remains with the author.
