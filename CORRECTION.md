# Correction of catalogue provenance and scientific interpretation

Hicham Boufourou — 8 September 2026

The earlier paper and reproducibility package treated an 81,088-pair input catalogue as observed Gaia DR3 wide binaries. Matching the redistributed columns by the two source identifiers identifies them numerically with `Newton_dr3_MSMS_d200pc_5.csv`, the official Newtonian mock in [Chae's Zenodo record 10652994](https://doi.org/10.5281/zenodo.10652994).

All 36 non-identifier columns agree across the 81,088 pairs at relative tolerance 1e-9, zero absolute tolerance and matching missing values. This is a numerical identity of columns; the CSV serializations are not byte-identical. The official mock has MD5 `997069b5635200853896d005a238dde9` and SHA-256 `4ed46795d71b507d1e6a4de9e7030b56240a8c29952f2b222bf4e2bae61e05cd`.

The associated [verification report](corrected-2026-09-08/work/final_paper/review/provenance_verification.json) and [verification source](corrected-2026-09-08/work/final_paper/verify_provenance.py) are provided. Retained source identifiers, magnitudes and sky positions do not make simulated proper motions into observations.

**The affected observational measurements and the gravity-exclusion claims based on them are withdrawn.** Renaming an input file cannot repair this inference. The replacement does not restore a confirmation of Newton or a general exclusion of MOND.

The real catalogue with corrected RUWE values is the 81,880-pair file from [Chae's Zenodo record 10986733](https://doi.org/10.5281/zenodo.10986733). Its MD5 is `1b6c5063163a4e6c07043d13aeb70f55`, and its SHA-256 is `f8de60c31865beddcf6b3be0f5d1a3a8d8b705f1a3df9318bce21d7e46dcb4f7`. In the corrected experiment it supplies empirical covariates under a declared selection, not a measured gravitational truth.

The [replacement manuscript](corrected-2026-09-08/ARTICLE_FINAL.pdf) constructs a conditional simulation experiment with distinct fitting, calibration and test roles. It quantifies correct, wrong and indeterminate decisions, including failed calibration-transfer targets. Its Appendix A explains the change in scientific scope.

Historical files remain available in [legacy/af8ddf07](legacy/af8ddf07) and in the original commit [af8ddf0724c881d1a7720a8c4cb4255195ace4ac](https://github.com/HBoufourou/paperI-wide-binaries/tree/af8ddf0724c881d1a7720a8c4cb4255195ace4ac). This notice changes their scientific interpretation; it does not erase the historical record.
