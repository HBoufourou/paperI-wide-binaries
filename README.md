# Estimator forensics for the wide-binary gravity test

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22114321.svg)](https://doi.org/10.5281/zenodo.22114321)
[![arXiv](https://img.shields.io/badge/arXiv-2608.24556-b31b1b.svg)](https://arxiv.org/abs/2608.24556)

Hicham Boufourou — Independent Researcher, Brussels, Belgium

## Important data-provenance correction (2026-09-02)

Repository versions 1.0 and 2.0 incorrectly treated an 81,088-row synthetic
Newtonian catalogue (`Newton_dr3_MSMS_d200pc_5.csv`) as the observational Gaia
DR3 sample. The affected Section 5 measurements, figures, and claims — including
the quoted `gamma = 1.05`, `Delta ln L = 129`, approximately `16 sigma`
rejection, and the deep-bin exclusion — are withdrawn.

The synthetic injection–recovery experiment remains a synthetic result and was
rerun after the correction. The corrected application uses the complete
81,880-row real Gaia catalogue with corrected RUWE from Chae's Zenodo record
10986733. Its outcome is model-dependent and does not support the previous
strong observational verdict. See [CORRECTION_NOTICE.md](CORRECTION_NOTICE.md)
for the numerical audit and scope of the correction.

## Corrected reproduction

Python 3.10+ is recommended.

```bash
python -m pip install -r requirements.txt
python code/data_source.py
python code/audit_catalog_provenance.py
python code/moteur_population.py
python code/phase_G.py
python code/corrected_dr3_analysis.py
python code/build_pdf.py
```

`data_source.py` downloads the authoritative table, checks its MD5 digest and
row count, and refuses the known `Newton_*` filenames and the 81,088-row
fingerprint. Downloaded external data are not committed.

## Repository map

- `paperI.tex`, `manuscript_content.tex` — corrected arXiv/preprint source.
- `paperI_mnras.tex` — MNRAS wrapper using the same manuscript body.
- `code/data_source.py` — verified real-catalogue loader.
- `code/audit_catalog_provenance.py` — independently reproduces the provenance finding.
- `code/corrected_dr3_analysis.py` — corrected observational analysis and figures.
- `code/phase_G.py` — corrected synthetic injection–recovery grid.
- `code/moteur_population.py` — population and hierarchical-triple engine.
- `results/` — machine-readable audit and corrected numerical results.
- `data/legacy_*`, `code/legacy/`, `figures/legacy/`, and `legacy/` —
  quarantined v1/v2 inputs, scripts, figures, and manuscript; retained only to
  reproduce the error, never to measure the sky.
- `protocole/PREREGISTRATION_DR4.md` — original frozen DR4 protocol.
- `protocole/PREREGISTRATION_DR4_AMENDMENT_2026-09-02.md` — transparent amendment.

## Version history

- **v3.0-correction (prepared 2026-09-02)** — provenance audit, corrected real
  Gaia input, reanalysis, withdrawal of affected observational conclusions.
- **v2.0 (2026-08-26)** — archived at
  [10.5281/zenodo.22114321](https://doi.org/10.5281/zenodo.22114321); affected by
  the catalogue error described above.
- **v1.0 (2026-07-15)** — archived at
  [10.5281/zenodo.22073431](https://doi.org/10.5281/zenodo.22073431); affected by
  the same error.

The Zenodo concept DOI is
[10.5281/zenodo.22073430](https://doi.org/10.5281/zenodo.22073430). The v3 DOI
must be inserted here after the new Zenodo version is published.

## Computational tools and responsibility

Large language models were used as computational and analytical assistants to
write, debug, and audit code and text. Numerical claims are tied to committed
scripts and machine-readable outputs. The author remains responsible for the
scientific assumptions, interpretation, and publication decisions.

## Licence

See [LICENSE](LICENSE).
