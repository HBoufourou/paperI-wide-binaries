
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22073431.svg)](https://doi.org/10.5281/zenodo.22073431)
# Estimator forensics for the wide-binary gravity test

**arXiv:2608.24556** (astro-ph.GA, August 2026)  
Preprint: https://arxiv.org/abs/2608.24556  
PDF: https://arxiv.org/pdf/2608.24556
**Cite as:** Boufourou (2026), *Estimator forensics for the wide-binary gravity test*, Paper I v1.0, Zenodo DOI 10.5281/zenodo.22073431 — code & données : github.com/HBoufourou/paperI-wide-binaries
>>>>>>> 0b29cb3c014ddfc7cb4fe232168cfec60c9538ff
# Reproducibility package — Boufourou (2026), Paper I

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22114321.svg)](https://doi.org/10.5281/zenodo.22114321)
[![arXiv](https://img.shields.io/badge/arXiv-2608.24556-b31b1b.svg)](https://arxiv.org/abs/2608.24556)

## "Estimator forensics for the wide-binary gravity test: the eccentricity–triple coupling manufactures a pseudo-signal, and a pre-registered protocol for Gaia DR4"

Hicham Boufourou — Independent Researcher, Brussels, Belgium — hicham.boufourou@hotmail.com

---

## Versions

- **v2.0 (2026-08-26)** — updated reproducibility package
  - Zenodo: [DOI 10.5281/zenodo.22114321](https://doi.org/10.5281/zenodo.22114321)
  - arXiv: [arXiv:2608.24556](https://arxiv.org/abs/2608.24556)
  - MNRAS: manuscript MN-26-2659-P (under review)
  - **Updates since v1.0:**
    - Threshold sensitivity scan (`code/etape_E_seuil.py`) showing stability of γ across T ∈ [1.2, 1.8]
    - Complete Table 1 with deep-bin measurements
    - Makarov (2026, AJ 171, 79) citation as independent Newtonian verdict
    - Explicit definition of "oracles" terminology (Section 3)
    - Acknowledgments section
    - MNRAS-formatted LaTeX source (`paperI_mnras.tex`)

- **v1.0 (2026-07-15)** — frozen protocol version
  - Zenodo: [DOI 10.5281/zenodo.22073431](https://doi.org/10.5281/zenodo.22073431)
  - Pre-registered protocol for Gaia DR4 analysis (unchanged)

**Cite as:** Boufourou (2026), *Estimator forensics for the wide-binary gravity test*, Paper I v2.0, Zenodo DOI 10.5281/zenodo.22114321 — code & données : github.com/HBoufourou/paperI-wide-binaries

---

## Contents

- `Boufourou_2026_PaperI_wide_binaries.pdf` — the article (preprint, arXiv version).
- `paperI.tex` — LaTeX source (arXiv version).
- `paperI_mnras.tex` — LaTeX source (MNRAS format, submitted 2026-08-26).
- `code/` — the full analysis pipeline:
  - `moteur_population.py` — synthetic wide-binary population engine (Keplerian orbits,
    El-Badry-type selection, Gaia-like noise, triple contamination).
  - `etape_AB_reduction.py` — independent reduction of the Chae (2024) sample
    (v_2D, v_c, quality cuts).
  - `etape_C_robustesse.py` — robustness grid over cuts and estimator variants.
  - `etape_D_verdict.py` — pre-registered decision criteria applied to the data.
  - `etape_E_seuil.py` — threshold sensitivity scan: stability of γ across T ∈ [1.2, 1.8] (Section 5).
  - `phase_G.py`, `phase_H.py` — injection–recovery experiments: the competing estimator
    families confronted with the same synthetic universes of known truth (Newtonian and
    MOND-EFE), producing the bias map and the eccentricity–triple coupling measurement.
  - `figs_final.py` — regenerates all figures.
  - `build_pdf.py` — regenerates the article PDF.
- `data/` — `chae_colonnes_utiles.csv`, `chae_complement.csv` (columns extracted from the
  public Chae 2024 sample, 81,088 pairs, Gaia DR3 / El-Badry et al. 2021),
  `vtilde_reduit_v1.csv` (independent reduction), `phase_G_resultats.csv` (injection–recovery
  results grid), `etape_E_seuil_resultats.csv` (threshold scan results).
- `figures/` — all article figures, including the master figure (six synthetic universes).
- `protocole/PREREGISTRATION_DR4.md` — the frozen, pre-registered protocol for Gaia DR4:
  estimators, cuts, decision criteria and STOP rule, written in advance of the DR4 release.

## Reproduction

Python 3 with numpy/scipy/matplotlib. Run, in order:
`etape_AB_reduction.py` → `etape_C_robustesse.py` → `etape_E_seuil.py` → `phase_G.py` → `phase_H.py` →
`etape_D_verdict.py` → `figs_final.py` → `build_pdf.py`.
Each script prints its numbers; figures and the PDF are written alongside.

## Data provenance

The parent sample is the public wide-binary catalogue of Chae (2024, ApJ 960, 114),
drawn from El-Badry et al. (2021, MNRAS 506, 2269), Gaia DR3. Only derived columns
required for the analysis are redistributed here.

## Acknowledgments

The author thanks Andrei Tokovinin for a careful reading of an early version and for pointing to the Makarov (2026) analysis, and the researchers who responded to the arXiv endorsement request. This work received no external funding.

## Declaration on computational tools

This work was carried out by the author alone, without institutional affiliation. Large language models were used as computational and analytical assistants throughout: to write and debug the numerical scripts, to run the scans and minimisations, to produce the figures, to check algebra and dimensional consistency, to search and summarise the literature, and to audit drafts for internal contradictions. Several distinct systems were used and cross-checked against one another; no single system's output was accepted without an independent numerical or analytical verification.

The author is solely responsible for the physical hypotheses, the interpretation of every result, the epistemic labels attached to each claim, and the decision of what to publish and what to withdraw. All quantitative statements in this paper are reproduced by the scripts of the reproducibility package; a reader who runs them obtains the numbers printed here without needing any language model.
