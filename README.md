# Reproducibility package — Boufourou (2026), Paper I
## "Estimator forensics for the wide-binary gravity test: the eccentricity–triple coupling manufactures a pseudo-signal, and a pre-registered protocol for Gaia DR4"
Hicham Boufourou — Independent Researcher, Brussels, Belgium — hicham.boufourou@hotmail.com

## Contents
- `Boufourou_2026_PaperI_wide_binaries.pdf` — the article (preprint).
- `code/` — the full analysis pipeline:
  - `moteur_population.py` — synthetic wide-binary population engine (Keplerian orbits,
    El-Badry-type selection, Gaia-like noise, triple contamination).
  - `etape_AB_reduction.py` — independent reduction of the Chae (2024) sample
    (v_2D, v_c, quality cuts).
  - `etape_C_robustesse.py` — robustness grid over cuts and estimator variants.
  - `etape_D_verdict.py` — pre-registered decision criteria applied to the data.
  - `phase_G.py`, `phase_H.py` — injection–recovery experiments: the competing estimator
    families confronted with the same synthetic universes of known truth (Newtonian and
    MOND-EFE), producing the bias map and the eccentricity–triple coupling measurement.
  - `figs_final.py` — regenerates all figures.
  - `build_pdf.py` — regenerates the article PDF.
- `data/` — `chae_colonnes_utiles.csv`, `chae_complement.csv` (columns extracted from the
  public Chae 2024 sample, 81,088 pairs, Gaia DR3 / El-Badry et al. 2021),
  `vtilde_reduit_v1.csv` (independent reduction), `phase_G_resultats.csv` (injection–recovery
  results grid).
- `figures/` — all article figures, including the master figure (six synthetic universes).
- `protocole/PREREGISTRATION_DR4.md` — the frozen, pre-registered protocol for Gaia DR4:
  estimators, cuts, decision criteria and STOP rule, written in advance of the DR4 release.

## Reproduction
Python 3 with numpy/scipy/matplotlib. Run, in order:
`etape_AB_reduction.py` → `etape_C_robustesse.py` → `phase_G.py` → `phase_H.py` →
`etape_D_verdict.py` → `figs_final.py` → `build_pdf.py`.
Each script prints its numbers; figures and the PDF are written alongside.

## Data provenance
The parent sample is the public wide-binary catalogue of Chae (2024, ApJ 960, 114),
drawn from El-Badry et al. (2021, MNRAS 506, 2269), Gaia DR3. Only derived columns
required for the analysis are redistributed here.
