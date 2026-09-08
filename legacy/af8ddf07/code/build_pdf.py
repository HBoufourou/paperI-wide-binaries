#!/usr/bin/env python3
"""build_pdf.py — version publiable (preprint) de l'article, PDF A4."""
from reportlab.lib.pagesizes import A4
from reportlab.lib.units import cm
from reportlab.lib.enums import TA_JUSTIFY, TA_CENTER
from reportlab.platypus import (SimpleDocTemplate, Paragraph, Spacer, Image,
                                HRFlowable)
from reportlab.lib.styles import ParagraphStyle
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.lib import colors

for name, path in [("DejaVu", "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf"),
                   ("DejaVu-Bold", "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf"),
                   ("DejaVu-Italic", "/usr/share/fonts/truetype/dejavu/DejaVuSans-Oblique.ttf")]:
    pdfmetrics.registerFont(TTFont(name, path))

S = dict(
 title=ParagraphStyle("t", fontName="DejaVu-Bold", fontSize=15.5, leading=19,
                      alignment=TA_CENTER, spaceAfter=10),
 author=ParagraphStyle("a", fontName="DejaVu", fontSize=10.5, leading=14,
                       alignment=TA_CENTER, spaceAfter=2),
 date=ParagraphStyle("d", fontName="DejaVu-Italic", fontSize=9, leading=12,
                     alignment=TA_CENTER, spaceAfter=14, textColor=colors.grey),
 h1=ParagraphStyle("h1", fontName="DejaVu-Bold", fontSize=11.5, leading=15,
                   spaceBefore=12, spaceAfter=5),
 abs=ParagraphStyle("ab", fontName="DejaVu", fontSize=8.8, leading=11.8,
                    alignment=TA_JUSTIFY, leftIndent=1*cm, rightIndent=1*cm,
                    spaceAfter=6),
 body=ParagraphStyle("b", fontName="DejaVu", fontSize=9.6, leading=12.8,
                     alignment=TA_JUSTIFY, spaceAfter=6),
 cap=ParagraphStyle("c", fontName="DejaVu", fontSize=8.2, leading=10.5,
                    alignment=TA_JUSTIFY, spaceBefore=3, spaceAfter=10,
                    textColor=colors.HexColor("#333333")),
 ref=ParagraphStyle("r", fontName="DejaVu", fontSize=8.6, leading=11,
                    alignment=TA_JUSTIFY, spaceAfter=2, leftIndent=0.4*cm,
                    firstLineIndent=-0.4*cm),
)

VT = "ṽ"; GM = "γ"; SG = "σ"; DL = "Δ"
story = []
P = lambda t, st="body": story.append(Paragraph(t, S[st]))

P("Estimator forensics for the wide-binary gravity test: the eccentricity–triple "
  "coupling manufactures a pseudo-signal, and a pre-registered protocol for Gaia DR4", "title")
P("<b>Hicham Boufourou</b>", "author")
P("Independent Researcher, Brussels, Belgium — hicham.boufourou@hotmail.com", "author")
P("Preprint — July 2026 — Paper I of a series", "date")
story.append(HRFlowable(width="100%", thickness=0.6, color=colors.grey, spaceAfter=8))

P("<b>ABSTRACT</b> — Analyses of Gaia DR3 wide binaries reach opposite verdicts on a claimed "
  f"low-acceleration gravitational anomaly: boost factors {GM} = G<sub>eff</sub>/G ≈ 1.4 "
  "(Chae 2023–2026; Hernandez et al. 2024) versus Newtonian gravity at high significance "
  "(Banik et al. 2024; Pittordis et al. 2025), from the same parent catalogue. We localize the "
  "disagreement experimentally. A generative population model — Keplerian binaries plus "
  "hierarchical triples with photocentre and time-averaging effects — is validated without "
  f"tuning against the published statistics of Pittordis et al. (2025): triple median {VT} "
  "reproduced to ≤0.04 in all eight separation bins, cut-survival 40.8% vs 40.5%. Synthetic "
  "catalogues of known truth, including the intrinsic velocity truncation of the El-Badry et al. "
  "(2021) selection and measurement noise resampled from the real data, are then confronted "
  "with three estimator families. A full-median estimator self-calibrated on the Newtonian "
  f"regime, even with a <i>perfect</i> eccentricity correction, recovers {GM} = 1.08–1.13 from "
  "purely Newtonian universes containing 20–30% residual triples — half-way to the claimed "
  f"anomaly — while a tail-truncated estimator and a mixture likelihood recover {GM} ≤ 1.04 "
  f"and {GM} = 1.00. An injected MOND-like boost ({GM} = 1.4 above 2 kau) is recovered at "
  f"{GM} ≥ 1.26 by all three: the hypotheses never overlap under robust estimation, so the "
  "controversy is decidable. The exercise exposed two noise pathologies relevant to all "
  "median-based analyses: a Rice-type bias requiring zone-matched noisy templates, and "
  "parallax-noise injection when perspective corrections use individual distances. Applying "
  "the calibrated estimators to the 81,088-pair Chae (2024) sample we measure "
  f"{GM}<sub>test</sub> = 1.045 [1.025, 1.068] (2–30 kau; f(e) systematic band 0.96–1.05) and, "
  f"from a joint ({GM}, f<sub>trip</sub>) mixture fit, {GM} = 1.05 [1.04, 1.07] with "
  "f<sub>trip</sub> = 0.15, consistent with the independent 0.17 of Pittordis et al. (2025); "
  f"{GM} = 1.4 is rejected at {DL}ln L = 129 (≈16{SG}). The deep-MOND bin "
  f"(g<sub>N</sub> &lt; 0.3 a<sub>0</sub>) is perspective-sensitive — {GM} moves from 0.97 to "
  "1.11 ± 0.07 after the angular correction, partly radial-velocity noise injection — and its "
  "full systematic bracket remains consistent with Newton while excluding 1.35–1.40. We freeze "
  "a pre-registered protocol (cuts, zones, both calibrated estimators, mandatory f(e) band, "
  "decision criteria written in advance) to be time-stamped before the Gaia DR4 release, where "
  f"{SG}<sub>stat</sub>({GM}) ≲ 0.005 will make estimator systematics fully dominant.", "abs")
story.append(HRFlowable(width="100%", thickness=0.6, color=colors.grey, spaceAfter=8))

P("1. Introduction", "h1")
P("Wide binaries probe gravity at internal accelerations below a<sub>0</sub> ≈ 1.2×10<super>−10</super> "
  "m s<super>−2</super>, where MOND-type theories with the Galactic external field effect (EFE) predict "
  f"an effective boost {GM} ≈ 1.35–1.40 of the gravitational constant, i.e. +18% on relative velocities "
  "(√1.4 = 1.183), while GR + dark matter predicts pure Newtonian dynamics (Banik & Zhao 2018; "
  "Pittordis & Sutherland 2018). Gaia DR3 made the test statistically possible — and produced a stark "
  "contradiction. Chae (2023, 2024), Chae et al. (2026) and Hernandez et al. (2024) report the anomaly "
  "at up to ~5σ; Banik et al. (2024, 19σ for Newton), Pittordis & Sutherland (2023), Pittordis et al. "
  "(2025, hereafter PS25) and the Quality Framework of Banik et al. (2026) find no deviation — from the "
  "same El-Badry et al. (2021) parent catalogue. The disagreement is therefore purely methodological — and the most recent exchange sharpens it: Chae & Yoon (2026) revisit precisely the two contested levers, data quality control and multiple-star modeling, and reconfirm the anomaly from the same catalogue. Independently of the statistical route, Pasquini et al. (2026) find that three of twelve rigorously vetted VLT-ESPRESSO systems admit no bound Newtonian orbit, keeping the per-system route open as well.")
P("Recent work has begun to localize it. On the 36-pair 3D sample, hierarchical Bayesian reanalysis of Saad & Ting (2026) showed the answer flips between {GM} = 1.12 and {GM} = 1.56 with the deprojection "
  f"choice alone. Banik et al. (2026) argued qualitatively that a loose upper limit on {VT} can generate "
  "a MOND-like signal, and PS25 criticized the calibration of the triple fraction at small separations. "
  "What has been missing is a controlled experiment: the competing estimator families confronted with the "
  "<i>same</i> synthetic universes of known truth, with the coupling between eccentricity modelling and "
  "tail treatment quantified. This paper provides that experiment (Sect. 4), an independent robust "
  "measurement on DR3 (Sect. 5), and a protocol frozen before Gaia DR4 (Sect. 6).")

P("2. Sample, observable and pipeline", "h1")
P("We use the public sample of Chae (2024): 81,088 pairs within 200 pc drawn from El-Badry et al. (2021), "
  "with proper motions and uncertainties, distances, masses, RUWE, chance-alignment probability, per-pair "
  "inferred eccentricities, and Gaia radial velocities for at least one component in 71% of pairs. The "
  f"observable is {VT} = {DL}v<sub>⊥</sub>/√(G<sub>N</sub>M<sub>tot</sub>/s<sub>2D</sub>), with Newtonian "
  f"normalization: an observer does not know G<sub>eff</sub>, so under G<sub>eff</sub> = {GM}G<sub>N</sub> "
  f"the whole {VT} distribution is multiplied by √{GM}.")
P(f"{DL}v<sub>⊥</sub> is computed from proper-motion differences with an angular perspective correction: "
  "the systemic velocity (mean proper motion plus available systemic radial velocity) is projected at each "
  "component's position at the common mean distance, and the predicted perspective proper-motion difference "
  "is subtracted (El-Badry 2019). Two implementation traps are worth recording. First, missing radial "
  "velocities in the source catalogue are encoded as sentinels (−10,000/−20,000 km s<super>−1</super>) and must be "
  "masked. Second, using individual parallax distances in the projection injects parallax noise amplified "
  "by the systemic tangential velocity (~30 km s<super>−1</super>), overwhelming the ~0.5 km s<super>−1</super> orbital "
  "signal; the correction must be purely angular, at the common distance. Pairs lacking any radial velocity "
  "receive the tangential part of the correction only.")
P("Pre-registered zones: validation 0.2–2 kau (Newtonian regime; mandatory calibration), test 2–30 kau, and "
  "a deep bin g<sub>N</sub>/a<sub>0</sub> ∈ [0.03, 0.3]. Fiducial cuts: RUWE &lt; 1.4 (both components), "
  f"R<sub>chance</sub> &lt; 0.01, {SG}<sub>{VT}</sub> &lt; 0.10 (41,760 pairs); robustness variants use "
  "RUWE &lt; 1.2 (37,197) and an HRD lobster-body proxy. A pipeline control requires the validation-zone "
  "median to match the Newtonian template to &lt;2%; we measure +1.4%.")

P("3. Generative population model, validated on published oracles", "h1")
P("Binaries are Keplerian orbits with isotropic orientation, uniform-in-time phase and switchable f(e) "
  "(flat; thermal; superthermal α = 1.3 of Hwang et al. 2022; Tokovinin & Kiyaeva 2016; or empirical per "
  "zone). Hierarchical triples follow PS25: Kroupa present-day mass function, log-flat outer orbits, "
  "Offner et al. (2023) lognormal inner semi-major axes with the Tokovinin (2014) stability limit, "
  "photocentre–barycentre factor with resolved/unresolved logic at 1″, inner velocities time-averaged over "
  "the 34-month DR3 baseline, the apparent-mass bias, and simulated RUWE / image-peak / lobster cuts.")
P("Without any tuning, the model reproduces the published PS25 statistics: binary P90 = 0.933–0.941 "
  f"(oracle 0.94 ± 0.01); fractions {VT} ≥ 0.8 of 23.1/21.0/20.5% for flat/thermal/superthermal (oracle "
  "23.2/21.0/20.6%); triple survival of the quality cuts 40.8% (oracle 40.5%, RUWE rejecting 36.5% vs "
  f"~37%); and triple median {VT} per separation bin from 0.87 to 1.60 versus the published 0.87–1.62, all "
  "eight bins agreeing to ≤0.04. Triple P90 values run 10–14% above PS25 — a heavier contaminating tail, "
  "conservative for our purpose, traced to a simplified L(M) relation. All synthetic catalogues include the "
  f"intrinsic El-Badry velocity truncation ({VT} ≤ 2.23/√M<sub>tot</sub>) and 2D measurement noise "
  "resampled pair-by-pair from the real error columns.")

P("4. Injection–recovery: the pseudo-signal map", "h1")
P(f"We generate catalogues with truth {GM} ∈ {{1.0; 1.4 above 2 kau (EFE-saturated)}} × f(e) ∈ {{thermal; "
  "superthermal; s-dependent}} × residual triple fraction f<sub>trip</sub> ∈ {0, 0.1, 0.2, 0.3}, and "
  "evaluate three estimators, each given the <i>correct</i> per-zone f(e), the same noise model and the "
  "same selection — the only remaining degree of freedom is the treatment of the tail. E1 (full-median, "
  "Chae-like): squared ratio of full medians test/validation, corrected by the expected Newtonian ratio. "
  f"E2 (truncated, this work): medians truncated at {VT} &lt; √2 with the identical truncation applied to "
  f"zone-matched noisy templates, inverted on a {GM} grid. E3 (mixture, PS-like): Poisson likelihood of the "
  f"test-zone histogram as (1−f)·binaries<sub>{GM}</sub> + f·triple template, ({GM}, f) free.")
story.append(Image("figure_maitresse_G.png", width=16.4*cm, height=8.85*cm))
P("<b>Figure 1.</b> The pseudo-signal map: three estimator families confronted with the same synthetic "
  "universes (El-Badry selection + real Gaia noise). Top row: Newtonian truth; bottom row: MOND-EFE truth "
  f"({GM} = 1.4 above 2 kau). Columns: eccentricity distributions. On Newtonian truth the full-median "
  f"estimator (red) manufactures {GM} = 1.08–1.13 at realistic contamination, identically in all f(e) "
  "panels; the truncated (blue) and mixture (green) estimators stay at unity. On MOND truth all three "
  "recover ≥1.26. Single-seed grid; statistical jitter ≈ ±0.02.", "cap")
P(f"On Newtonian truth, E1 manufactures {GM} = 1.08–1.13 at f<sub>trip</sub> = 0.2–0.3, identically in all "
  "three f(e) panels — the culprit is the residual tail, not the eccentricity correction when the latter is "
  "correct; E2 stays ≤1.04 and E3 at 1.00. On MOND truth all three recover ≥1.26 (E2: 1.31–1.36, the ~5% "
  "dilution being a↔r<sub>p</sub> mixing at the 2 kau boundary). The two regimes never overlap under E2/E3: "
  "the controversy is decidable with DR3-class data. Two methodological by-products: (i) a Rice-type bias — "
  f"2D noise inflates |{VT}| more in the test zone where v<sub>c</sub> is small; before zone-matched noisy "
  "templates even E2 returned 1.18 on pure Newtonian truth; the same-cuts-on-simulations rule must extend "
  f"to noise, per zone; (ii) the {GM}↔f<sub>trip</sub> degeneracy of mixture fits as f<sub>trip</sub> → 0.")

P("5. DR3 measurement", "h1")
P(f"With the calibrated estimators: {GM}<sub>test</sub>(E2) = 1.045 [1.025, 1.068] (empirical per-zone e), "
  "stable under RUWE &lt; 1.2 (1.047) and the HRD variant (1.059); f(e) systematic band {empirical 1.045; "
  "thermal 0.958; superthermal 0.957} — the empirical-eccentricity correction alone pushes the estimate up "
  "by ~9%, and the full band brackets unity. The mixture likelihood on the real test zone gives "
  f"{GM} = 1.05 [1.04, 1.07] with f<sub>trip</sub> = 0.15, in striking external agreement with the "
  f"independent PS25 value of 0.17; {GM} = 1.4 is rejected at {DL}ln L = 129 (≈16{SG}), while "
  f"{DL}ln L({GM}=1) = 4.8 lies within the f(e) systematic band (Fig. 2).")
story.append(Image("fig2_mixture_fit.png", width=13.4*cm, height=8.2*cm))
P(f"<b>Figure 2.</b> Test-zone (2–30 kau) {VT} histogram of the fiducial Gaia DR3 sample versus the joint "
  f"({GM}, f<sub>trip</sub>) mixture fit and its components. The grey dash-dotted curve shows the same "
  f"model at the MOND-EFE prediction {GM} = 1.4, rejected at {DL}ln L = 129.", "cap")
P(f"The deep bin (N = 341) is perspective-sensitive: {GM}(E2) = 0.97 without and 1.11 ± 0.07 with the "
  f"angular correction, part of the shift being radial-velocity noise injection ({SG}<sub>RV</sub>·θ folded "
  "into the modulus where v<sub>c</sub> ~ 0.1 km s<super>−1</super>); with the thermal/superthermal e-models the "
  "corrected value drops to ~1.01. The full deep-bin bracket [≈0.97, 1.19] is consistent with Newton and "
  f"excludes 1.35–1.40 at ≈3{SG}; the global test zone excludes it at ~16{SG} (statistical), with "
  "systematics bounded by the published band (Fig. 3). Reading the map of Fig. 1 backwards, a full-median "
  f"analysis of this same sample with empirical-e correction would report {GM} ≈ 1.12 (global) to 1.20 "
  "(deep) — reproducing the amplitude of the claimed anomaly from a data set that robust estimation shows "
  "to be Newtonian.")
story.append(Image("fig3_gamma_vs_gN.png", width=12.6*cm, height=8.15*cm))
P(f"<b>Figure 3.</b> Recovered {GM} versus internal Newtonian acceleration (E2, perspective-corrected, "
  "fiducial cuts). Circles: empirical per-zone f(e) with 68% bootstrap intervals; open squares: thermal "
  "f(e) (lower edge of the systematic band). The MOND-AQUAL + EFE prediction (shaded) is excluded across "
  "the full acceleration range.", "cap")

P("6. Pre-registered protocol for Gaia DR4", "h1")
P("DR4 (2026) will multiply proper-motion precision by 2–4 and the usable sample by 3–5, driving "
  f"{SG}<sub>stat</sub>({GM}) below 0.005: estimator systematics will fully dominate, and post-hoc analysis "
  "choices will be able to produce any verdict. We therefore freeze, before the release, a complete "
  "protocol — sample and cuts; the three zones; estimators E2 and E3 exactly as calibrated in Sect. 4 with "
  "their measured biases; the mandatory f(e) band; noise-matched, selection-matched templates — and four "
  "decision criteria written in advance, including a STOP-and-publish rule if the validation-zone control "
  "fails and an explicit grey-zone outcome. The protocol document and code hash will be deposited on Zenodo "
  "prior to DR4; the present DR3 measurement serves as its calibration standard.")

P("7. Discussion", "h1")
P("The map quantifies both prior qualitative claims: the statement of Banik et al. (2026) that a loose "
  f"{VT} limit can create a MOND-like signal, and the PS25 criticism of small-separation triple "
  "calibration, merge into a single measured number — the leverage of the residual triple tail on a full "
  f"median, ≈ +0.08–0.13 in {GM} at realistic contamination, insensitive to the eccentricity model when the "
  "latter is correct, and additive with e-model mismatch when it is not. Combined with the deprojection "
  "sensitivity demonstrated on the 3D sample by Saad & Ting (2026), the Chae–Banik disagreement is now "
  "localized at every front where it appears. Two caveats temper the DR3 numbers rather than the "
  "conclusion: the per-pair eccentricities of Chae (2024) are prior-dependent inferences, and their use "
  "defines the upper edge of our systematic band; and the deep bin couples perspective correction, "
  "radial-velocity availability and small numbers — a warning for any analysis that leans on the widest "
  "separations. Our triple model, while anchored on PS25 oracles, is one family; the injection–recovery "
  "framework is precisely the tool that lets any alternative triple model be substituted and the map "
  "redrawn.")

P("8. Conclusions", "h1")
P("Confronted with identical synthetic universes, the estimator families in use in the literature diverge "
  "exactly as the literature does; the divergence is manufactured by the residual-triple tail acting on "
  f"full medians, not by gravity. Robust estimators calibrated on that map measure {GM} = 1.00–1.05 on "
  "Gaia DR3, exclude the MOND-EFE prediction at high significance, and — being frozen and pre-registered — "
  "turn Gaia DR4 into a decisive, tamper-proof test.")

P("Data and code availability", "h1")
P("The full pipeline (population engine, injection–recovery grid, estimators), the reduced catalogue, all "
  "validation outputs and the pre-registration document are assembled in a complete reproducibility "
  "package (scripts, data, per-phase logs) released with this preprint. Papers II and III of this series "
  "apply the measurement of this paper to, respectively, the elasticity of the dark sector and the "
  "coherence structure of superfluid dark matter.")

P("References", "h1")
for r in [
 "Banik I., Zhao H., 2018, MNRAS, 480, 2660",
 "Banik I., Pittordis C., Sutherland W., et al., 2024, MNRAS, 527, 4573",
 "Banik I., et al., 2026, arXiv:2602.24035 (Quality Framework)",
 "Chae K.-H., 2023, ApJ, 952, 128",
 "Chae K.-H., 2024, ApJ, 960, 114",
 "Chae K.-H., 2025, arXiv:2502.09373",
 "Chae K.-H., Lee B.-C., Hernandez X., et al., 2026, arXiv:2601.21728",
 "Chae K.-H., Yoon Y., 2026, arXiv:2607.14450 (data quality control and multiple-star modeling)",
 "El-Badry K., 2019, MNRAS, 482, 5018",
 "El-Badry K., Rix H.-W., Heintz T. M., 2021, MNRAS, 506, 2269",
 "Hernandez X., et al., 2024, MNRAS, 533, 729",
 "Hwang H.-C., Ting Y.-S., Zakamska N. L., 2022, MNRAS, 512, 3383",
 "Offner S. S. R., Moe M., Kratter K. M., et al., 2023, ASPC, 534, 275",
 "Pittordis C., Sutherland W., 2018, MNRAS, 480, 1778",
 "Pittordis C., Sutherland W., 2023, OJAp, 6, 4",
 "Pasquini L., Saglia R., Patat F., et al., 2026, arXiv:2602.04661 (wide binaries without viable bound Newtonian orbits)",
 "Pittordis C., Sutherland W., Shepherd P., 2025, arXiv:2504.07569 (PS25)",
 "Tokovinin A., 2014, AJ, 147, 87",
 "Tokovinin A., Kiyaeva O., 2016, MNRAS, 456, 2070",
 "Saad S. M., Ting Y.-S., 2026, arXiv:2603.11015 (orbital-modeling sensitivity of the 3D sample)",
]:
    P(r, "ref")

def footer(canvas, doc):
    canvas.saveState()
    canvas.setFont("DejaVu", 8)
    canvas.setFillColor(colors.grey)
    canvas.drawCentredString(A4[0]/2, 1.1*cm, f"{doc.page}")
    canvas.drawString(1.9*cm, 1.1*cm, "Boufourou 2026 — preprint")
    canvas.restoreState()

doc = SimpleDocTemplate("Boufourou_2026_PaperI_wide_binaries.pdf", pagesize=A4,
                        leftMargin=1.9*cm, rightMargin=1.9*cm,
                        topMargin=1.8*cm, bottomMargin=1.8*cm,
                        title="Estimator forensics for the wide-binary gravity test",
                        author="Hicham Boufourou")
doc.build(story, onFirstPage=footer, onLaterPages=footer)
print("PDF construit")
