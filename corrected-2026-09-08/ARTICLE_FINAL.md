# Conditional gravity decisions for wide binaries: orbital populations, shared astrometry and calibration-library transfer

H. Boufourou

## Abstract

We construct a conditional Newtonian/QUMOND decision experiment linking a conservative two-source potential, an invariant-defined orbital population, unresolved companions and shared astrometric fits over simulated 34- and 66-month baselines. Five information variants and an envelope over three calibration libraries are evaluated on identical catalogue draws. Two fresh test-library realizations provide 540 scenario-specific cells, each containing 1,000 catalogues of 1,000 systems. The full later-baseline envelope meets the observed classification targets in all 540 cells, but the stronger marginal Wilson classification and true-model-retention targets hold in 531 and 501 cells, respectively; minimum true-model retention is 98.1 per cent. Removing diagnostics satisfies the classification target in all 540 cells under the stronger criterion, showing that their inclusion is not a uniform improvement in this finite representation. Covariates alone meet none of the classification targets, while removing their direct marginal likelihood contribution exposes additional weak cells. Enlarging the calibration envelope recovers 126 correct catalogue-view decisions but also creates 15 wrong exclusive decisions from previous double rejections. These controlled results identify information and calibration trade-offs within a restricted simulator. They support reproducible conditional experimentation, without an observational gravity verdict or a guarantee of decisive Gaia DR4 inference.

Keywords: methods: statistical - binaries: visual - astrometry - gravitation

## 1. Introduction

Wide binaries offer a setting in which predictions for low-acceleration dynamics meet uncertain orbital populations, stellar masses, unresolved companions and astrometric selection. These ingredients enter the same observables. A high classification rate under a chosen simulator is therefore insufficient unless the contribution of its information blocks and the calibration of its decisions are separately examined. We develop and test a conditional inference chain addressing that methodological question. The hypotheses are Newtonian dynamics and one specified quasi-linear formulation of modified Newtonian dynamics, QUMOND. No inference about which hypothesis describes the actual sky is made here.

Unresolved companions and kinematic truncation are established concerns in this field. Clarke demonstrated that inner photocentre motion can produce a high relative-velocity tail under Newtonian dynamics. Banik et al., Chae, and Pittordis, Sutherland & Shepherd developed simulations, population comparisons and contamination models for wide-binary gravity tests. Penoyre, Belokurov & Evans studied astrometric fits and their changes between observation baselines. These precedents motivate a comparison on identical simulated parents; they also prevent treating any one of these ingredients as a newly discovered mechanism. [Clarke (2020)](https://academic.oup.com/mnrasl/article/491/1/L72/5613964); [Banik et al. (2024)](https://arxiv.org/html/2311.03436v2); [Chae (2024a)](https://arxiv.org/html/2402.05720v4); [Pittordis et al. (2025)](https://arxiv.org/html/2504.07569v2); [Penoyre et al. (2022)](https://arxiv.org/html/2111.10380v3)

We link a conservative interaction potential for two finite prescribed sources, an invariant-defined stationary orbital parent, inner photocentre motion and shared astrometry to a decision that profiles nuisance parameters and can abstain. Component-removal experiments isolate the information supplied by diagnostics, the later simulated astrometric product and the marginal covariate distribution. A separate library-envelope experiment evaluates calibration transfer to fresh physical libraries. These experiments test the utility of this particular construction, without claiming the first simulator, the first simulation calibration, or a universal ranking of published methods.

Quality cuts and the detectability of inner companions also have dedicated observational and simulation studies. The present diagnostics are therefore evaluated as information within a specified observer, rather than as a complete contamination census or as an assumed guarantee supplied by a longer mission. [Cookson et al. (2026)](https://arxiv.org/html/2602.24035v1); [Manchanda et al. (2023)](https://arxiv.org/html/2210.07781v3)

The numerical finite-mass treatment builds on the QUMOND action and existing work on the external-field effect and the full two-mass problem. It is an implementation within a deliberately restricted population model, not a new theory of gravity. The restrictions are necessary to interpret the measured classification rates: changing the distribution function, mass ratios, external-field geometry or observing process defines another experiment. [Milgrom (2010)](https://arxiv.org/html/0911.5464); [Banik & Zhao (2018)](https://arxiv.org/html/1509.08457v3); [Pflamm-Altenburg (2025)](https://arxiv.org/html/2509.01493v1)

This manuscript substantially replaces an earlier version whose purported observational input was a Newtonian mock catalogue. The earlier observational estimates and gravity-exclusion claims are withdrawn. Appendix A records the provenance correction. The new experiment is explicitly synthetic and has independent generation roles; it does not reinterpret the mock as observations or use the new method to restore the invalidated conclusion.

## 2. Physical and observational model

### 2.1 Conservative two-source interaction

Use dimensionless units G = M_total = a_0 = 1. The prescribed source mass ratio is q in {0.1, 0.3, 1}, with m_1 = 1/(1+q), m_2 = q/(1+q), and source centres R_1 = -m_2 r and R_2 = m_1 r. Each source is a rigid Plummer distribution of scale b_i = 0.00125 sqrt(m_i). Newtonian and QUMOND calculations use the same prescribed sources. The external Newtonian acceleration is uniform, has magnitude a_0, and points along a synthetic celestial z-axis. It is not a reconstructed Galactic-centre direction.

For the QUMOND interpolation function and its action primitive, we use

$$ν(y) = (1 + sqrt(1 + 4/y))/2,     Q′(z) = ν(sqrt(z)).$$

Let u_e = -e_z be the external potential gradient, opposite to the external acceleration, u_i the Newtonian gradient of source i, and Q_ph(y) = Q(y²) - y². Define the mixed inclusion/exclusion term

$$ΔQ_ph = Q_ph(|u_e + u_1 + u_2|) - Q_ph(|u_e + u_1|) - Q_ph(|u_e + u_2|) + Q_ph(|u_e|).$$

The interaction energy and relative potential are

$$U = U_N - (1/8π) ∫ ΔQ_ph d³x,     Ψ = U/(m_1 m_2),     a_relative = -∇Ψ.$$

Here U_N is the exact interaction of the two Plummer sources, evaluated by a one-dimensional quadrature. The nonlinear function is applied after summing the Newtonian source gradients. The construction consequently differs from superposing two independent MOND solutions. Its zero is Ψ at infinity equal to zero, permitted by the nonzero external field. An isolated logarithmic MOND potential is outside the declared domain. [Milgrom (2010)](https://arxiv.org/html/0911.5464)

The production quadrature uses a stable mixed-Hessian evaluation for weak internal fields to limit cancellation in the four-term difference. Its leading spatial-tail correction is proportional to -m_1 m_2 [ν_e(1+K_e/3)-1]/R, where K_e is the logarithmic derivative of ν at external Newtonian field unity. This correction concerns the spatial integration boundary; it is not substituted for the angular dependence of the orbital external-field asymptote. Quadrature refinement, boundary doubling, independent prolate-coordinate energy calculations and direct source-force comparisons test the implementation. [Banik & Zhao (2018)](https://arxiv.org/html/1509.08457v3)

The final grid stores 1,683 energies: 33 radii over 0.03 <= r <= 30, 17 signed angular cosines and three q values. We interpolate the scalar A = -rΨ in log radius and signed angular cosine with bicubic polynomials. Acceleration is the analytic derivative of the same interpolant. The sign of the angular cosine is retained: unequal source masses do not justify reflection symmetrization. Intermediate q values and out-of-domain radii are refused rather than extrapolated.

An initial-grid development criterion of 0.5 per cent triggered the prescribed refinement. The retained grid has a maximum development force discrepancy of 0.03799 per cent. Its maximum discrepancy on 21 reserved positions is 0.01025 per cent for QUMOND and 0.0002183 per cent for Newtonian Plummer sources, below a 1 per cent criterion at each position. The maximum over three additional earlier fine references is 0.03059 per cent. These are pointwise comparisons, not a uniform error certificate for continuum QUMOND. An independent reconstruction verifies 49,152 outward Bernstein intervals exactly for the encoded polynomial; this bound has a different scope from physical quadrature accuracy.

### 2.2 Invariant-defined orbital parent

The static axisymmetric interpolant conserves the relative specific energy and axial angular momentum,

$$E = |v|²/2 + Ψ,     L_z = (r × v)_z.$$

At each q, the parent phase-space distribution is specified by

$$F(E,L_z,q) ∝ (-E)^(β+1/2) I_admissible(E,L_z),     0.1 <= -E <= 5.$$

The three q values have equal coefficients in this phase-space definition. They are proposed uniformly, but their normalized parent and selected fractions need not be equal. The dependence on the conserved quantities establishes stationarity within the static interpolated model. It does not establish that this distribution describes the observed eccentricity or anisotropy distribution of wide binaries.

Two-integral distribution functions have a long history in axisymmetric stellar dynamics. Here stationarity follows directly because the specified distribution and mask remain constant along a trajectory conserving E and L_z; no originality is claimed for that general principle. [Hunter & Qian (1993)](https://academic.oup.com/mnras/article/262/2/401/1161204)

The common bound C_max = 1.63417323970216 satisfies Ψ >= -C_max/r for both gravities and all three q values. The admissibility mask requires E < -C_max/30 and

$$r_lo = L_z² / [C_max + sqrt(C_max² + 2 E L_z²)] > 0.035.$$

Because E >= L_z²/(2r²) - C_max/r, an orbit conserving E and L_z within this mask cannot cross either tabulated radial boundary. Invalid discriminants and vanishing L_z are rejected. These restrictions retain a stationary parent but exclude some orbital families, especially low-|L_z| states. The common numerical barrier is not an astrophysically neutral prior.

We propose log radius uniformly, position and velocity directions isotropically, and the binding fraction u = (-E)/(-Ψ) from a Beta(5/2,3/2) distribution. The speed is sqrt[2(-Ψ)(1-u)]. The importance weight for β = 1 is r³(-Ψ)³; other β values multiply it by (-E)^(β-1). This follows from the phase-space measure and is not a prescription obtaining velocities from the square root of a local force. The invariant mask generally makes the retained directions anisotropic despite isotropic proposals.

The preceding development experiment used β = 1 and 1.6 for fitting and β = 1.3 as a reserved stress case. Its failures motivated the prospective expansion to β in {1,1.3,1.6} in this manuscript. New β = 1.15 and 1.45 tests remain outside that grid. Energy enters these population weights only; it is never a measured feature used to classify an individual simulated system.

### 2.3 Empirical covariates, masses and companions

The empirical covariate reservoir is derived from Chae's corrected-RUWE observed catalogue, Zenodo record 10986733. It supplies distances, sky positions, total photometric masses and astrometric error scales, not the new simulated orbital truth. Its inherited selection includes distance, chance-alignment probability, quality, proper-motion precision, radial-velocity availability and a transverse-velocity condition. We do not reconstruct an unselected Gaia parent or claim to have measured the future selection function. The exact reduced array, raw catalogue checksum and selection provider are supplied. [Chae (2024b)](https://zenodo.org/records/10986733)

The fixed reservoir contains 29,732 pairs. The inherited cuts require d<200 pc, R_chance<0.01, both RUWE values<1.4, the larger flat-frame normalized PM error<0.1, at least one usable radial velocity with a finite error, and flat-frame transverse speed divided by v_N at most 2.23/sqrt(M_phot/M_sun). They further require either 0.2≤s/kAU<30 or 0.03≤g_N/a_0<0.3, with g_N=GM_phot/s². These cuts select the empirical reservoir before generating new orbits. Its original measured velocity is not used as the new simulated velocity, but the selection it influenced remains a conditioning choice.

The total photometric mass proposed from this reservoir is fixed before assigning the simulated true mass. A pure system has photometric masses equal to its true source masses. For an inner companion, the host is chosen equiprobably and q_in is uniform over [0.02,1]. Under luminosity proportional to mass to the fourth power,

$$g = (1 + q_in⁴)^(1/4)/(1 + q_in),     M_true,total = M_phot,total / [1 - h(1-g)],$$

where h is the true outer-source mass fraction of the host. The photometric host mass is its true total mass multiplied by g. Thus hidden mass enters the orbital dynamics while the proposed total photometric mass remains anchored. The outer mass ratio is not silently changed after selecting the potential. Detailed correlations between the new component mass split and the empirical errors remain approximate; mass-luminosity uncertainty, extinction and photometric mass errors are not fitted.

At most one source contains one inner companion. Its log10 semimajor axis is normally distributed with standard deviation 1.5 and location log10(a_prior), truncated to 0.01--300 AU. Its eccentricity density is 0.4+1.2e for 0 <= e <= 1, with uniform mean anomaly and isotropic orientation. Reject semimajor axis and eccentricity until a(1+e) < r_lo r_M/20, where r_M = sqrt(GM_true,total/a_0). This is a conservative imposed hierarchy, not a theorem of triple stability. The internal orbit is Newtonian Keplerian and the external force uses the host's total monopole mass. Internal quadrupole, secular coupling, tides and a complete QUMOND three-body solution are absent.

An initially resolved inner pair, at projected separation at least one arcsecond, rejects the entire proposed parent. Its phase is not redrawn to manufacture an unresolved system. Later resolution fractions are recorded without applying a further selection. The outer sample requires the initial projected photocentre separation, including the inner offset, to lie between 2 and 30 kAU. This instantaneous selection does not preserve stationarity of the selected catalogue. Position-measurement error is omitted from this separation. No cut on the newly generated velocity or diagnostic score is imposed.

The physical inner-prior locations are 40 and 80 AU for templates; 60 and 160 AU are additional tests. These labels describe the law before truncation and hierarchy/resolution selection, not the median of the retained companions. The mixture parameter f below is the fraction of accompanied systems among the selected parents, not an unselected multiplicity fraction.

### 2.4 Joint astrometric observation

Outer trajectories are integrated using the gradient of the same scalar potential defining their initial energy. At every observing epoch, the relative displacement is split between the sources according to their true mass fractions, projected into a common sky tangent plane and added to the inner photocentre displacement and astrometric noise before fitting. Source offsets and systemic linear motion lie in the fitted design space. Perspective acceleration, systemic radial-velocity effects, distance uncertainty and a varying Galactic field are omitted.

The simulated 34- and 66-month products use nested observations of the same systems. Three fixed scenarios represent uniform cadence with white noise, clustered cadence with correlated noise, and sparse cadence with a persistent slope floor. Five- and seven-parameter generalized least-squares fits produce proper motions and source residual/acceleration diagnostics. The covariance of the cross-baseline proper-motion difference uses the shared observations. Its identity in this nested simulation is not asserted for real Gaia releases. Noise scaling matches an empirical average component variance, not every real source covariance.

For each source and baseline, let χ5² and χ7² be the whitened residual sums of squares for the five- and seven-parameter fits, divided by the source noise variance, and let d=n_epoch-5. The residual score uses the cube-root transformation

$$z_score = [(max(χ5²,10^-100)/d)^(1/3) - (1 - 2/(9d))] / sqrt(2/(9d)).$$

With barΦ the standard-normal upper tail, the residual and acceleration scores are

$$r_score = -ln[max(barΦ(z_score),10^-300)],     a_score = max(χ5²-χ7²,0)/2.$$

The proper-motion-change score is Δμᵀ C_Δ^-1 Δμ, for each source and for the relative motion. These transformations define diagnostic coordinates, not final hypothesis-rejection probabilities. No independence between their components is required by the calibrated decision, and these simulator scores are not public RUWE values.

The total fitted proper motion decomposes into the fitted outer, inner and noise contributions for software validation, while only the observed sum enters inference. The wide trajectory is not counted again as an additional linear velocity. The production observer passed 174 implementation checks; population and integrator fixtures passed 47 and 11 checks, respectively. A separate 42-criterion physical-chain validation tested trajectory refinement, energy and L_z conservation, shared covariance and the joint observation. Every new library also checks the full epoch-state invariants and fitted decomposition. The construction is related to earlier simulated astrometric fits across baselines, while serving a different conditional gravity-decision experiment here. [Penoyre et al. (2022)](https://arxiv.org/html/2111.10380v3); [Pittordis et al. (2025)](https://arxiv.org/html/2504.07569v2)

## 3. Statistical decision

### 3.1 Observable representation and nuisance family

For each system, define the Newtonian circular normalization from the total photometric mass and observed projected separation s,

$$v_N = sqrt(G M_phot / s),     r_M,phot = sqrt(G M_phot / a_0).$$

We use a_0 = 1.2 × 10^-10 m s^-2, GM_sun = 1.3271244 × 10^20 m³ s^-2 and AU = 149597870700 m. Proper motions are converted with 4.740470463533349 km s^-1 per arcsecond per year at one parsec. Both signed projected velocity components, parallel and perpendicular to the observed separation, are divided by v_N. The formal precision coordinate is the square root of half the trace of the relative proper-motion covariance, converted to the same dimensionless units.

The full later-baseline representation has four separation classes in s/r_M,phot, with boundaries 0.3,1,3; eight classes per signed velocity component, with boundaries -sqrt(2),-1,-0.5,0,0.5,1,sqrt(2); four diagnostic classes with boundaries 3,10,100; and two precision classes split at 0.05. The diagnostic is the maximum of the source residual/acceleration scores at both baselines and source/relative proper-motion-change Mahalanobis scores using their shared covariance. The resulting histogram has 2,048 joint bins. Overflow bins make the score computable but do not establish physical support outside the declared model.

Within each gravity, β and fitting component, normalize the population weights and record n_eff = 1/sum(w_i²). A symmetric Dirichlet total concentration of one gives p = (n_eff p_empirical + 1/B)/(n_eff+1), with B=2048. This prevents zero logarithms; it does not create new physical evidence. The projections below preserve this smoothing by aggregating the full probabilities.

The selected-parent mixture is

$$p_h,θ = (1-f) p_h,pure + f[(1-η)p_h,40 + η p_h,80],$$

with β in {1,1.3,1.6}, f in {0,0.075,0.15,0.225,0.30}, and η in {0,0.5,1}; only η=0 is retained when f=0. This gives 39 nuisance points per gravity. η is a mixture fraction of two selected-component laws. A physical prior 80 AU corresponds to η=1; a physical prior 60 or 160 AU is not assigned a true η merely because an implementation interface contains a placeholder.

### 3.2 Information-removal experiments

All variants use the same parent catalogue indices, fitting populations, nuisance grid and calibration rules. They are diagnostic variants of one inference construction, not complete implementations of other authors' published pipelines.

| Variant | Information retained | Intended interpretation |
| --- | --- | --- |
| Full66 | Complete 2,048-bin joint representation | Proposed information set |
| NoDiagnostics66 | Exact 512-bin projection removing diagnostics | Contribution of the diagnostic block |
| Full34 | Proper motions/covariance at 34 months and source diagnostics at 34 months only | Simulated earlier-baseline reference |
| CovariatesOnly66 | Eight bins of separation and formal precision only | Contribution from the assumed covariate distribution |
| Conditional66 | Full66 mixture normalized within separation/precision strata | Remove the direct marginal likelihood term |

Full34 never reads later-baseline proper motions, covariance, diagnostics or difference fields. It remains a simulated reference whose diagnostics come from known epoch fits; it is not a replication of all public Gaia DR3 products. Full66 versus Full34 therefore measures an added information package, including difference diagnostics, rather than isolating observation duration alone.

CovariatesOnly66 can discriminate if the assumed orbital parent and selection predict different separation or precision distributions. That is not evidence of leaked latent labels. It measures a source of information whose transport to the real sky depends on the parent model. Conditional66 forms the mixture first and then normalizes by its marginal stratum probability. With z denoting the eight separation/precision strata and y the remaining coordinates, its score is

$$D_h,θ,conditional = 2 sum_(z,y) n_z,y log[n_z,y / (n_z p_h,θ(y|z))].$$

This score is computed for every θ before minimization. Subtracting two separately minimized deviances would not be equivalent. Removing the marginal term does not equalize the distributions of n_z between hypotheses or make power independent of the parent. References are not generated conditional on every possible observed stratum-count vector. The experiment therefore tests a conditional likelihood score, without claiming an exactly conditional calibration at fixed n_z.

### 3.3 Scores, calibration and abstention

For the ordinary joint or projected representation, use multinomial deviance D_h,θ = 2 sum_b n_b log[n_b/(N p_h,θ,b)], with zero-count terms zero and N=1000. Profile separately under each gravity to obtain D_h = min_θ D_h,θ and R_h = D_h-D_other. No asymptotic chi-square interpretation is imposed. Best nuisance points, near-best alternatives and the full profiles remain available.

For each calibration-library realization j and each nuisance point, simulate 1,000 catalogues. For A in {D,R}, the upper-tail rank value is

$$p_h,j,A = max_θ { [1 + count(A_cal,h,j,θ >= A_observed,h)] / 1001 }.$$

Within each library, combine the two statistics by p_h,j = min[1,2 min(p_h,j,D,p_h,j,R)]. The maximum over nuisance points and Bonferroni factor two have their finite-family interpretation under a calibrating empirical law. The baseline uses j=0. The prospective library-envelope method is

$$p_h,envelope = max_(j=0,1,2) p_h,j.$$

The maximum is taken after combining statistics within each library. Taking different library maxima separately for D and R would define another rule. Reject a gravity only when its chosen p value is at most 0.01. Select Newton only when it is retained and QUMOND rejected, select QUMOND for the converse, and otherwise report indeterminate with both-retained or both-rejected reason. The catalogue interface flags invalid or unsupported inputs as unsupported and indeterminate; it produces no Newtonian or QUMOND selection for them.

Monte Carlo rank tests, nuisance maximization and the plus-one treatment of sampled tail probabilities have established statistical antecedents. Their validity depends on the simulated laws and exchangeability assumptions; assembling them does not confer a new guarantee for unmodelled physical populations. [Dufour (2006)](https://jeanmariedufour.research.mcgill.ca/Dufour_2006_JE_MCT.pdf); [Phipson & Smyth (2010)](https://gksmyth.github.io/pubs/PermPValuesPreprint.pdf)

The envelope can only decrease rejection of each individual hypothesis relative to j=0. It need not monotonically decrease exclusive wrong classifications: a double rejection can become retention of only the wrong hypothesis. Correct decisions, wrong decisions, true-hypothesis retention and both reasons for abstention must therefore be reported separately.

The rank construction is finite-sample valid under exchangeability for the corresponding empirical law; nuisance maximization and the Bonferroni combination can be conservative. Validity does not mean a rejection rate exactly equal to the nominal threshold and does not automatically transfer to a separately generated finite library. The envelope includes three empirical calibration realizations. It does not include every possible realization, continuous nuisance distribution or real Gaia observing process. Fresh test realizations evaluate this transfer without redefining it as a universal coverage guarantee.

## 4. Experimental design and verification

### 4.1 Development history and prospective test

The preceding experiment generated 22 physical libraries and found strong discrimination in 216 labelled cells, but only 208 met the stronger Wilson classification target. All eight exceptions involved Newton with β=1.3, outside its fitting grid. Minimum retention of the true model within the fitted family was 98.3 per cent, below a nominal 99 per cent level. Those observations motivated the present protocol. The old outputs are preserved as development evidence, not relabelled as new held-out tests.

The present fitting libraries are the same six fixed populations. The original calibration libraries define j=0; two fresh independent sets define j=1 and j=2. Each fresh calibration set has Newton/QUMOND times pure,comp40,comp80. Each of two fresh test sets has both gravities times pure,comp40,comp80,comp60,comp160. Thus 32 new physical libraries supplement 12 fixed historical inputs. Each new library starts with 20,000 accepted weighted parents; a prespecified effective-size rule permits 40,000 with the same seed prefix if needed before decision evaluation. The producer, protocol, inputs and all streams are hashed before scientific execution.

This design varies calibration and test realizations while conditioning on the fixed fitting libraries. It is not a full repetition of the entire fitting/calibration/test pipeline. Empirical covariates are reused across physical roles, while new physical states and astrometric noise use separate streams. Repeated catalogues sample these finite libraries; increasing the number of draws does not increase their physical support.

Calibration comprises three library realizations, two gravities, 39 nuisance points and 1,000 catalogues, each under three paired scenarios: 702,000 catalogue views scored by the five variants. In each fresh test realization, each gravity has 39 cases at the three fitting β values: one pure case per β and three positive contamination fractions crossed with four physical companion priors. Six additional cases per gravity use β=1.15/1.45 with pure systems or 30 per cent contamination under priors 60/160. There are 90 cases per realization,180 in total and 540 scenario-specific cells, with 1,000 catalogue draws per cell. The 540,000 views are paired across variants and calibration rules.

The test grid does not repeat every calibration nuisance point: f=0.15 and mixtures with η=0.5 have no dedicated transfer-test cells. Success in the tested cells therefore does not separately certify transfer at all 39 fitted nuisance points.

### 4.2 Predeclared targets and interpretation

For every cell, the observed classification target is correct decisions at least 95 per cent and wrong decisions at most 1 per cent. A separate true-model-retention target is at least 99 per cent. Stronger descriptive targets require the lower 95 per cent Wilson bound on correct decisions to reach 95 per cent, the upper bound on wrong decisions to be at most 1 per cent, and the lower bound on retention to reach 99 per cent. Classification and retention gates are never conflated. Each variant and calibration rule is reported in every tested cell.

These are marginal Monte Carlo intervals conditional on the realized libraries, without simultaneous correction over all cells or propagation over an unlimited population of fitting-library rebuilds. Selecting a worst cell does not turn its interval into a simultaneous test. Comparisons report paired discordant counts and differences in correct, wrong and indeterminate rates. The protocol selects the envelope prospectively, not by choosing whichever rule performs better in the final tests. No threshold, bin boundary or test population is changed in response to those outcomes.

The design distinguishes objectives, data-generating mechanisms, methods and performance measures, and retains Monte Carlo uncertainty in its reporting, following established guidance for simulation studies. [Morris et al. (2019)](https://onlinelibrary.wiley.com/doi/full/10.1002/sim.8086)

### 4.3 Computational verification

The upstream energy, interpolation, trajectories and observation operator have separate implementation and independent reconstruction controls. The earlier statistical reconstruction compared every observable histogram, weighted template, catalogue index, profiled score, calibrated probability and decision without importing the production inference engine. The new experiment additionally checks its fresh physical libraries, early-baseline isolation, exact projection variants, mixture-before-conditioning, the envelope rule and new seed allocation. Complete new verification results are reported below and recorded with hashes.

Independent here means a separate numerical implementation or a separate simulation stream, as specified for each check. The development and checks used AI assistance; they are not described as independent human peer review or as validation by an external astrophysics group. The released material enables such scientific scrutiny.

## 5. Results

### 5.1 Completed physical and computational experiment

All 32 fresh libraries met the effective-size rule at 20,000 accepted parents; no 40,000-parent regeneration was required. They contain 640,000 fresh systems and 1,920,000 paired scenario views. Together with the 12 historical fitting/calibration inputs, the inference experiment uses 44 physical libraries and 880,000 accepted parents. The minimum effective size across their five β weight sets is 16,320.7, and the maximum normalized parent weight is 0.000193509. These global diagnostics do not bound support in every observable bin or turn resampled catalogues into independent new physical systems.

The independent physical-library check passed 6,143 grouped checks covering all accepted initial rows and the 132 saved observation files. Its scope includes reconstructed invariants, masses, selection, weights and recorded observation algebra. Full regeneration of every rejected proposal, epoch trajectory, noise realization and astrometric fit was not part of this check.

The complete separate statistical reconstruction passed 28,842 grouped checks. It regenerated the measured codes, weighted templates, all 414 million saved parent indices, count hashes, profiles, rank values, decisions and descriptive tables without importing the production inference engine. Calibrated probabilities and decisions were exactly reproduced.

The public catalogue interface was also replayed on four saved catalogues in all five variants, both calibration rules and three scenarios:120 calls and 497 checks passed, with exactly matching probabilities and decisions. These checks establish numerical agreement within the stated scopes. They do not establish the empirical adequacy of the physical population.

### 5.2 Classification and true-model retention

Table 1 gives cell extrema and Table 2 keeps the classification and true-model-retention targets separate. Full66 with the prespecified envelope meets the observed classification target in all 540 cells: correct decisions are at least 98.1 per cent and wrong decisions at most 1.0 per cent. The stronger Wilson classification gate passes in 531 cells. The observed 99 per cent retention target passes in 537 cells, while its stronger Wilson gate passes in 501. Thus the observed classification success does not establish the intended retention reliability across the tested cells.

| Variant / calibration | Minimum correct (%) | Maximum wrong (%) | Minimum true retained (%) |
| --- | --- | --- | --- |
| Full66 / base | 97.5 | 0.9 | 97.5 |
| Full66 / envelope | 98.1 | 1.0 | 98.1 |
| NoDiagnostics66 / base | 98.5 | 0.1 | 98.5 |
| NoDiagnostics66 / envelope | 98.7 | 0.1 | 98.7 |
| Full34 / base | 92.6 | 0.9 | 98.6 |
| Full34 / envelope | 92.7 | 0.7 | 98.8 |
| CovariatesOnly66 / base | 0.0 | 2.2 | 97.7 |
| CovariatesOnly66 / envelope | 0.0 | 1.4 | 98.4 |
| Conditional66 / base | 90.1 | 2.5 | 97.3 |
| Conditional66 / envelope | 90.1 | 2.1 | 97.9 |

Table 1. Extrema over 540 cells for each fixed variant/rule. Columns may attain their extrema in different cells. These are conditional Monte Carlo summaries, not population-wide bounds. Every cell and its marginal interval is supplied with the data.

| Variant / calibration | Class. count | Class. Wilson | Ret. count | Ret. Wilson | Both Wilson |
| --- | --- | --- | --- | --- | --- |
| Full66 / base | 540 | 532 | 532 | 488 | 488 |
| Full66 / envelope | 540 | 531 | 537 | 501 | 501 |
| NoDiagnostics66 / base | 540 | 540 | 534 | 469 | 469 |
| NoDiagnostics66 / envelope | 540 | 540 | 538 | 493 | 493 |
| Full34 / base | 536 | 489 | 533 | 449 | 439 |
| Full34 / envelope | 535 | 495 | 538 | 485 | 468 |
| CovariatesOnly66 / base | 0 | 0 | 506 | 375 | 0 |
| CovariatesOnly66 / envelope | 0 | 0 | 528 | 436 | 0 |
| Conditional66 / base | 528 | 502 | 535 | 501 | 480 |
| Conditional66 / envelope | 528 | 508 | 537 | 514 | 492 |

Table 2. Numbers of 540 cells satisfying the criteria in section 4.2. “Class.” concerns correct≥95 per cent and wrong≤1 per cent; “Ret.” concerns true-model retention≥99 per cent. Wilson columns apply the corresponding marginal 95 per cent interval bounds. They are not simultaneous guarantees over all cells.

FIGURE:gates

Figure 1. Counts of tested cells meeting the stronger classification and retention criteria, plotted separately. Open circles use the original calibration library; filled squares use the prespecified three-library envelope. Catalogue draws are paired across methods and rules. The figure shows that additional information and additional calibration libraries do not improve every decision criterion uniformly.

The three Full66-envelope cells missing the observed retention target all occur under Newton with β=1 in test realization 2. The pure sparse-floor cell retains the true model in 986/1000 draws, with 10 wrong exclusive decisions and 4 double rejections. At f=0.30 with the 160-AU prior, the uniform scenario retains it in 987/1000 draws, while the sparse-floor scenario retains it in 981/1000, with 2 wrong decisions and 17 double rejections. The marginal retention interval for that last cell is[97.052,98.780] per cent. The pure example is inside the fitted physical family; the shortfall is not confined to outside-prior tests.

Under Full66-envelope, Wilson retention passes in 257/270 cells of realization 1 and 244/270 of realization 2; classification passes in 266/270 and 265/270, respectively. This difference is direct evidence of sensitivity within the two realized test libraries. Selecting a worst cell and quoting its marginal interval does not establish a simultaneous global calibration test or a distribution over arbitrary library rebuilds.

### 5.3 Information supplied by the later product and by diagnostics

With the envelope, Full66 has a higher correct-decision rate than Full34 in 424 cells, a lower rate in 9, and an unchanged rate in 107. The paired difference ranges from−0.4 to+7.1 percentage points. The mean difference over the equally weighted cells is+1.0128 points, with positive means in both test realizations. This supports the value of the later simulated information package within this experiment. It does not isolate baseline duration, reproduce public DR3-to-DR4 products, or establish a real-survey power forecast.

Diagnostic inclusion has a much less favourable interpretation. Relative to NoDiagnostics66 under the envelope, Full66 improves the correct rate in 122 cells and degrades it in 103; its change ranges from−1.5 to+1.0 percentage points. Across the 540,000 paired catalogue views,451 gain a correct decision and 440 lose one, a net gain of only 11. Wrong decisions increase by 148 and decrease by 2. The mean correct-rate difference has opposite signs between the two test realizations. NoDiagnostics66 satisfies the stronger classification gate in all 540 cells, compared with 531 for Full66, although its stronger retention count is 493 rather than 501.

The diagnostic block therefore has no demonstrated uniform advantage in this fixed representation and calibration scheme. This is not evidence that astrometric diagnostics contain no physical information. Adding the block changes the histogram dimension, support and calibrated score; the present experiment does not separately attribute the trade-off to each of those mechanisms. The predeclared Full66 result remains the primary reported construction, and the ablation is not retrospectively promoted to a tuned optimum.

### 5.4 Covariate information and conditional likelihood

CovariatesOnly66 meets the observed classification target in none of 540 cells under either rule. Separation and formal precision alone are therefore insufficient for the targeted discrimination in this experiment. Full66 improves the envelope correct rate relative to that control in every cell, by 42.8–100 percentage points.

The conditional-likelihood experiment asks a different question. Conditional66 removes the direct marginal separation/precision contribution after mixing and before profiling. Its envelope meets the observed classification target in 528 cells and the stronger classification gate in 508. The minimum correct rate is 90.1 per cent for QUMOND,β=1.6,f=0.30,prior 40 in the clustered scenario of realization 1. Its maximum wrong rate is 2.1 per cent for pure Newton,β=1, in the sparse-floor scenario of realization 2, with marginal 95 per cent interval[1.378,3.189] per cent. Both examples use fitted β values and fitted or pure component laws.

Full66 exceeds Conditional66 in correct rate in 318 cells and falls below it in 11, with a difference range of−0.3 to+9.6 percentage points. Thus the marginal information contributes to robustness in this experiment, even though the marginal control is insufficient alone. Neither comparison proves that a real orbital parent has been identified or that conditional scoring eliminates dependence on the population assumptions.

### 5.5 Calibration transfer and exclusive decisions

For Full66, the envelope raises correct rates in 61 cells and lowers them in none, recovering 126 correct catalogue-view decisions. It also introduces 15 wrong exclusive decisions, one in each of 15 cells, without removing an earlier wrong decision. All 15 transitions start from double rejection and end with only the wrong hypothesis retained. Maximum wrong rate consequently rises from 0.9 to 1.0 per cent; stronger classification passes fall from 532 to 531 even as stronger retention passes rise from 488 to 501.

This is consistent with the envelope's monotonic increase in individual p-values. It protects retention of each hypothesis but does not monotonically protect an exclusive two-hypothesis classification. The effect also depends on representation: Full34 loses 850 correct decisions net while removing 100 wrong decisions net; Conditional66 loses 41 correct decisions net while removing 57 wrong decisions. All individual retention changes are nonnegative. Counts refer to paired catalogue views with reused physical libraries and cadences, not independent astrophysical events.

The enlarged calibration family therefore improves some transfer criteria, with residual failures and representation-dependent costs. It does not provide a uniformly better classifier or an established 99 per cent retention guarantee on fresh physical laws. All 5,400 cell/method/rule rows and 7,020 paired comparison rows, including both reasons for indeterminacy and every failed gate, remain available.


## 6. Discussion

The experiment identifies a conditional information gain from the later astrometric product and a nonuniform trade-off from diagnostic inclusion. It also separates two calibration questions that can otherwise be conflated: whether a true hypothesis is retained and whether an exclusive decision is correct. The three-library envelope improves the former in the measured comparisons, while sometimes converting an abstention into an error. These are results of the specified inference chain, not a ranking of other published pipelines.

The component tests make the role of the parent distribution visible. A marginal-only control cannot meet the classification objective, while discarding the direct marginal term weakens some full-model decisions. Kinematic and diagnostic distributions still depend on the stipulated parent within each stratum. Extending flexibility in that parent, fitting library size, binning or companion law would constitute a new development experiment requiring new reserved tests. The present work retains its frozen choices and reports their weaknesses.

For future Gaia use, the immediate utility is to test proposed inference choices against a known simulated truth, preserve abstention and quantify the consequences of adding observations or nuisance freedom. A decisive observational analysis additionally requires an empirically constrained contamination and selection model, a justified link to the actual astrometric products, and verified error behaviour when those models vary. Those requirements are not established by the present synthetic success rates.


The physical assumptions delimit every rate in this study. The three true mass ratios, prescribed Plummer softening, fixed external-field magnitude and synthetic orientation do not span an astrophysical Galaxy. The E/L_z mask deliberately removes some orbit families. A flexible real population might absorb gravity differences that remain clear under this restricted parent. Even the conditional likelihood score retains model-dependent velocity and diagnostic distributions within the conditioning strata.

Contamination is represented by one inner companion at most, selected from specified distributions and observed through a known noise operator. Unbound associations, chance alignments, general higher multiplicity, realistic photometric and distance errors, time-dependent resolution/selection, Galactic tides and perspective effects require additional treatment. The internal Kepler-plus-external-monopole construction is not a full QUMOND triple integration, and its approximation error is not uniformly bounded here. A favourable test under these assumptions cannot certify all physical sources of a false verdict.

The astrometric diagnostics are defined by the simulator. Their usefulness must be checked against real source-level information and the actual product definitions before deployment. Likewise, an exact shared covariance in a nested synthetic observing design does not establish the covariance between real Gaia releases. The additional information in the simulated 66-month product should not be converted into a promise that DR4 will settle the gravity question.

The resulting experiment makes three properties separately inspectable: discrimination within a defined physical family, the information accounting for that discrimination, and transfer of calibration across finite libraries. The saved inputs and decisions support independent reuse and criticism, including further tests of the distribution function, contamination and selection. Computational agreement is evidence about implementation; empirical adequacy of the model remains a separate question.

## 7. Conclusions

A conservative two-source interaction, restricted invariant-defined orbital parent, shared astrometric observer and simulation-calibrated decision have been linked in a reproducible conditional experiment. The five fixed information variants use identical parent draws, and two fresh test-library realizations expose calibration transfer beyond the empirical laws used for calibration.

The later full product meets the observed classification target in all 540 tested cells under the envelope, but meets both stronger classification and retention targets in 501. Removing diagnostics can improve classification criteria, while removing the marginal covariate term weakens some cells. Additional calibration libraries improve true-model retention but can introduce wrong exclusive decisions from previous double rejections. No variant/rule satisfies all combined stronger targets in every tested cell.

The supported contribution is a quantified account of information and calibration trade-offs in this specified wide-binary simulator, together with the material needed to reproduce and challenge it. The invalid earlier observational interpretation is withdrawn. The experiment supplies neither a new measurement of gravity on the sky nor a promise that Gaia DR4 will decide between Newtonian dynamics and all forms of MOND.


## Data and software availability

The accompanying reproducibility dossier contains the physical model, input covariate array, source provenance, fixed and fresh libraries, protocols, statistical sources, saved catalogue indices, profiles, calibration references, decisions, all cell-level results and independent-verification products. A file-level SHA-256 manifest and relocation checks identify the released execution. The archive is provided with this replacement manuscript; no public deposit of the new archive is claimed before it actually occurs. The original repository and archive remain historical sources and must not be mistaken for the new corrected implementation. [Chae (2024b)](https://zenodo.org/records/10986733)

Python 3.12.14 and NumPy 2.3.5 were used for numerical calculations; pandas 3.0.1 is an import dependency of the inherited covariate-selection provider. Rendering dependencies are separated from scientific execution. Commands for inspecting saved results and for regenerating the experiment in a clean copy are included. Historical caches, outputs and failed tests are preserved rather than overwritten to manufacture a successful reproduction.

## Acknowledgements and author responsibility

The author thanks Austin Mander for identifying the original data-provenance problem. OpenAI Codex assisted with code development, numerical checks, literature retrieval, analysis and manuscript preparation. Separate AI-assisted implementations were used for computational cross-checking; this does not substitute for external scientific peer review. The exact tool/version metadata available to the author and its stated limitations are documented in the dossier. The author is responsible for reviewing and approving the scientific content and its claims before submission. This document does not assert that an external human specialist has already completed that review. [OpenAI (2026)](https://openai.com/codex/)

## Appendix A. Correction of the earlier observational claim

The earlier manuscript, “Estimator forensics for the wide-binary gravity test: the eccentricity-triple coupling manufactures a pseudo-signal, and a pre-registered protocol for Gaia DR4”, used an 81,088-pair input as if it were an observed Gaia catalogue. The provenance investigation identified it with the public Newtonian realization numbered 5 at Zenodo 10652994: all 36 non-identifier columns agree for the 81,088 identifier-matched pairs at relative tolerance 10^-9 and zero absolute tolerance, with matching missing values. This is numerical catalogue identity, not byte identity of differently formatted CSV files. Identifiers, positions and other preserved catalogue fields can remain observational when a generator replaces orbital-distance and proper-motion columns and adds hidden-companion masses. Matching preserved identifiers alone therefore does not establish observational provenance. [Chae (2024c)](https://zenodo.org/records/10652994)

The numerical estimates from that mock cannot support the previous measurement-on-the-sky interpretation or the reported gravity exclusion. Those conclusions are withdrawn. The corrected observed catalogue at Zenodo 10986733 has 81,880 pairs and a corrected RUWE column; its file has MD5 1b6c5063163a4e6c07043d13aeb70f55 and SHA-256 f8de60c31865beddcf6b3be0f5d1a3a8d8b705f1a3df9318bce21d7e46dcb4f7. Here it supplies explicitly conditional empirical covariates. No new observational estimate is inferred from the present simulation classification rates. [Chae (2024b)](https://zenodo.org/records/10986733)

The correction and development are distinct logical steps. Identifying the mock invalidates the old observational interpretation regardless of the new method's performance. The present physical and statistical construction is a substantial replacement with different claims, rather than a changed filename followed by the same conclusion. The accompanying editorial note explains this change of scope and leaves the handling of the ongoing submission to the editor.

## Appendix B. Reproduction boundaries

The physical interaction API rejects unsupported radii and mass ratios and requires matching model/evaluator validation hashes. The statistical execution uses the declared catalogue size, observation scenarios and calibrated family. Execution checks reject missing provenance or changed hashes. The measured-catalogue interface returns unsupported and indeterminate for invalid fields, covariance, size, scenario or unconfirmed scope. Its scope assertion is a caller declaration, not a test establishing that an external catalogue obeys the model. The interface and templates are supplied for reproducible conditional use, not as an already validated Gaia ingestion service. Agreeing code paths do not override these physical limitations.

Saved parent indices allow every count vector to be reconstructed without keeping repeated copies of the same catalogue for each variant. Canonical count hashes, score surfaces and decisions support independent numerical replay. Physical-library checks reconstruct accepted-row identities and invariant quantities; they do not replay every rejected proposal or every original along-scan fit. Each verification report states which calculations were repeated and which were checked through saved inputs and hashes.

The frozen source and protocol predate the current reserved tests, but follow earlier development results. They should be described as a prospective test within an iterative project, not as a publicly registered or human-blinded experiment. No earlier test is relabelled as fresh. New seeds, unchanged historical hashes and the disclosed development sequence make this boundary inspectable.

## References

Indranil Banik, Hongsheng Zhao, 2018, [The External Field Dominated Solution In QUMOND & AQUAL: Application To Tidal Streams](https://arxiv.org/html/1509.08457v3). ScieFed Journal of Astrophysics, 1, 1000008.

Indranil Banik, Charalambos Pittordis, Will Sutherland, Benoit Famaey, Rodrigo Ibata, Steffen Mieske, Hongsheng Zhao, 2024, [Strong constraints on the gravitational law from Gaia DR3 wide binaries](https://doi.org/10.1093/mnras/stad3393). Monthly Notices of the Royal Astronomical Society, 527, 4573--4615.

Kyu-Hyun Chae, 2024a, [Measurements of the Low-acceleration Gravitational Anomaly from the Normalized Velocity Profile of Gaia Wide Binary Stars and Statistical Testing of Newtonian and Milgromian Theories](https://doi.org/10.3847/1538-4357/ad61e9). The Astrophysical Journal, 972, 186.

Kyu-Hyun Chae, 2024b, [Python scripts to test gravity with the dynamics of wide binary stars](https://doi.org/10.5281/zenodo.10986733). Zenodo; record 10986733, corrected-RUWE catalogue.

Kyu-Hyun Chae, 2024c, [Python scripts to test gravity with the dynamics of wide binary stars](https://doi.org/10.5281/zenodo.10652994). Zenodo; record 10652994, Newtonian realization 5.

C. J. Clarke, 2020, [The distribution of relative proper motions of wide binaries in Gaia DR2: MOND or multiplicity?](https://doi.org/10.1093/mnrasl/slz161). Monthly Notices of the Royal Astronomical Society: Letters, 491, L72--L75.

Stephen A. Cookson, Indranil Banik, Kareem El-Badry, Will Sutherland, Zephyr Penoyre, Charalambos Pittordis, Cathie J. Clarke, 2026, [A quality framework for testing gravity with wide binaries: no evidence for MOND](https://doi.org/10.1093/mnras/stag342). Monthly Notices of the Royal Astronomical Society, 547, stag342.

Jean-Marie Dufour, 2006, [Monte Carlo tests with nuisance parameters: A general approach to finite-sample inference and nonstandard asymptotics](https://doi.org/10.1016/j.jeconom.2005.06.007). Journal of Econometrics, 133, 443--477.

C. Hunter, Edward Qian, 1993, [Two-integral distribution functions for axisymmetric galaxies](https://doi.org/10.1093/mnras/262.2.401). Monthly Notices of the Royal Astronomical Society, 262, 401--428.

Dhruv Manchanda, Will Sutherland, Charalambos Pittordis, 2023, [Wide Binaries as a Modified Gravity test: prospects for detecting triple-system contamination](https://doi.org/10.21105/astro.2210.07781). The Open Journal of Astrophysics.

Mordehai Milgrom, 2010, [Quasi-linear formulation of MOND](https://doi.org/10.1111/j.1365-2966.2009.16184.x). Monthly Notices of the Royal Astronomical Society, 403, 886.

Tim P. Morris, Ian R. White, Michael J. Crowther, 2019, [Using simulation studies to evaluate statistical methods](https://doi.org/10.1002/sim.8086). Statistics in Medicine.

OpenAI, 2026, [Codex](https://openai.com/codex/). OpenAI; software used 7–8 September 2026; historical version identifiers not established.

Zephyr Penoyre, Vasily Belokurov, N. Wyn Evans, 2022, [Astrometric identification of nearby binary stars -- I. Predicted astrometric signals](https://doi.org/10.1093/mnras/stac959). Monthly Notices of the Royal Astronomical Society, 513, 2437--2456.

J. Pflamm-Altenburg, 2025, [Numerical solutions of the complete two-body system in QUMOND](https://doi.org/10.1051/0004-6361/202555656). Astronomy & Astrophysics, 703, A68.

Belinda Phipson, Gordon K. Smyth, 2010, [Permutation P-values Should Never Be Zero: Calculating Exact P-values When Permutations Are Randomly Drawn](https://doi.org/10.2202/1544-6115.1585). Statistical Applications in Genetics and Molecular Biology, 9.

Charalambos Pittordis, Will Sutherland, Paul Shepherd, 2025, [Wide Binaries from Gaia DR3: testing GR vs MOND with realistic triple modelling](https://doi.org/10.33232/001c.142887). The Open Journal of Astrophysics.
