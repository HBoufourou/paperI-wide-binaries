# Joint trajectory and photocentre observer — declared before new libraries

8 September 2026. This is new conditional method development. Earlier releases and data remain immutable. No actual DR4 data, human blinding or astrophysical population completeness is claimed.

## Operator

The fixed 34/66-month along-scan operator from dr4_binary is extended in a separate file. A relative outer trajectory, evaluated at every shared observing epoch, is projected into the same tangent plane as the inner photocentre offsets. Source displacements are −m2/M times the relative displacement and +m1/M times it. Both displacements, the internal offsets and the declared correlated noise are added BEFORE the five/seven-parameter GLS fits. The source masses are true barycentric masses; photocentre luminosity follows the explicitly stated L∝m^4 relation.

The full wide motion is supplied as a trajectory. Its tangent velocity is not also added as a straight slope. Separate fitted wide, inner and noise contributions are returned as simulation truth for validation, never as decision features. Inner periods use the total inner mass. Additive source offsets and systemic proper motion lie in the fitted design space; distance, tangent plane and parallax factors remain known, and perspective/systemic-RV effects are omitted in this conditional experiment.

The old cadence, temporal/pair noise, source formal-error scaling, persistent slope floors and shared-release covariance are preserved. Full cross-release covariance is used for proper-motion differences. Its exact identity in this nested design is not asserted for real Gaia releases. Noise scaling matches mean component variance, not every measured marginal covariance.

## Checks before production

Implementation-only seeds 269000001–269000004. Zero outer trajectories must recover every legacy relative output to 1e−11 absolute/relative, with identical random draws. Artificial straight outer trajectories must recover their analytically prescribed source and relative slopes, leave residual/acceleration diagnostics invariant to 1e−8, and agree with the old straight-slope observer for relative outputs. For a curved outer trajectory, fitted PM must equal outer+inner+noise to 1e−10 and source2−source1. The combined residuals must be computed from the sum, not combined after fits as independent flags. Shared covariance matrices must remain unchanged by deterministic signals and positive definite. Invalid trajectory shape/nonfinite values and invalid mass fractions must be rejected.

## Trajectories and scope

Initial states will be drawn from an explicitly energy-defined phase-space parent and an instantaneous observational selection. That selected sample is not itself claimed stationary. Epoch propagation must use the gradient of the SAME scalar potential used for the initial energies. Newton and QUMOND use the same units and prescribed true masses. Integrator refinements, energy/angular-momentum conservation over diagnostic trajectories and absence of out-of-domain evaluations are required before scientific libraries. No square-root-of-force velocity prescription is permitted.

The finite-mass QUMOND potential and its validation are owned by the physics protocol; this document does not pre-certify them. A failed potential or trajectory check stops its dependent scientific use. Full decision libraries, population hypotheses, selection, seed allocation and thresholds will be frozen separately before their generation or evaluation. Current code fixtures are not scientific decision evidence.
