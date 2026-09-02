# Amendment to the Gaia DR4 pre-registration

Date: 2026-09-02, before access to Gaia DR4 data.

The original `PREREGISTRATION_DR4.md` remains unchanged as a frozen historical
record. Its estimator definitions and decision rules are not rewritten here.
This amendment records a pre-DR4 provenance and implementation correction.

1. The DR3 calibration paragraph in Section 7 of the original document is
   invalid. Those numbers came from Chae's synthetic Newtonian v5 catalogue,
   not the real Gaia DR3 catalogue. They must not be used as a DR4 benchmark.
2. The verified input rule is strengthened: every catalogue must have a public
   source URL/DOI, filename, digest, schema, and row count checked before any
   inference. A `Newton_*` file is a simulation by definition in this data
   family and cannot enter an observational branch.
3. Trial gravity parameters must be applied before measurement noise and
   catalogue/estimator selection. For hierarchical triples, gamma scales the
   low-acceleration outer-orbit velocity only; compact inner photocentre motion
   remains Newtonian. Common random numbers are used across the gamma grid.
4. Radial-velocity uncertainty is propagated through the angular perspective
   correction.
5. Applying rule D1 to corrected Gaia DR3 gives validation offsets of +10.99%
   (empirical eccentricities; fail), +1.77% (thermal; pass), and +4.06%
   (superthermal; fail). The corrected DR3 exercise is therefore not a new
   calibration standard and is reported as model-dependent/non-conclusive.

No DR4 outcome has been inspected. The purpose of this amendment is to prevent
a known data-classification error and known implementation inconsistencies from
being frozen into the future analysis.
