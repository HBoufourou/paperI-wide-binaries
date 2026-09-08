# Before-production amendment: condition on photometric total mass

8 September 2026. This amendment precedes every physical validation population, every scientific bank and every decision calibration/test. The potential grid and software-only fixtures were available; no gravity-decision result was inspected.

The first population draft treated the empirical photometric total as a prior for simulated true mass, then reduced simulated photometric mass in companion systems. That choice could give the decision engine an avoidable difference in photometric-mass distributions, while its formal error covariates stayed anchored to the original stars.

The executed design instead conditions pure and companion comparisons on the same empirical total photometric mass. For prescribed outer mass fractions and a drawn inner host/mass ratio, the exact L∝m^4 mapping gives g=(1+qin^4)^(1/4)/(1+qin), b=1−fraction_host*(1−g), and Mtrue_total=Mphot_total_anchor/b. Draw and retain host/qin before computing rM, velocities, selection and hierarchy. Do not draw them again in the inner-orbit helper.

The observed total mass now agrees exactly with its anchor while the true orbital mass includes the companion. Per-source component masses and error correlations remain conditional approximations; this change does not validate a Gaia mass/selection model. It changes a prospective population specification, not a failed result or a decision threshold. The original pre-amendment software fixtures were not scientific populations; the updated algebra and selection fixtures are rerun before production. Keep the anchor and b in every latent library for independent verification.
