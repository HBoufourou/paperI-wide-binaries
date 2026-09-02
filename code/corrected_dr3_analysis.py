#!/usr/bin/env python3
"""Corrected observational analysis of the real Gaia DR3 catalogue.

This supersedes the observational part of ``phase_H.py`` from repository
versions 1/2.  It uses Chae's 81,880-pair Gaia catalogue with corrected RUWE
(Zenodo 10986733), propagates radial-velocity uncertainty in the perspective
correction, reports the mandatory validation-zone control for every
eccentricity family, and varies both binary and triple outer orbits with gamma
in the mixture likelihood.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from data_source import CATALOG_MD5, CATALOG_ROWS, CATALOG_URL, REPO_ROOT, load_catalog
from moteur_population import (
    orbit_plane,
    orient_matrix,
    to_sky,
    triple_vtilde,
    triples,
)

K = 4.74047e-3
SQRT2 = np.sqrt(2.0)
GAMMA_GRID = np.linspace(0.7, 1.9, 121)
FRACTION_GRID = np.linspace(0.0, 0.5, 51)
A0 = 1.2e-10


def sky_basis(ra_deg, dec_deg):
    ra, dec = np.radians(ra_deg), np.radians(dec_deg)
    radial = np.stack(
        [np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)]
    )
    alpha = np.stack([-np.sin(ra), np.cos(ra), np.zeros_like(ra)])
    delta = np.stack(
        [-np.sin(dec) * np.cos(ra), -np.sin(dec) * np.sin(ra), np.cos(dec)]
    )
    return radial, alpha, delta


def prepare_observables(frame):
    d1 = frame["d1[pc]"].to_numpy()
    d2 = frame["d2[pc]"].to_numpy()
    distance = 0.5 * (d1 + d2)
    separation = frame["s[kau]"].to_numpy()
    total_mass = (frame["M1[Msun]"] + frame["M2[Msun]"]).to_numpy()
    gn_a0 = 5.930e-3 * total_mass / (separation * 1000.0) ** 2 / A0

    rv1 = frame["RV1[km/s]"].to_numpy()
    rv2 = frame["RV2[km/s]"].to_numpy()
    erv1 = frame["RV1_err[km/s]"].to_numpy()
    erv2 = frame["RV2_err[km/s]"].to_numpy()
    ok1 = np.isfinite(rv1) & np.isfinite(erv1) & (np.abs(rv1) < 500)
    ok2 = np.isfinite(rv2) & np.isfinite(erv2) & (np.abs(rv2) < 500)
    both = ok1 & ok2
    rv_system = np.where(
        both,
        0.5 * (rv1 + rv2),
        np.where(ok1, rv1, np.where(ok2, rv2, 0.0)),
    )
    erv_system = np.where(
        both,
        0.5 * np.hypot(erv1, erv2),
        np.where(ok1, erv1, np.where(ok2, erv2, 0.0)),
    )

    dra_deg = (frame["RA2[deg]"].to_numpy() - frame["RA1[deg]"].to_numpy() + 180) % 360 - 180
    ra_mean = frame["RA1[deg]"].to_numpy() + 0.5 * dra_deg
    dec_mean = 0.5 * (frame["DEC1[deg]"] + frame["DEC2[deg]"]).to_numpy()
    r_mean, a_mean, d_mean = sky_basis(ra_mean, dec_mean)
    mua_mean = 0.5 * (frame["mu1ra[mas/yr]"] + frame["mu2ra[mas/yr]"]).to_numpy()
    mud_mean = 0.5 * (frame["mu1dec[mas/yr]"] + frame["mu2dec[mas/yr]"]).to_numpy()
    v_system = rv_system * r_mean + K * distance * (
        mua_mean * a_mean + mud_mean * d_mean
    )

    _, a1, q1 = sky_basis(frame["RA1[deg]"].to_numpy(), frame["DEC1[deg]"].to_numpy())
    _, a2, q2 = sky_basis(frame["RA2[deg]"].to_numpy(), frame["DEC2[deg]"].to_numpy())
    mu1a = np.sum(v_system * a1, axis=0) / (K * distance)
    mu1d = np.sum(v_system * q1, axis=0) / (K * distance)
    mu2a = np.sum(v_system * a2, axis=0) / (K * distance)
    mu2d = np.sum(v_system * q2, axis=0) / (K * distance)

    raw_dra = (frame["mu2ra[mas/yr]"] - frame["mu1ra[mas/yr]"]).to_numpy()
    raw_ddec = (frame["mu2dec[mas/yr]"] - frame["mu1dec[mas/yr]"]).to_numpy()
    corrected_dra = raw_dra - (mu2a - mu1a)
    corrected_ddec = raw_ddec - (mu2d - mu1d)
    va = K * distance * corrected_dra
    vd = K * distance * corrected_ddec
    v2d_corrected = np.hypot(va, vd)
    v2d_flat = K * distance * np.hypot(raw_dra, raw_ddec)

    # Proper-motion errors plus the radial-velocity contribution to the angular
    # perspective correction.  The latter is sigma_RV times the change of sky
    # basis between the two component positions.
    sigma_mua = K * distance * np.hypot(
        frame["mu1ra_err[mas/yr]"], frame["mu2ra_err[mas/yr]"]
    ).to_numpy()
    sigma_mud = K * distance * np.hypot(
        frame["mu1dec_err[mas/yr]"], frame["mu2dec_err[mas/yr]"]
    ).to_numpy()
    rv_coeff_a = np.sum(r_mean * (a2 - a1), axis=0)
    rv_coeff_d = np.sum(r_mean * (q2 - q1), axis=0)
    sigma_va = np.hypot(sigma_mua, erv_system * rv_coeff_a)
    sigma_vd = np.hypot(sigma_mud, erv_system * rv_coeff_d)
    sigma_v2d = np.sqrt((va * sigma_va) ** 2 + (vd * sigma_vd) ** 2) / np.maximum(
        v2d_corrected, 1e-12
    )

    circular_velocity = 0.94179 * np.sqrt(total_mass / separation)
    vtilde = v2d_corrected / circular_velocity
    vtilde_flat = v2d_flat / circular_velocity
    sigma_vtilde = sigma_v2d / circular_velocity

    base = (
        np.isfinite(vtilde)
        & np.isfinite(sigma_vtilde)
        & (distance < 200)
        & (frame["R_chance"].to_numpy() < 0.01)
    )
    fiducial = (
        base
        & (frame["ruwe1"].to_numpy() < 1.4)
        & (frame["ruwe2"].to_numpy() < 1.4)
        & (sigma_vtilde < 0.10)
    )
    strict = (
        base
        & (frame["ruwe1"].to_numpy() < 1.2)
        & (frame["ruwe2"].to_numpy() < 1.2)
        & (sigma_vtilde < 0.10)
    )

    return {
        "frame": frame,
        "d1": d1,
        "d2": d2,
        "distance": distance,
        "separation": separation,
        "total_mass": total_mass,
        "gn_a0": gn_a0,
        "vtilde": vtilde,
        "vtilde_flat": vtilde_flat,
        "sigma_vtilde": sigma_vtilde,
        "has_rv": ok1 | ok2,
        "fiducial": fiducial,
        "strict": strict,
    }


def zone_masks(state, selection):
    separation = state["separation"]
    gn_a0 = state["gn_a0"]
    return (
        selection & (separation >= 0.2) & (separation < 2),
        selection & (separation >= 2) & (separation < 30),
        selection & (gn_a0 >= 0.03) & (gn_a0 < 0.3),
    )


def binaries_scalefree(n, eccentricities, rng):
    values = eccentricities[np.isfinite(eccentricities)]
    e = np.clip(rng.choice(values, n), 0, 0.995)
    anomaly = rng.uniform(0, 2 * np.pi, n)
    x, y, vx, vy, _ = orbit_plane(np.ones(n), e, anomaly, 1.0)
    rotation = orient_matrix(n, rng)
    sx, sy, _ = to_sky(x, y, rotation)
    ux, uy, _ = to_sky(vx, vy, rotation)
    return np.hypot(ux, uy) * np.sqrt(np.hypot(sx, sy))


class ZoneTemplate:
    def __init__(self, state, zone, e_model, n=400_000, seed=3):
        rng = np.random.default_rng(seed)
        frame = state["frame"]
        if e_model == "empirical":
            eccentricities = frame["e"].to_numpy()[zone]
        elif e_model == "thermal":
            eccentricities = np.sqrt(rng.uniform(0, 1, 200_000))
        elif e_model == "superthermal":
            eccentricities = rng.uniform(0, 1, 200_000) ** (1 / 2.3)
        else:
            raise ValueError(e_model)
        self.signal = binaries_scalefree(n, eccentricities, rng)
        indices = rng.choice(np.flatnonzero(zone), n)
        self.sigma = state["sigma_vtilde"][indices]
        self.mass = state["total_mass"][indices]
        self.noise_x = rng.normal(size=n)
        self.noise_y = rng.normal(size=n)

    def observed(self, gamma=1.0, threshold=None):
        values = np.hypot(
            np.sqrt(gamma) * self.signal + self.noise_x * self.sigma,
            self.noise_y * self.sigma,
        )
        values = values[values <= 2.23 / np.sqrt(self.mass)]
        if threshold is not None:
            values = values[values < threshold]
        return values


def e2_estimate(state, selection, e_model, *, zone="test", threshold=SQRT2, seed=50):
    validation, test, deep = zone_masks(state, selection)
    target = test if zone == "test" else deep
    template_validation = ZoneTemplate(state, validation, e_model, seed=3)
    template_target = ZoneTemplate(state, target, e_model, seed=4)
    denominator = np.median(template_validation.observed(threshold=threshold))
    prediction = np.array(
        [
            np.median(template_target.observed(gamma, threshold=threshold))
            for gamma in GAMMA_GRID
        ]
    ) / denominator
    v_validation = state["vtilde"][validation]
    v_target = state["vtilde"][target]
    observed_ratio = np.median(v_target[v_target < threshold]) / np.median(
        v_validation[v_validation < threshold]
    )
    estimate = np.interp(observed_ratio, prediction, GAMMA_GRID)

    rng = np.random.default_rng(seed)
    bootstrap = []
    for _ in range(400):
        sample_validation = rng.choice(v_validation, len(v_validation))
        sample_target = rng.choice(v_target, len(v_target))
        ratio = np.median(sample_target[sample_target < threshold]) / np.median(
            sample_validation[sample_validation < threshold]
        )
        bootstrap.append(np.interp(ratio, prediction, GAMMA_GRID))
    low, high = np.quantile(bootstrap, [0.16, 0.84])
    return {
        "gamma": float(estimate),
        "low68": float(low),
        "high68": float(high),
        "n": int(target.sum()),
        "observed_ratio": float(observed_ratio),
    }


def validation_control(state, selection, e_model):
    validation, _, _ = zone_masks(state, selection)
    template = ZoneTemplate(state, validation, e_model, seed=3)
    observed = np.median(state["vtilde"][validation][state["vtilde"][validation] < SQRT2])
    expected = np.median(template.observed(threshold=SQRT2))
    delta = observed / expected - 1
    return {
        "observed_median": float(observed),
        "template_median": float(expected),
        "fractional_delta": float(delta),
        "passes_2_percent": bool(abs(delta) < 0.02),
    }


def mixture_fit(state, selection, e_model, triple_sample):
    _, test, _ = zone_masks(state, selection)
    bins = np.linspace(0, 2.4, 49)
    counts, _ = np.histogram(state["vtilde"][test], bins)
    binary_template = ZoneTemplate(state, test, e_model, n=800_000, seed=9)

    triple_zone = (
        triple_sample["survive"]
        & (triple_sample["r_p_kau"] >= 2)
        & (triple_sample["r_p_kau"] < 30)
    )
    rng = np.random.default_rng(10)
    indices = rng.choice(np.flatnonzero(test), triple_zone.sum())
    sigma = state["sigma_vtilde"][indices]
    noise_x = rng.normal(size=triple_zone.sum())
    noise_y = rng.normal(size=triple_zone.sum())
    mass = triple_sample["Mapp"][triple_zone]

    likelihood = np.empty((len(GAMMA_GRID), len(FRACTION_GRID)))
    binary_histograms = []
    triple_histograms = []
    for i, gamma in enumerate(GAMMA_GRID):
        binary = np.histogram(binary_template.observed(gamma), bins)[0].astype(float)
        binary /= binary.sum()

        triple_signal = triple_vtilde(triple_sample, gamma)[triple_zone]
        triple_observed = np.hypot(
            triple_signal + noise_x * sigma,
            noise_y * sigma,
        )
        triple_observed = triple_observed[triple_observed <= 2.23 / np.sqrt(mass)]
        triple_hist = np.histogram(triple_observed, bins)[0].astype(float)
        triple_hist /= triple_hist.sum()

        models = counts.sum() * (
            (1 - FRACTION_GRID[:, None]) * binary
            + FRACTION_GRID[:, None] * triple_hist
        ) + 1e-9
        likelihood[i] = np.sum(counts * np.log(models) - models, axis=1)
        binary_histograms.append(binary)
        triple_histograms.append(triple_hist)

    maximum = np.unravel_index(np.argmax(likelihood), likelihood.shape)
    gamma = GAMMA_GRID[maximum[0]]
    fraction = FRACTION_GRID[maximum[1]]
    profile = likelihood.max(axis=1)
    interval = GAMMA_GRID[profile >= profile.max() - 0.5]
    best_binary = binary_histograms[maximum[0]]
    best_triple = triple_histograms[maximum[0]]
    best_model = counts.sum() * ((1 - fraction) * best_binary + fraction * best_triple)
    return {
        "gamma": float(gamma),
        "low68": float(interval.min()),
        "high68": float(interval.max()),
        "f_trip": float(fraction),
        "delta_loglike_gamma_1": float(
            profile.max() - profile[np.argmin(np.abs(GAMMA_GRID - 1.0))]
        ),
        "delta_loglike_gamma_1_4": float(
            profile.max() - profile[np.argmin(np.abs(GAMMA_GRID - 1.4))]
        ),
        "bins": bins,
        "counts": counts,
        "best_model": best_model,
    }


def hrd_diagnostic(state):
    frame = state["frame"]
    magnitude1 = frame["MagG1"].to_numpy() - 5 * np.log10(state["d1"] / 10) - frame[
        "A_G1[mag]"
    ].to_numpy()
    magnitude2 = frame["MagG2"].to_numpy() - 5 * np.log10(state["d2"] / 10) - frame[
        "A_G2[mag]"
    ].to_numpy()
    colour = np.concatenate([frame["bp_rp1"].to_numpy(), frame["bp_rp2"].to_numpy()])
    magnitude = np.concatenate([magnitude1, magnitude2])
    valid = np.isfinite(colour) & np.isfinite(magnitude) & (colour > 0.1) & (colour < 3.2)
    polynomial = np.polyfit(colour[valid], magnitude[valid], 5)
    lobster1 = magnitude1 - np.polyval(polynomial, frame["bp_rp1"].to_numpy())
    lobster2 = magnitude2 - np.polyval(polynomial, frame["bp_rp2"].to_numpy())
    selected = (np.abs(lobster1 - lobster2) <= 0.40) & (np.minimum(lobster1, lobster2) >= -0.75)
    return selected


def make_figures(summary, state, mixture_models, acceleration_results):
    figures = REPO_ROOT / "figures"
    figures.mkdir(exist_ok=True)

    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.1), sharey=True)
    for axis, e_model in zip(axes, ("thermal", "empirical")):
        fit = mixture_models[e_model]
        bins = fit["bins"]
        centres = 0.5 * (bins[:-1] + bins[1:])
        axis.step(bins[:-1], fit["counts"], where="post", color="black", lw=1.2, label="Gaia DR3")
        axis.plot(centres, fit["best_model"], color="#1e8449", lw=2, label="best mixture")
        axis.set_title(
            f"{e_model}: $\\gamma$={fit['gamma']:.2f}, $f_{{\\rm trip}}$={fit['f_trip']:.2f}"
        )
        axis.set_xlabel(r"$\tilde{v}$")
        axis.grid(alpha=0.25)
    axes[0].set_ylabel("pairs per bin")
    axes[0].legend(fontsize=8)
    fig.suptitle("Corrected Gaia DR3 mixture fits: eccentricity-model dependence")
    fig.tight_layout()
    fig.savefig(figures / "fig2_mixture_fit_corrected.png", dpi=180)
    plt.close(fig)

    labels = ["empirical", "thermal", "superthermal"]
    colours = ["#2471a3", "#7d3c98", "#1e8449"]
    x = np.array([np.sqrt(low * high) for low, high in acceleration_results["bins"]])
    fig, axis = plt.subplots(figsize=(6.8, 4.4))
    axis.axhspan(1.35, 1.40, color="#c0392b", alpha=0.18, label="MOND-AQUAL + EFE")
    axis.axhline(1.0, color="black", ls="--", lw=1.1, label="Newton")
    for label, colour in zip(labels, colours):
        values = acceleration_results[label]
        y = np.array([item["gamma"] for item in values])
        low = y - np.array([item["low68"] for item in values])
        high = np.array([item["high68"] for item in values]) - y
        axis.errorbar(x, y, yerr=[low, high], marker="o", capsize=3, color=colour, label=label)
    axis.set_xscale("log")
    axis.invert_xaxis()
    axis.set_xlabel(r"$g_{\rm N}/a_0$  (deep regime $\rightarrow$)")
    axis.set_ylabel(r"recovered $\gamma$")
    axis.set_ylim(0.82, 1.55)
    axis.grid(alpha=0.25, which="both")
    axis.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(figures / "fig3_gamma_vs_gN_corrected.png", dpi=180)
    plt.close(fig)


def run_analysis():
    frame = load_catalog()
    state = prepare_observables(frame)
    validation, test, deep = zone_masks(state, state["fiducial"])
    e_models = ("empirical", "thermal", "superthermal")

    perspective = {}
    for name, mask in (("validation", validation), ("test", test), ("deep", deep)):
        difference = state["vtilde"][mask] - state["vtilde_flat"][mask]
        perspective[name] = {
            "median_delta": float(np.median(difference)),
            "p90_abs_delta": float(np.quantile(np.abs(difference), 0.9)),
            "rv_availability": float(state["has_rv"][mask].mean()),
        }

    validation_results = {
        model: validation_control(state, state["fiducial"], model) for model in e_models
    }
    e2 = {
        zone: {
            model: e2_estimate(
                state,
                state["fiducial"],
                model,
                zone=zone,
                seed=100 + 10 * index + (zone == "deep"),
            )
            for index, model in enumerate(e_models)
        }
        for zone in ("test", "deep")
    }
    e2_strict = {
        model: e2_estimate(state, state["strict"], model, seed=200 + index)
        for index, model in enumerate(e_models)
    }

    thresholds = {}
    for index, model in enumerate(e_models):
        thresholds[model] = {}
        for threshold in (1.2, 1.4, SQRT2, 1.6, 1.8):
            thresholds[model][f"{threshold:.3f}"] = e2_estimate(
                state,
                state["fiducial"],
                model,
                threshold=threshold,
                seed=300 + index,
            )

    triple_rng = np.random.default_rng(77)
    triple_sample = triples(1_200_000, triple_rng, cutset="sample", d_max=200.0)
    mixture_models = {
        model: mixture_fit(state, state["fiducial"], model, triple_sample)
        for model in e_models
    }

    acceleration_bins = [(3, 10), (1, 3), (0.3, 1), (0.03, 0.3)]
    acceleration_results = {"bins": acceleration_bins}
    for index, model in enumerate(e_models):
        acceleration_results[model] = []
        for bin_index, (low, high) in enumerate(acceleration_bins):
            selection = state["fiducial"] & (state["gn_a0"] >= low) & (state["gn_a0"] < high)
            # Reuse the E2 machinery by making this bin the temporary deep zone.
            original = state["gn_a0"]
            scaled = np.full_like(original, np.nan)
            scaled[selection] = 0.1
            state["gn_a0"] = scaled
            result = e2_estimate(
                state,
                state["fiducial"],
                model,
                zone="deep",
                seed=400 + 20 * index + bin_index,
            )
            state["gn_a0"] = original
            acceleration_results[model].append(result)

    hrd = hrd_diagnostic(state)
    summary = {
        "status": "model-dependent / non-conclusive",
        "scope": "corrected real-catalogue analysis; no high-significance gravity verdict",
        "catalogue": {
            "url": CATALOG_URL,
            "md5": CATALOG_MD5,
            "rows": CATALOG_ROWS,
        },
        "counts": {
            "parent": int(len(frame)),
            "fiducial": int(state["fiducial"].sum()),
            "strict_ruwe_1_2": int(state["strict"].sum()),
            "validation": int(validation.sum()),
            "test": int(test.sum()),
            "deep": int(deep.sum()),
            "hrd_proxy_retained": int((state["strict"] & hrd).sum()),
        },
        "perspective": perspective,
        "validation_control": validation_results,
        "e2": e2,
        "e2_strict": e2_strict,
        "threshold_scan": thresholds,
        "e3": {
            model: {key: value for key, value in fit.items() if key not in {"bins", "counts", "best_model"}}
            for model, fit in mixture_models.items()
        },
        "acceleration_bins": {
            "bins": acceleration_bins,
            **{model: acceleration_results[model] for model in e_models},
        },
    }

    results = REPO_ROOT / "results"
    results.mkdir(exist_ok=True)
    output = results / "corrected_dr3_summary.json"
    output.write_text(json.dumps(summary, indent=2) + "\n")
    make_figures(summary, state, mixture_models, acceleration_results)
    print(json.dumps(summary, indent=2))
    print(f"\nwrote {output}")
    return summary


if __name__ == "__main__":
    run_analysis()
