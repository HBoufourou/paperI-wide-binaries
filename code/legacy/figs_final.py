#!/usr/bin/env python3
"""figs_final.py — Fig.2 (donnees vs fit de melange) et Fig.3 (gamma par bin g_N/a0)."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

exec(open("phase_H.py").read())          # reutilise tout le pipeline H (namespace complet)

# ---------- Fig.2 : histogramme zone test, donnees vs modele ----------
fig, ax = plt.subplots(figsize=(7.2, 4.4))
ctr = 0.5*(bins[1:] + bins[:-1])
hb_ml = np.histogram(tT.observed(g_ml), bins)[0].astype(float); hb_ml /= hb_ml.sum()
tot = hd.sum()*((1-f_ml)*hb_ml + f_ml*ht)
ax.step(bins[:-1], hd, where="post", color="k", lw=1.4, label="Gaia DR3 data (test zone, 2–30 kau)")
ax.plot(ctr, tot, color="#1e8449", lw=2, label=f"best fit: $\\gamma$={g_ml:.2f}, $f_{{\\rm trip}}$={f_ml:.2f}")
ax.plot(ctr, hd.sum()*(1-f_ml)*hb_ml, "--", color="#2471a3", lw=1.5, label="binaries component")
ax.plot(ctr, hd.sum()*f_ml*ht, ":", color="#c0392b", lw=1.8, label="triples component")
hb14 = np.histogram(tT.observed(1.4), bins)[0].astype(float); hb14 /= hb14.sum()
ax.plot(ctr, hd.sum()*((1-f_ml)*hb14 + f_ml*ht), "-.", color="gray", lw=1.3,
        label="same model, $\\gamma$=1.4 (MOND-EFE)")
ax.set_xlabel(r"$\tilde{v}$"); ax.set_ylabel("pairs per bin")
ax.set_xlim(0, 2.4); ax.legend(fontsize=8); ax.grid(alpha=0.25)
fig.tight_layout(); fig.savefig("fig2_mixture_fit.png", dpi=170); plt.close(fig)

# ---------- Fig.3 : gamma par bin g_N/a0 (E2, e empirique + bande thermique) ----------
def E2_mask(mZ, e_model):
    mV, _, _ = zone_masks(fid)
    tV = ZoneTemplate(mV, e_model, seed=3); tZ = ZoneTemplate(mZ, e_model, seed=4)
    denom = np.median(tV.observed(Tcut=SQ2))
    pred = np.array([np.median(tZ.observed(g, Tcut=SQ2)) for g in GGRID])/denom
    vV, vZ = vt_corr[mV], vt_corr[mZ]
    r_obs = np.median(vZ[vZ < SQ2])/np.median(vV[vV < SQ2])
    ghat = np.interp(r_obs, pred, GGRID)
    gb = []
    for _ in range(250):
        bV = rng.choice(vV, len(vV)); bZ = rng.choice(vZ, len(vZ))
        gb.append(np.interp(np.median(bZ[bZ < SQ2])/np.median(bV[bV < SQ2]), pred, GGRID))
    gb = np.array(gb)
    return ghat, np.quantile(gb, 0.16), np.quantile(gb, 0.84)

gbins = [(3, 10), (1, 3), (0.3, 1), (0.03, 0.3)]
xs, ys, lo_, hi_, yth = [], [], [], [], []
for gl, gh in gbins:
    mZ = fid & (gN_a0 >= gl) & (gN_a0 < gh)
    g, lo, hi = E2_mask(mZ, "emp")
    gt, _, _  = E2_mask(mZ, "thermal")
    xs.append(np.sqrt(gl*gh)); ys.append(g); lo_.append(g-lo); hi_.append(hi-g); yth.append(gt)
fig, ax = plt.subplots(figsize=(6.8, 4.4))
ax.axhspan(1.35, 1.40, color="#c0392b", alpha=0.18, label="MOND-AQUAL + EFE prediction")
ax.axhline(1.0, color="k", ls="--", lw=1.2, label="Newton")
ax.errorbar(xs, ys, yerr=[lo_, hi_], fmt="o", color="#2471a3", capsize=4, ms=7,
            label="E2, empirical per-zone $f(e)$ (68%)")
ax.plot(xs, yth, "s", color="#7d3c98", ms=6, mfc="none",
        label="E2, thermal $f(e)$ (systematic band edge)")
ax.set_xscale("log"); ax.invert_xaxis()
ax.set_xlabel(r"$g_{\rm N}/a_0$  (deep-MOND regime $\rightarrow$)")
ax.set_ylabel(r"recovered $\gamma = G_{\rm eff}/G_{\rm N}$")
ax.set_ylim(0.82, 1.5); ax.legend(fontsize=8, loc="upper left"); ax.grid(alpha=0.25, which="both")
fig.tight_layout(); fig.savefig("fig3_gamma_vs_gN.png", dpi=170); plt.close(fig)
print("\nfigures 2 et 3 sauvees")
for x, y, l, h, t in zip(xs, ys, lo_, hi_, yth):
    print(f"  gN/a0~{x:5.2f} : gamma_emp={y:.3f} +{h:.3f}/-{l:.3f} | thermal={t:.3f}")
