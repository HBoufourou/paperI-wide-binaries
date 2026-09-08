#!/usr/bin/env python3
"""
phase_H.py — MESURE DR3 FINALE pour l'article.
  H-1  correction de perspective (vitesse relative transverse 3D via RV)
  H-2  E2 corrige du biais de Rice : gabarits par zone, BRUITES avec les
       sigma reels de la zone, tronques El-Badry (N2 complet)
  H-3  grille de cuts : fiduciaire / strict-2026 (ruwe<1.2) / +HRD-lobster
  H-4  E3 sur les vraies donnees : vraisemblance jointe (gamma, f_trip),
       gabarit de triples du moteur (cuts echantillon, bruit zone test)
Bandes systematiques : f(e) ∈ {empirique(Chae), thermique, super 1.3}.
"""
import numpy as np, pandas as pd
from moteur_population import (orbit_plane, orient_matrix, to_sky,
                               sample_e, triples)

rng = np.random.default_rng(77)
K = 4.74047e-3            # km/s par (mas/an * pc)
SQ2 = np.sqrt(2)
GGRID = np.linspace(0.7, 1.9, 121)

# ---------------- donnees ----------------
import os
A = pd.read_csv("data/chae_colonnes_utiles.csv")
CPATH = None
for p in ["/mnt/user-data/uploads/chae_complement.csv", "data/chae_complement.csv"]:
    if os.path.exists(p): CPATH = p
if CPATH:
    B = pd.read_csv(CPATH)
    df = A.merge(B, on=["source_id1", "source_id2"], validate="1:1")
    HAS_C = True
else:
    df = A; HAS_C = False
    print("!! complement absent : H-1 (perspective) et HRD sautes, v2d plat utilise")
assert len(df) == 81088, len(df)

d1, d2 = df["d1[pc]"].values, df["d2[pc]"].values
d = 0.5*(d1 + d2)
s_kau = df["s[kau]"].values
Mtot = (df["M1[Msun]"] + df["M2[Msun]"]).values
gN_a0 = 5.930e-3*Mtot/(s_kau*1000)**2/1.2e-10

def sky_basis(ra_deg, dec_deg):
    a, de = np.radians(ra_deg), np.radians(dec_deg)
    rhat = np.stack([np.cos(de)*np.cos(a), np.cos(de)*np.sin(a), np.sin(de)])
    ahat = np.stack([-np.sin(a), np.cos(a), np.zeros_like(a)])
    dhat = np.stack([-np.sin(de)*np.cos(a), -np.sin(de)*np.sin(a), np.cos(de)])
    return rhat, ahat, dhat

if HAS_C:
    # RV : sentinelles -10000/-20000 = manquantes (DR3 : etoiles brillantes seulement)
    RV1, RV2 = df["RV1[km/s]"].values, df["RV2[km/s]"].values
    ok1, ok2 = np.abs(RV1) < 500, np.abs(RV2) < 500
    RVsys = np.where(ok1 & ok2, 0.5*(RV1+RV2), np.where(ok1, RV1, np.where(ok2, RV2, 0.0)))
    has_rv = ok1 | ok2
    # vitesse systemique : mu moyen + RVsys, a la position moyenne
    dRA = (df["RA2[deg]"].values - df["RA1[deg]"].values + 180) % 360 - 180
    RAm = df["RA1[deg]"].values + 0.5*dRA
    DEm = 0.5*(df["DEC1[deg]"].values + df["DEC2[deg]"].values)
    rm, am, dmh = sky_basis(RAm, DEm)
    mua_m = 0.5*(df["mu1ra[mas/yr]"] + df["mu2ra[mas/yr]"]).values
    mud_m = 0.5*(df["mu1dec[mas/yr]"] + df["mu2dec[mas/yr]"]).values
    v_sys = RVsys*rm + K*d*(mua_m*am + mud_m*dmh)
    # mu predit par la vitesse systemique a la position/distance de chaque etoile
    r1h, a1h, d1h = sky_basis(df["RA1[deg]"].values, df["DEC1[deg]"].values)
    r2h, a2h, d2h = sky_basis(df["RA2[deg]"].values, df["DEC2[deg]"].values)
    # distance COMMUNE (moyenne) : ne garder que le terme ANGULAIRE de la
    # perspective ; utiliser d1,d2 individuels injecterait l'erreur de
    # parallaxe amplifiee par v_t systemique (~30 km/s) — piege connu.
    mu1pa = np.sum(v_sys*a1h, 0)/(K*d); mu1pd = np.sum(v_sys*d1h, 0)/(K*d)
    mu2pa = np.sum(v_sys*a2h, 0)/(K*d); mu2pd = np.sum(v_sys*d2h, 0)/(K*d)
    dmu_ca = (df["mu2ra[mas/yr]"] - df["mu1ra[mas/yr]"]).values - (mu2pa - mu1pa)
    dmu_cd = (df["mu2dec[mas/yr]"] - df["mu1dec[mas/yr]"]).values - (mu2pd - mu1pd)
    v2d_corr = K*d*np.hypot(dmu_ca, dmu_cd)                          # H-1
    print(f"RV valides : etoile1 {ok1.sum()}, etoile2 {ok2.sum()}, "
          f"au moins une {has_rv.sum()} ({has_rv.mean()*100:.0f}%)")

# version plate (phase B) pour mesurer l'effet de perspective
dra  = (df["mu1ra[mas/yr]"] - df["mu2ra[mas/yr]"]).values
ddec = (df["mu1dec[mas/yr]"] - df["mu2dec[mas/yr]"]).values
dmu  = np.hypot(dra, ddec)
v2d_flat = K*dmu*d

vc = 0.94179*np.sqrt(Mtot/s_kau)
vt_flat = v2d_flat/vc
vt_corr = (v2d_corr/vc) if HAS_C else vt_flat.copy()
s_dmu = np.sqrt((dra*np.hypot(df["mu1ra_err[mas/yr]"], df["mu2ra_err[mas/yr]"]))**2 +
                (ddec*np.hypot(df["mu1dec_err[mas/yr]"], df["mu2dec_err[mas/yr]"]))**2)/np.maximum(dmu, 1e-12)
sig_vt = K*s_dmu*d/vc

print("="*66)
print("H-1 — CORRECTION DE PERSPECTIVE (RV 3D)%s :" % ("" if HAS_C else " [SAUTEE - complement absent]"))
for nom, m in [("valid [0.2-2]", (s_kau >= 0.2) & (s_kau < 2)),
               ("test  [2-30]", (s_kau >= 2) & (s_kau < 30)),
               ("profond gN<0.3a0", (gN_a0 >= 0.03) & (gN_a0 < 0.3))]:
    dd = vt_corr[m] - vt_flat[m]
    print(f"  {nom:18s} med(dvt)={np.median(dd):+.4f}  P90(|dvt|)={np.quantile(np.abs(dd),0.9):.4f}"
          f"  med(vt): {np.median(vt_flat[m]):.4f} -> {np.median(vt_corr[m]):.4f}")

# ---------------- cuts ----------------
base = np.isfinite(vt_corr) & (d < 200) & (df.R_chance < 0.01).values
fid    = base & (df.ruwe1 < 1.4).values & (df.ruwe2 < 1.4).values & (sig_vt < 0.10)
strict = base & (df.ruwe1 < 1.2).values & (df.ruwe2 < 1.2).values & (sig_vt < 0.10)
if HAS_C:
    MG1 = df.MagG1.values - 5*np.log10(d1/10) - df["A_G1[mag]"].values
    MG2 = df.MagG2.values - 5*np.log10(d2/10) - df["A_G2[mag]"].values
    col = np.concatenate([df.bp_rp1.values, df.bp_rp2.values])
    mg  = np.concatenate([MG1, MG2])
    okc = np.isfinite(col) & np.isfinite(mg) & (col > 0.1) & (col < 3.2)
    cf = np.polyfit(col[okc], mg[okc], 5)
    lob1 = MG1 - np.polyval(cf, df.bp_rp1.values)
    lob2 = MG2 - np.polyval(cf, df.bp_rp2.values)
    hrd = (np.abs(lob1 - lob2) <= 0.40) & (np.minimum(lob1, lob2) >= -0.75)
    print(f"retention HRD-lobster : {hrd.mean()*100:.0f}%% (PS25 : ~58%%)")
else:
    hrd = np.ones(len(df), bool)
strict_hrd = strict & hrd
print(f"\nH-3 — EFFECTIFS : fiduciaire {fid.sum()} | strict ruwe<1.2 {strict.sum()}"
      f" | strict+HRD {strict_hrd.sum()}")

# ---------------- gabarits binaires par zone (bruites, N2) ----------------
def zone_masks(sel):
    return (sel & (s_kau >= 0.2) & (s_kau < 2),
            sel & (s_kau >= 2) & (s_kau < 30),
            sel & (gN_a0 >= 0.03) & (gN_a0 < 0.3))

def binaries_scalefree(n, e_vals, rng):
    e = np.clip(rng.choice(e_vals[np.isfinite(e_vals)], n), 0, 0.995)
    Ma = rng.uniform(0, 2*np.pi, n)
    x, y, vx, vy, _ = orbit_plane(np.ones(n), e, Ma, 1.0)
    R = orient_matrix(n, rng)
    sx, sy, _ = to_sky(x, y, R); ux, uy, _ = to_sky(vx, vy, R)
    s2 = np.hypot(sx, sy)
    return np.hypot(ux, uy)*np.sqrt(s2)

class ZoneTemplate:
    """gabarit binaire : signal scale-free + bruit et masses RESAMPLES de la zone."""
    def __init__(self, zmask, e_model, n=400000, seed=3):
        r = np.random.default_rng(seed)
        if e_model == "emp":  evals = df["e"].values[zmask]
        elif e_model == "thermal": evals = np.sqrt(r.uniform(0, 1, 200000))
        else:                 evals = r.uniform(0, 1, 200000)**(1/2.3)
        self.sig = binaries_scalefree(n, evals, r)
        idx = r.choice(np.where(zmask)[0], n)
        self.svt = sig_vt[idx]; self.M = Mtot[idx]
        self.g1, self.g2 = r.normal(size=n), r.normal(size=n)
    def observed(self, gamma=1.0, Tcut=None):
        v = np.hypot(np.sqrt(gamma)*self.sig + self.g1*self.svt, self.g2*self.svt)
        v = v[v <= 2.23/np.sqrt(self.M)]                      # troncature El-Badry
        return v if Tcut is None else v[v < Tcut]

def E2_final(sel, e_model, Tcut=SQ2, deep=False):
    mV, mT, mD = zone_masks(sel)
    mZ = mD if deep else mT
    tV = ZoneTemplate(mV, e_model, seed=3)
    tZ = ZoneTemplate(mZ, e_model, seed=4)
    denom = np.median(tV.observed(Tcut=Tcut))
    pred = np.array([np.median(tZ.observed(g, Tcut=Tcut)) for g in GGRID])/denom
    vV, vZ = vt_corr[mV], vt_corr[mZ]
    r_obs = np.median(vZ[vZ < Tcut])/np.median(vV[vV < Tcut])
    g_hat = np.interp(r_obs, pred, GGRID)
    gb = []
    for _ in range(300):
        bV = rng.choice(vV, len(vV)); bZ = rng.choice(vZ, len(vZ))
        gb.append(np.interp(np.median(bZ[bZ < Tcut])/np.median(bV[bV < Tcut]), pred, GGRID))
    gb = np.array(gb)
    return g_hat, np.quantile(gb, 0.16), np.quantile(gb, 0.84), mZ.sum()

print("\nH-2 — E2 FINAL (perspective + gabarits bruites par zone, troncature sqrt2) :")
for nom, sel in [("fiduciaire", fid), ("strict 1.2", strict), ("strict+HRD", strict_hrd)]:
    for deep, lab in [(False, "test"), (True, "profond")]:
        g, lo, hi, n = E2_final(sel, "emp", deep=deep)
        print(f"  {nom:11s} {lab:8s} (e empirique) : gamma={g:.3f} 68%[{lo:.3f},{hi:.3f}] N={n}")
print("  bande f(e) (fiduciaire, test) :", end=" ")
for em in ["emp", "thermal", "super"]:
    g, lo, hi, _ = E2_final(fid, em)
    print(f"{em}={g:.3f}", end="  ")
print()

# ---------------- H-4 : vraisemblance jointe (gamma, f_trip) ----------------
print("\nH-4 — FIT DE MELANGE (gamma, f_trip) sur la zone test (fiduciaire) :")
mV, mT, mD = zone_masks(fid)
bins = np.linspace(0, 2.4, 49)                  # la donnee est tronquee El-Badry
hd, _ = np.histogram(vt_corr[mT], bins)
tT = ZoneTemplate(mT, "emp", n=800000, seed=9)
tr = triples(1200000, rng, cutset="sample", d_max=200.0)
ok = tr["survive"] & (tr["r_p_kau"] >= 2) & (tr["r_p_kau"] < 30)
rtr = np.random.default_rng(10)
idx = rtr.choice(np.where(mT)[0], ok.sum())
vtr = tr["vt"][ok] + 0  # signal
vtr = np.hypot(vtr + rtr.normal(size=len(vtr))*sig_vt[idx], rtr.normal(size=len(vtr))*sig_vt[idx])
vtr = vtr[vtr <= 2.23/np.sqrt(tr["Mapp"][ok])]
ht = np.histogram(vtr, bins)[0].astype(float); ht /= ht.sum()
best = (-np.inf, None, None)
FGRID = np.linspace(0, 0.5, 51)
LL = np.zeros((len(GGRID), len(FGRID)))
for i, g in enumerate(GGRID):
    hb = np.histogram(tT.observed(g), bins)[0].astype(float)
    hb /= max(hb.sum(), 1)
    for j, f in enumerate(FGRID):
        model = hd.sum()*((1-f)*hb + f*ht) + 1e-9
        LL[i, j] = np.sum(hd*np.log(model) - model)
imax = np.unravel_index(np.argmax(LL), LL.shape)
g_ml, f_ml = GGRID[imax[0]], FGRID[imax[1]]
prof = LL.max(axis=1)
ci = GGRID[prof >= prof.max() - 0.5]
print(f"  gamma_ML = {g_ml:.3f}  [68% profil : {ci.min():.3f}, {ci.max():.3f}]   f_trip_ML = {f_ml:.3f}")
print(f"  Delta lnL (gamma=1.4 vs ML) = {prof.max() - prof[np.argmin(np.abs(GGRID-1.4))]:.1f}")
print(f"  Delta lnL (gamma=1.0 vs ML) = {prof.max() - prof[np.argmin(np.abs(GGRID-1.0))]:.1f}")
np.save("phase_H_LL.npy", LL)
print("\nFIN PHASE H")
