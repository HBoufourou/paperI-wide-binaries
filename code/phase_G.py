#!/usr/bin/env python3
"""
phase_G.py — L'EXPÉRIENCE DÉCISIVE : injection-récupération croisée.

On fabrique des catalogues synthétiques à VÉRITÉ CONNUE :
  gamma ∈ {1.0, 1.4} × f(e) ∈ {thermique, super(1.3), dépendant de s}
  × f_triple ∈ {0, 0.1, 0.2, 0.3}
avec sélection El-Badry (vtilde ≤ 2.23/√M), cuts « notre échantillon »
(ruwe<1.4 simulé) et bruit de mouvement propre RÉEL (resamplé des
erreurs du CSV de Chae). Puis on applique LOYALEMENT trois estimateurs :

  E1 « type Chae » : rapport² des médianes pleines test/validation,
     corrigé du rapport newtonien attendu selon les excentricités
     par zone que l'analyste MESURE (il connaît f(e) vrai par zone).
     Sa seule cécité : il ne modélise pas la queue résiduelle de triples.
  E2 « le nôtre »  : médianes tronquées à vtilde<√2 (troncature
     identique appliquée aux gabarits) + inversion de gamma sur grille.
  E3 « type PS »   : fit de mélange (binaires_gamma + gabarit de
     triples) par maximum de vraisemblance Poisson sur l'histogramme
     de la zone test, gamma et f_trip libres.

Chaque estimateur reçoit les gabarits construits avec le BON f(e)
par zone et le MÊME bruit et la MÊME sélection (règle N2 respectée
pour tous) : le seul degré de liberté qui varie est le TRAITEMENT
DE LA QUEUE. La carte de sortie mesure donc exactement le levier N1
(triples) et son couplage avec N3 (excentricités).
"""
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from moteur_population import (orbit_plane, orient_matrix, to_sky,
                               sample_e, triples, KMS)

rng = np.random.default_rng(2026)

# ---------- données réelles : masses et bruit ----------
dat = pd.read_csv("/mnt/user-data/uploads/chae_colonnes_utiles.csv")
Mtot_pool = (dat["M1[Msun]"] + dat["M2[Msun]"]).values
d = 0.5*(dat["d1[pc]"] + dat["d2[pc]"]).values
dra  = (dat["mu1ra[mas/yr]"] - dat["mu2ra[mas/yr]"]).values
ddec = (dat["mu1dec[mas/yr]"] - dat["mu2dec[mas/yr]"]).values
dmu  = np.hypot(dra, ddec)
s_dmu = np.sqrt((dra*np.hypot(dat["mu1ra_err[mas/yr]"], dat["mu2ra_err[mas/yr]"]))**2 +
                (ddec*np.hypot(dat["mu1dec_err[mas/yr]"], dat["mu2dec_err[mas/yr]"]))**2)/np.maximum(dmu,1e-12)
sigv_pool = 4.74047e-3*s_dmu*d          # sigma(v2d) reel [km/s]

def vc_kms(M, r_kau):
    return 0.94179*np.sqrt(M/r_kau)

# ---------- binaires : (s_hat, vtilde) scale-free ----------
def binaries_sv(n, edist, rng, gamma=1.0, e_given=None):
    e = e_given if e_given is not None else sample_e(edist, n, rng)
    Ma = rng.uniform(0, 2*np.pi, n)
    x, y, vx, vy, _ = orbit_plane(np.ones(n), e, Ma, gamma*1.0)
    R = orient_matrix(n, rng)
    sx, sy, _ = to_sky(x, y, R); ux, uy, _ = to_sky(vx, vy, R)
    s_hat = np.hypot(sx, sy)
    return s_hat, np.hypot(ux, uy)*np.sqrt(s_hat)

def add_noise_and_select(r_kau, vt, M, rng, elbadry=True):
    """bruit 2D resample du vrai catalogue + troncature El-Badry."""
    sig_vt = rng.choice(sigv_pool, len(vt)) / vc_kms(M, r_kau)
    vt_obs = np.hypot(vt + sig_vt*rng.normal(size=len(vt)),
                      sig_vt*rng.normal(size=len(vt)))
    keep = vt_obs <= (2.23/np.sqrt(M) if elbadry else np.inf)
    return vt_obs[keep], keep

# ---------- pools de triples (gamma=1 et 1.4), cuts « échantillon » ----------
print("generation des pools de triples (gamma=1 et 1.4)...")
TR = {}
for g in (1.0, 1.4):
    t = triples(1500000, rng, cutset="sample", d_max=200.0, gamma_g=g)
    ok = t["survive"] & (t["r_p_kau"] > 0.05) & (t["r_p_kau"] < 40)
    TR[g] = dict(r=t["r_p_kau"][ok], vt=t["vt"][ok], M=t["Mapp"][ok])
    print(f"  gamma={g}: pool = {ok.sum()} triples post-cuts")

# ---------- construction d'un catalogue synthetique ----------
ZV, ZT = (0.2, 2.0), (2.0, 30.0)
def make_mock(n, f_trip, gamma, fe_truth, rng):
    """MOND-EFE sature : gamma s'applique aux orbites a >= 2 kau
    (basse acceleration), la zone de validation reste newtonienne —
    conforme a la question pre-enregistree."""
    nb = int(n*(1-f_trip)); nt = n - nb
    a = 10**rng.uniform(np.log10(0.1), np.log10(60), nb)      # kau
    if fe_truth == "sdep":
        e = np.where(a < 2.0, sample_e("thermal", nb, rng), sample_e("super", nb, rng))
        s_hat, vt = binaries_sv(nb, None, rng, 1.0, e_given=e)
    else:
        s_hat, vt = binaries_sv(nb, fe_truth, rng, 1.0)
    vt = vt*np.where(a >= 2.0, np.sqrt(gamma), 1.0)           # boost EFE
    r_b = a*s_hat
    M_b = rng.choice(Mtot_pool, nb)
    vtb, kb = add_noise_and_select(r_b, vt, M_b, rng)
    r_b = r_b[kb]
    # triples : pool gamma=1 pour r_p<2 kau, pool gamma pour r_p>=2 kau
    ii1 = rng.choice(len(TR[1.0]["r"]), nt)
    r1, v1, m1 = (TR[1.0][k][ii1] for k in ("r", "vt", "M"))
    if gamma != 1.0:
        wide = r1 >= 2.0
        pool_w = np.where(TR[gamma]["r"] >= 2.0)[0]
        jj = rng.choice(pool_w, wide.sum())
        r1[wide] = TR[gamma]["r"][jj]
        v1[wide] = TR[gamma]["vt"][jj]
        m1[wide] = TR[gamma]["M"][jj]
    vtt, kt = add_noise_and_select(r1, v1, m1, rng)
    r_t = r1[kt]
    return np.concatenate([r_b, r_t]), np.concatenate([vtb, vtt])

# ---------- gabarits binaires par f(e) et zone (bruites, tronques) ----------
def template(fe, rng, n=500000, boost=1.0, Tcut=None, zone=(2.0, 30.0)):
    _, vt = binaries_sv(n, fe, rng)
    vt = boost*vt
    M = rng.choice(Mtot_pool, n)
    r = 10**rng.uniform(np.log10(zone[0]), np.log10(zone[1]), n)
    vto, _ = add_noise_and_select(r, vt, M, rng)
    if Tcut: vto = vto[vto < Tcut]
    return vto

MEDB = {}   # medianes de gabarit, cache
def medB(fe, boost=1.0, Tcut=None, zone=(2.0, 30.0)):
    key = (fe, round(boost, 4), Tcut, zone)
    if key not in MEDB:
        MEDB[key] = np.median(template(fe, np.random.default_rng(5), boost=boost,
                                       Tcut=Tcut, zone=zone))
    return MEDB[key]

# ---------- les trois estimateurs ----------
GGRID = np.linspace(0.7, 1.9, 61)
SQ2 = np.sqrt(2)

def zones(r, vt):
    mV = (r >= ZV[0]) & (r < ZV[1]); mT = (r >= ZT[0]) & (r < ZT[1])
    return vt[mV], vt[mT]

def E1_chae(vV, vT, feV, feT):
    r_obs = np.median(vT)/np.median(vV)
    r_N   = medB(feT, zone=ZT)/medB(feV, zone=ZV)     # gabarits apparies par zone (N2)
    return (r_obs/r_N)**2

def E2_notre(vV, vT, feV, feT):
    r_obs = np.median(vT[vT < SQ2])/np.median(vV[vV < SQ2])
    pred = np.array([medB(feT, boost=np.sqrt(g), Tcut=SQ2, zone=ZT)
                     for g in GGRID])/medB(feV, Tcut=SQ2, zone=ZV)
    return np.interp(r_obs, pred, GGRID)

def E3_ps(vT, feT, rng):
    bins = np.linspace(0, 3.2, 41)
    hd, _ = np.histogram(vT, bins)
    # gabarit binaire par gamma (boost exact) + gabarit triple (zone test)
    base = template(feT, np.random.default_rng(7), n=800000, zone=ZT)
    mtr = (TR[1.0]["r"] >= ZT[0]) & (TR[1.0]["r"] < ZT[1])
    trt, _ = add_noise_and_select(TR[1.0]["r"][mtr], TR[1.0]["vt"][mtr],
                                  TR[1.0]["M"][mtr], np.random.default_rng(8))
    ht = np.histogram(trt, bins)[0].astype(float); ht /= ht.sum()
    best, gbest = -np.inf, 1.0
    for g in GGRID:
        hb = np.histogram(np.sqrt(g)*base, bins)[0].astype(float)
        hb /= max(hb.sum(), 1)
        for f in np.linspace(0, 0.5, 26):
            model = hd.sum()*((1-f)*hb + f*ht) + 1e-9
            ll = np.sum(hd*np.log(model) - model)
            if ll > best: best, gbest = ll, g
    return gbest

# ---------- la grille ----------
print("\ngrille d'injection-recuperation :")
rows = []
FE = {"thermal": ("thermal", "thermal"), "super": ("super", "super"),
      "sdep":   ("thermal", "super")}     # f(e) (zoneV, zoneT) que mesure l'analyste
for fe_truth, (feV, feT) in FE.items():
    for gamma in (1.0, 1.4):
        for ft in (0.0, 0.1, 0.2, 0.3):
            r, vt = make_mock(300000, ft, gamma, fe_truth, rng)
            vV, vT = zones(r, vt)
            g1 = E1_chae(vV, vT, feV, feT)
            g2 = E2_notre(vV, vT, feV, feT)
            g3 = E3_ps(vT, feT, rng)
            rows.append(dict(fe=fe_truth, gamma=gamma, f_trip=ft,
                             E1=g1, E2=g2, E3=g3, NV=len(vV), NT=len(vT)))
            print(f"  fe={fe_truth:8s} verite gamma={gamma:.1f} f_trip={ft:.1f} -> "
                  f"E1(Chae)={g1:5.3f}  E2(notre)={g2:5.3f}  E3(PS)={g3:5.3f}")
res = pd.DataFrame(rows)
res.to_csv("phase_G_resultats.csv", index=False)

# ---------- figure maitresse ----------
fig, axes = plt.subplots(2, 3, figsize=(13, 7), sharex=True, sharey="row")
titles = {"thermal": "f(e) thermique (les 2 zones)",
          "super": "f(e) superthermique (les 2 zones)",
          "sdep": "f(e) dépendant de s (therm.→super)"}
for j, fe in enumerate(["thermal", "super", "sdep"]):
    for i, gamma in enumerate([1.0, 1.4]):
        ax = axes[i, j]
        sub = res[(res.fe == fe) & (res.gamma == gamma)]
        for est, lab, c, mk in [("E1", "E1 médiane pleine (type Chae)", "#c0392b", "o"),
                                ("E2", "E2 tronqué √2 (le nôtre)", "#2471a3", "s"),
                                ("E3", "E3 fit de mélange (type PS)", "#1e8449", "^")]:
            ax.plot(sub.f_trip, sub[est], mk+"-", color=c, label=lab)
        ax.axhline(gamma, color="k", ls="--", lw=1, label=f"vérité γ={gamma}")
        ax.axhline(1.35, color="gray", ls=":", lw=1)
        if i == 0: ax.set_title(titles[fe], fontsize=10)
        if j == 0: ax.set_ylabel(f"γ récupéré (vérité {gamma})")
        if i == 1: ax.set_xlabel("fraction de triples résiduelle")
        ax.grid(alpha=0.3)
axes[0, 0].legend(fontsize=7, loc="upper left")
fig.suptitle("Carte du pseudo-signal : trois estimateurs face au même univers synthétique "
             "(sélection El-Badry + bruit Gaia réel)", fontsize=11)
fig.tight_layout()
fig.savefig("figure_maitresse_G.png", dpi=160)
fig.savefig("figure_maitresse_G.pdf")
print("\nfigure sauvee : figure_maitresse_G.png/pdf ; table : phase_G_resultats.csv")
