# =============================================================
# EXP1 — PHASE D : VERDICT (criteres C2 complet + C4)
# A executer dans Colab (Drive monte). Colle TOUTE la sortie.
#   D-1  med(e) par zone (validation vs test vs bin profond)
#   D-2  MC Newton a excentricites EMPIRIQUES (colonne e de Chae,
#        tiree par zone) -> rapport newtonien attendu r_N
#   D-3  gamma corrige N3 : (r_obs/r_N)^2, bootstrap 68%
#   D-4  estimateur TRONQUE anti-triples (vt<1.2 et vt<sqrt(2)),
#        troncature identique donnees/MC (N2), inversion de gamma
#        sur grille -> gamma_hat par zone et par troncature
# =============================================================
import pandas as pd, numpy as np, os

RAW = "/content/drive/MyDrive/EXP1_binaires_gaia/data/raw"
df = pd.read_csv(os.path.join(RAW, "Newton_dr3_MSMS_d200pc_5.csv"))
rng_g = np.random.default_rng(42); SQ2 = np.sqrt(2)

# ---------- reduction + coupures fiduciaires (phases B/C) ----------
d    = 0.5*(df["d1[pc]"] + df["d2[pc]"]).values
dra  = (df["mu1ra[mas/yr]"]  - df["mu2ra[mas/yr]"]).values
ddec = (df["mu1dec[mas/yr]"] - df["mu2dec[mas/yr]"]).values
dmu  = np.hypot(dra, ddec)
v2d  = 4.74047e-3*dmu*d
Mtot = (df["M1[Msun]"] + df["M2[Msun]"]).values
s_kau = df["s[kau]"].values
vc   = 29.784*np.sqrt(Mtot/(s_kau*1000.0))
vt   = v2d/vc
gN_a0 = 5.930e-3*Mtot/(s_kau*1000.0)**2 / 1.2e-10
s_dra  = np.hypot(df["mu1ra_err[mas/yr]"],  df["mu2ra_err[mas/yr]"]).values
s_ddec = np.hypot(df["mu1dec_err[mas/yr]"], df["mu2dec_err[mas/yr]"]).values
sig_vt = 4.74047e-3*np.sqrt((dra*s_dra)**2+(ddec*s_ddec)**2)/np.maximum(dmu,1e-12)*d/vc
e_ch = df["e"].values

fid = (np.isfinite(vt) & (d<200) & (df.ruwe1<1.4) & (df.ruwe2<1.4)
       & (df.R_chance<0.01) & (sig_vt<0.10)).values

ZV = fid & (s_kau>=0.2) & (s_kau<2)                # validation
ZT = fid & (s_kau>=2)   & (s_kau<30)               # test global
ZD = fid & (gN_a0>=0.03) & (gN_a0<0.3)             # bin profond

print("="*64)
print("D-1 — EXCENTRICITES DE CHAE PAR ZONE :")
for nom, z in [("validation [0.2-2 kau]",ZV),("test [2-30 kau]",ZT),
               ("profond g_N<0.3 a0",ZD)]:
    ee = e_ch[z]
    print(f"  {nom:24s} N={z.sum():6d}  med(e)={np.median(ee):.3f} "
          f" q16={np.quantile(ee,0.16):.3f} q84={np.quantile(ee,0.84):.3f}")

# ---------- moteur MC (identique phase C, e en argument) ----------
def kepler_E(Ma, e, it=45):
    E = Ma.copy()
    for _ in range(it):
        E -= (E - e*np.sin(E) - Ma)/(1 - e*np.cos(E))
    return E

def mc_vt(e_emp, n=400000, seed=7):
    r = np.random.default_rng(seed)
    e = r.choice(e_emp[np.isfinite(e_emp)], n)
    e = np.clip(e, 0, 0.995)
    Ma = r.uniform(0,2*np.pi,n); E = kepler_E(Ma, e)
    rr = 1 - e*np.cos(E)
    x, y   = np.cos(E)-e, np.sqrt(1-e**2)*np.sin(E)
    vx, vy = -np.sin(E)/rr, np.sqrt(1-e**2)*np.cos(E)/rr
    inc = np.arccos(r.uniform(-1,1,n)); w = r.uniform(0,2*np.pi,n)
    def proj(px,py):
        return px*np.cos(w)-py*np.sin(w), (px*np.sin(w)+py*np.cos(w))*np.cos(inc)
    sx,sy = proj(x,y); ux,uy = proj(vx,vy)
    return np.hypot(ux,uy)*np.sqrt(np.hypot(sx,sy))         # vtilde Newton

vtN_V = mc_vt(e_ch[ZV], seed=7)      # Newton, e empirique zone validation
vtN_T = mc_vt(e_ch[ZT], seed=8)      # Newton, e empirique zone test
vtN_D = mc_vt(e_ch[ZD], seed=9)      # Newton, e empirique bin profond

print("\nD-2 — RAPPORT NEWTONIEN ATTENDU (MC a e empiriques) :")
mV, mT, mD = np.median(vtN_V), np.median(vtN_T), np.median(vtN_D)
print(f"  med MC Newton : valid={mV:.4f}  test={mT:.4f}  profond={mD:.4f}")
rN_T, rN_D = mT/mV, mD/mV
print(f"  r_N attendu (Newton pur) : test/valid={rN_T:.4f}  profond/valid={rN_D:.4f}")

def boot_med(x, nb=600):
    bs = np.median(rng_g.choice(x, (nb, len(x))), axis=1)
    return np.median(x), bs

print("\nD-3 — GAMMA CORRIGE N3 (mediane pleine) :")
medV, bsV = boot_med(vt[ZV])
for nom, z, rN in [("test [2-30 kau]", ZT, rN_T), ("profond g_N<0.3a0", ZD, rN_D)]:
    medZ, bsZ = boot_med(vt[z])
    g  = (medZ/medV/rN)**2
    gb = (bsZ/np.median(bsV)/rN)**2
    print(f"  {nom:20s} r_obs={medZ/medV:.4f}  gamma_corr={g:.3f} "
          f"68%[{np.quantile(gb,0.16):.3f},{np.quantile(gb,0.84):.3f}]")

# ---------- D-4 : estimateur tronque + inversion ----------
print("\nD-4 — ESTIMATEUR TRONQUE (anti-triples, troncature identique MC) :")
gammas = np.linspace(0.7, 1.8, 111)
for Tcut in [1.2, SQ2]:
    # courbe de prediction : ratio des medianes tronquees en fonction de gamma
    medV_N = np.median(vtN_V[vtN_V < Tcut])
    for nom, z, vtN_Z in [("test", ZT, vtN_T), ("profond", ZD, vtN_D)]:
        pred = np.array([np.median((np.sqrt(g)*vtN_Z)[np.sqrt(g)*vtN_Z < Tcut])
                         for g in gammas]) / medV_N
        # observe (troncature identique)
        oV = vt[ZV]; oZ = vt[z]
        r_obs = np.median(oZ[oZ<Tcut]) / np.median(oV[oV<Tcut])
        g_hat = np.interp(r_obs, pred, gammas)
        # bootstrap
        gb = []
        for _ in range(400):
            bV = rng_g.choice(oV, len(oV)); bZ = rng_g.choice(oZ, len(oZ))
            rb = np.median(bZ[bZ<Tcut])/np.median(bV[bV<Tcut])
            gb.append(np.interp(rb, pred, gammas))
        gb = np.array(gb)
        print(f"  Tcut={Tcut:.3f}  zone {nom:8s}: r_obs={r_obs:.4f}  "
              f"gamma_hat={g_hat:.3f}  68%[{np.quantile(gb,0.16):.3f},"
              f"{np.quantile(gb,0.84):.3f}]")

print("\nRAPPEL DES HYPOTHESES EN CONCURRENCE :")
print("  Newton            : gamma = 1.00")
print("  MOND-AQUAL + EFE  : gamma ~ 1.35-1.40")
print("\nFIN PHASE D — copie-colle TOUT ; ensuite on redige le verdict C4.")
