# =============================================================
# EXP1 — PHASE C : GRILLE DE ROBUSTESSE (critere C3)
# A executer dans Colab (Drive monte). Colle TOUTE la sortie.
#   C-1  coupures N5 explicites (comptages sequentiels)
#   C-2  propagation des erreurs -> sigma_vtilde par paire
#   C-3  gamma par tranche (s ET g_N/a0) + bootstrap 68%
#   C-4  grille de coupures qualite (stabilite de gamma)
#   C-5  nerf N3 : mini-MC, mediane Newton pour 4 f(e)
#        + excentricites 'e' de Chae par tranche
# =============================================================
import pandas as pd, numpy as np, os

RAW = "/content/drive/MyDrive/EXP1_binaires_gaia/data/raw"
df = pd.read_csv(os.path.join(RAW, "Newton_dr3_MSMS_d200pc_5.csv"))
rng_g = np.random.default_rng(42)
SQ2 = np.sqrt(2)

# ---------- reduction (identique phase B) ----------
d    = 0.5*(df["d1[pc]"] + df["d2[pc]"]).values
dra  = (df["mu1ra[mas/yr]"]  - df["mu2ra[mas/yr]"]).values
ddec = (df["mu1dec[mas/yr]"] - df["mu2dec[mas/yr]"]).values
dmu  = np.hypot(dra, ddec)
v2d  = 4.74047e-3 * dmu * d
Mtot = (df["M1[Msun]"] + df["M2[Msun]"]).values
s_kau = df["s[kau]"].values
s_au  = s_kau * 1000.0
vc   = 29.784*np.sqrt(Mtot/s_au)
vt   = v2d/vc
gN_a0 = 5.930e-3*Mtot/s_au**2 / 1.2e-10        # g_N/a0

# ---------- C-2 : erreurs propagees ----------
s_dra  = np.hypot(df["mu1ra_err[mas/yr]"],  df["mu2ra_err[mas/yr]"]).values
s_ddec = np.hypot(df["mu1dec_err[mas/yr]"], df["mu2dec_err[mas/yr]"]).values
s_dmu  = np.sqrt((dra*s_dra)**2 + (ddec*s_ddec)**2)/np.maximum(dmu, 1e-12)
sig_vt = 4.74047e-3*s_dmu*d / vc               # erreur sur vtilde

# ---------- C-1 : coupures N5 sequentielles ----------
print("="*64)
print("C-1 — COUPURES N5 (comptages sequentiels) :")
m = np.isfinite(vt); print(f"  base finie                : {m.sum()}")
for lab, cut in [("d < 200 pc",            d < 200),
                 ("ruwe1<1.4 & ruwe2<1.4", (df.ruwe1<1.4)&(df.ruwe2<1.4)),
                 ("R_chance < 0.01",        df.R_chance<0.01)]:
    m &= cut.values if hasattr(cut,'values') else cut
    print(f"  + {lab:22s}: {m.sum()}")
base = m.copy()

# ---------- outils ----------
def boot_med(x, nb=500):
    med = np.median(x)
    bs = np.median(rng_g.choice(x, (nb, len(x))), axis=1)
    return med, np.quantile(bs,0.16), np.quantile(bs,0.84)

def table(mask, axis, edges, name):
    mv = mask & (s_kau>=0.2) & (s_kau<2)
    med_v,_,_ = boot_med(vt[mv])
    print(f"  ref validation [0.2-2 kau] : med={med_v:.4f} (N={mv.sum()})")
    print(f"  {name:14s} {'N':>6s} {'med':>7s} {'gamma':>7s} {'68% CI':>15s} {'f>sq2':>6s} {'med(e)':>7s}")
    for lo,hi in edges:
        mb = mask & (axis>=lo) & (axis<hi)
        if mb.sum() < 30: continue
        med,q16,q84 = boot_med(vt[mb])
        g  = (med/med_v)**2
        gl, gh = (q16/med_v)**2, (q84/med_v)**2
        fe = np.mean(vt[mb]>SQ2)*100
        me = np.median(df['e'].values[mb])
        print(f"  [{lo:6.2f},{hi:6.2f}) {mb.sum():6d} {med:7.4f} {g:7.3f} [{gl:6.3f},{gh:6.3f}] {fe:5.1f}% {me:7.3f}")
    return med_v

print("\nC-3 — GAMMA PAR TRANCHE (coupure fiduciaire sig_vt < 0.10) :")
fid = base & (sig_vt < 0.10)
print(f"  N apres coupure qualite : {fid.sum()}")
print("\n  --- axe s [kau] ---")
table(fid, s_kau, [(2,5),(5,10),(10,30)], "s [kau]")
print("\n  --- axe g_N/a0 (le vrai axe de la controverse) ---")
table(fid, gN_a0, [(3,10),(1,3),(0.3,1),(0.03,0.3)], "g_N/a0")

print("\nC-4 — STABILITE DE GAMMA GLOBAL [2-30 kau] SUR LA GRILLE :")
print(f"  {'coupure sig_vt':16s} {'N_test':>7s} {'gamma':>7s} {'68% CI':>15s}")
for sc in [0.05, 0.10, 0.20, np.inf]:
    mq = base & (sig_vt < sc)
    mv = mq & (s_kau>=0.2)&(s_kau<2); mt = mq & (s_kau>=2)&(s_kau<30)
    med_v,_,_ = boot_med(vt[mv]); med_t,q16,q84 = boot_med(vt[mt])
    lab = f"< {sc}" if np.isfinite(sc) else "aucune"
    print(f"  {lab:16s} {mt.sum():7d} {(med_t/med_v)**2:7.3f} "
          f"[{(q16/med_v)**2:6.3f},{(q84/med_v)**2:6.3f}]")

# ---------- C-5 : nerf N3, mediane Newton selon f(e) ----------
print("\nC-5 — MINI-MC : mediane vtilde NEWTON selon la distribution f(e)")
def kepler_E(Ma, e, it=45):
    E = Ma.copy()
    for _ in range(it):
        E -= (E - e*np.sin(E) - Ma)/(1 - e*np.cos(E))
    return E

def med_mc(edist, n=300000, seed=7):
    r = np.random.default_rng(seed)
    e = {"e=0":        np.zeros(n),
         "uniforme":   r.uniform(0,1,n),
         "thermique":  np.sqrt(r.uniform(0,1,n)),          # f(e)=2e
         "super a=1.3":r.uniform(0,1,n)**(1/2.3),          # f(e) prop e^1.3
        }[edist]
    Ma = r.uniform(0,2*np.pi,n); E = kepler_E(Ma, e)
    rr = 1 - e*np.cos(E)                                    # a=1, GM=1
    x, y   = np.cos(E)-e, np.sqrt(1-e**2)*np.sin(E)
    vx, vy = -np.sin(E)/rr, np.sqrt(1-e**2)*np.cos(E)/rr
    inc = np.arccos(r.uniform(-1,1,n)); w = r.uniform(0,2*np.pi,n)
    def proj(px,py):
        X = px*np.cos(w)-py*np.sin(w); Y = (px*np.sin(w)+py*np.cos(w))*np.cos(inc)
        return X, Y
    sx,sy = proj(x,y);  s2 = np.hypot(sx,sy)
    ux,uy = proj(vx,vy); v2 = np.hypot(ux,uy)
    return np.median(v2*np.sqrt(s2))                        # vtilde, GM=1

meds = {}
for ed in ["e=0","uniforme","thermique","super a=1.3"]:
    meds[ed] = med_mc(ed)
    print(f"  f(e) {ed:12s}: med vtilde = {meds[ed]:.4f}")
lo, hi = min(meds.values()), max(meds.values())
print(f"  bande N3 totale : [{lo:.4f}, {hi:.4f}]  soit +/-{(hi/lo-1)*50:.1f}% autour du centre")
print(f"  -> impact maximal sur gamma si f(e) change entre zones : "
      f"x{(hi/lo)**2:.3f} (a comparer au boost MOND attendu 1.35-1.40)")
print("\nFIN PHASE C — copie-colle TOUT.")
