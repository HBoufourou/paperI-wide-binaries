# =============================================================
# EXP1 — PHASE A+B : INSPECTION + REDUCTION vtilde
# A executer dans Colab (Drive monte). Colle TOUTE la sortie.
# Le script s'auto-adapte aux noms de colonnes de Chae ; s'il ne
# reconnait pas tout, il liste ce qui manque (et on ajuste).
# Unites verifiees :
#   v_2D [km/s] = 4.74047e-3 * dmu[mas/an] * d[pc]
#   v_c  [km/s] = 29.784 * sqrt(Mtot[Msun] / s[au])
# =============================================================
import pandas as pd, numpy as np, os

RAW = "/content/drive/MyDrive/EXP1_binaires_gaia/data/raw"
CSV = os.path.join(RAW, "Newton_dr3_MSMS_d200pc_5.csv")
OUT = os.path.join(RAW, "..", "vtilde_reduit_v1.csv")

df = pd.read_csv(CSV)
print("=" * 64)
print(f"A1 — {os.path.basename(CSV)} : {len(df)} lignes, {len(df.columns)} colonnes")
print("COLONNES :", list(df.columns))

# ---------- detection auto des colonnes ----------
def find(cands):
    lc = {c.lower(): c for c in df.columns}
    for cand in cands:
        if cand.lower() in lc:
            return lc[cand.lower()]
    # match partiel en dernier recours
    for cand in cands:
        hits = [c for c in df.columns if cand.lower() in c.lower()]
        if len(hits) == 1:
            return hits[0]
    return None

need = {
  "pmra1":  ["pmra_A","pmra1","pmraA","pmra_a","pm_ra1"],
  "pmdec1": ["pmdec_A","pmdec1","pmdecA","pm_dec1"],
  "pmra2":  ["pmra_B","pmra2","pmraB","pm_ra2"],
  "pmdec2": ["pmdec_B","pmdec2","pmdecB","pm_dec2"],
  "M1":     ["M_A","mass1","MassA","M1","mag_mass_A","MA"],
  "M2":     ["M_B","mass2","MassB","M2","MB"],
}
# distance : parallaxe OU colonne distance
opt = {
  "plx1":  ["parallax_A","parallax1","plx_A","plx1"],
  "plx2":  ["parallax_B","parallax2","plx_B","plx2"],
  "d":     ["d_A","dist","distance","d_pc","dA","r_med_A"],
  "s_au":  ["sep_AU","s_AU","sep_au","proj_sep","sep"],
  "theta": ["theta_arcsec","theta","sep_arcsec","ang_sep"],
  # au cas ou Chae fournit deja l'observable :
  "vt_chae": ["vtilde","v_tilde","vp_vc","x_obs"],
  "v2d_chae":["v_p","vp","dV","delta_v","v2d"],
}
col = {k: find(v) for k, v in {**need, **opt}.items()}
print("\nA2 — MAPPING DETECTE :")
for k, v in col.items():
    print(f"  {k:9s} -> {v}")

missing = [k for k in need if col[k] is None]
if missing:
    print("\n!!! MAPPING MANQUANT :", missing)
    print("Colle la sortie complete dans le chat, j'ajuste le script.")
else:
    # ---------- distance ----------
    if col["d"]:
        d = df[col["d"]].astype(float).values
        if np.nanmedian(d) < 50:  # probablement en parallaxe ou kpc ? garde-fou
            print("\n[garde-fou] mediane 'distance' < 50 : verifier l'unite !")
    elif col["plx1"]:
        p1 = df[col["plx1"]].astype(float).values
        p2 = df[col["plx2"]].astype(float).values if col["plx2"] else p1
        d = 1000.0 / (0.5 * (p1 + p2))      # mas -> pc
    else:
        d = None
        print("\n!!! ni distance ni parallaxe trouvees")

    # ---------- separation projetee en au ----------
    if col["s_au"]:
        s_au = df[col["s_au"]].astype(float).values
        if np.nanmedian(s_au) < 100:   # sans doute en arcsec, pas en au
            print("[garde-fou] 'sep' mediane < 100 : interpretee en ARCSEC")
            s_au = s_au * d
    elif col["theta"] and d is not None:
        s_au = df[col["theta"]].astype(float).values * d   # arcsec * pc = au
    else:
        s_au = None
        print("!!! separation introuvable")

    if d is not None and s_au is not None:
        dmu = np.sqrt((df[col["pmra1"]] - df[col["pmra2"]])**2 +
                      (df[col["pmdec1"]] - df[col["pmdec2"]])**2).values  # mas/an
        v2d = 4.74047e-3 * dmu * d                      # km/s
        Mtot = (df[col["M1"]] + df[col["M2"]]).astype(float).values
        vc = 29.784 * np.sqrt(Mtot / s_au)              # km/s
        vt = v2d / vc
        s_kau = s_au / 1000.0

        print("\n" + "=" * 64)
        print("B1 — SANITE PHYSIQUE (ordres de grandeur attendus) :")
        print(f"  d [pc]     : med={np.nanmedian(d):.1f}   (attendu < 200)")
        print(f"  s [kau]    : med={np.nanmedian(s_kau):.3f} min={np.nanmin(s_kau):.3f} max={np.nanmax(s_kau):.2f}")
        print(f"  Mtot [Msun]: med={np.nanmedian(Mtot):.2f}")
        print(f"  v2d [km/s] : med={np.nanmedian(v2d):.3f}  (attendu ~0.2-2)")
        print(f"  vtilde     : med={np.nanmedian(vt):.4f}  frac>sqrt(2)={np.nanmean(vt>np.sqrt(2))*100:.1f}%")

        print("\nB2 — vtilde PAR TRANCHE (ref MC Newton thermique : med=0.549, P95=1.038) :")
        print(f"  {'s [kau]':12s} {'N':>6s} {'med(vt)':>9s} {'P95':>7s} {'f>sqrt2':>8s}")
        for lo, hi in [(0.2,0.5),(0.5,1),(1,2),(2,5),(5,10),(10,30)]:
            m = (s_kau >= lo) & (s_kau < hi) & np.isfinite(vt)
            if m.sum() > 5:
                print(f"  [{lo:4.1f},{hi:4.1f}) {m.sum():6d} {np.median(vt[m]):9.4f} "
                      f"{np.quantile(vt[m],0.95):7.4f} {np.mean(vt[m]>np.sqrt(2))*100:7.1f}%")

        mv = (s_kau>=0.2)&(s_kau<2)&np.isfinite(vt)
        mt = (s_kau>=2)&(s_kau<30)&np.isfinite(vt)
        med_v, med_t = np.median(vt[mv]), np.median(vt[mt])
        print("\nB3 — ESTIMATEUR BRUT (auto-calibre, insensible au 1er ordre a f(e)) :")
        print(f"  med validation [0.2-2 kau] = {med_v:.4f}  (N={mv.sum()})")
        print(f"  med test      [2-30 kau]  = {med_t:.4f}  (N={mt.sum()})")
        print(f"  rapport = {med_t/med_v:.4f}  ->  gamma_brut = {(med_t/med_v)**2:.4f}")
        print(f"  (Newton attend 1.00 ; MOND-EFE attend ~1.35-1.4 ; SANS modele de triples = borne haute)")

        out = pd.DataFrame({"s_kau": s_kau, "vtilde": vt, "Mtot": Mtot,
                            "v2d_kms": v2d, "d_pc": d})
        out.to_csv(OUT, index=False)
        print(f"\nB4 — sauvegarde : {os.path.abspath(OUT)} ({len(out)} lignes)")

print("\nFIN — copie-colle TOUTE la sortie dans le chat.")
