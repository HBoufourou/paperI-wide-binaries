#!/usr/bin/env python3
"""
moteur_population.py — PHASE F : moteur de population pour l'article
« forensique des estimateurs » (binaires pures + triples hiérarchiques
génératives à la Pittordis-Sutherland 2025, sélection et cuts simulés).

VALIDATION (oracles EXTERNES, publiés — pas nos propres chiffres) :
  V-F1  binaires thermiques : médiane vtilde = 0,549 (notre moteur C1, recoupé)
  V-F2  P90(vtilde) binaires ≈ 0,94 ± 0,01, quasi indépendant de f(e)  [PS25 §3.1]
  V-F3  fraction vtilde ≥ 0,8 : flat 23,2 % ; thermique 21,0 % ; γ=1,3 20,6 %  [PS25 §3.1]
  V-F4  triples APRÈS cuts : médiane/P90 par bin r_p ≈ Table 2 de PS25
  V-F5  survie des cuts ≈ 40,5 % (ruwe le plus sévère ~37 % de rejet)  [PS25 §3.3.1]
"""
import numpy as np

G4   = 4*np.pi**2          # au^3 yr^-2 Msun^-1
KMS  = 4.74047             # au/yr -> km/s
T_DR3 = 34/12              # baseline DR3 en années

# ---------------- outils képlériens ----------------
def kepler_E(Ma, e, it=45):
    E = Ma.copy()
    for _ in range(it):
        E -= (E - e*np.sin(E) - Ma)/(1 - e*np.cos(E))
    return E

def orient_matrix(n, rng):
    """angles d'Euler isotropes ; renvoie les 9 coefficients de rotation."""
    Om = rng.uniform(0, 2*np.pi, n)
    w  = rng.uniform(0, 2*np.pi, n)
    ci = rng.uniform(-1, 1, n); si = np.sqrt(1-ci**2)
    cO, sO, cw, sw = np.cos(Om), np.sin(Om), np.cos(w), np.sin(w)
    return cO, sO, ci, si, cw, sw

def to_sky(a_, b_, R):
    """vecteur (a_,b_,0) du plan orbital -> composantes ciel (X,Y) et ligne de visée Z."""
    cO, sO, ci, si, cw, sw = R
    X1 = a_*cw - b_*sw
    Y1 = (a_*sw + b_*cw)*ci
    Z1 = (a_*sw + b_*cw)*si
    return X1*cO - Y1*sO, X1*sO + Y1*cO, Z1

def orbit_plane(a, e, Ma, GM):
    """position (au) et vitesse (au/yr) dans le plan orbital."""
    E  = kepler_E(Ma, e)
    x, y = a*(np.cos(E)-e), a*np.sqrt(1-e**2)*np.sin(E)
    r  = a*(1-e*np.cos(E))
    nmo = np.sqrt(GM/a**3)                       # moyen mouvement
    vx = -a*nmo*np.sin(E)/(1-e*np.cos(E))
    vy =  a*nmo*np.sqrt(1-e**2)*np.cos(E)/(1-e*np.cos(E))
    return x, y, vx, vy, r

# ---------------- distributions ----------------
def sample_e(kind, n, rng, alpha=1.3):
    if kind == "flat":      return rng.uniform(0, 1, n)
    if kind == "thermal":   return np.sqrt(rng.uniform(0, 1, n))
    if kind == "super":     return rng.uniform(0, 1, n)**(1/(1+alpha))
    if kind == "tokovinin": # f(e)=0.4+1.2e -> CDF 0.4e+0.6e^2
        u = rng.uniform(0, 1, n)
        return (-0.4 + np.sqrt(0.16 + 2.4*u))/1.2
    raise ValueError(kind)

def kroupa_pdmf(n, rng, lo=0.18, hi=2.0):
    """IMF Kroupa x (M/0.95)^-2.5 au-dessus de 0.95 (PDMF, PS25 §3.3), continue."""
    def pdf(M):
        p = np.where(M < 0.5, (M/0.5)**-1.3, (M/0.5)**-2.3)
        return p*np.where(M > 0.95, (M/0.95)**-2.5, 1.0)
    out = np.empty(0)
    pmax = pdf(np.array([lo]))[0]
    while len(out) < n:
        M = np.exp(rng.uniform(np.log(lo), np.log(hi), 2*n))
        keep = rng.uniform(0, pmax, 2*n) < pdf(M)/1.0
        out = np.concatenate([out, M[keep]])
    return out[:n]

def mag_G(M, d_pc):
    """G apparent approx : L=M^4 -> M_G ~ 4.83-10log10(M) (suffisant pour les cuts)."""
    return 4.83 - 10*np.log10(M) + 5*np.log10(d_pc/10)

# ---------------- binaires pures (scale-free) ----------------
def binaries_vtilde(n, edist, rng, alpha=1.3, gamma_g=1.0):
    """vtilde pour binaires pures, a=1, G_eff=gamma_g*G, normalisation newtonienne."""
    e  = sample_e(edist, n, rng, alpha)
    Ma = rng.uniform(0, 2*np.pi, n)
    x, y, vx, vy, _ = orbit_plane(np.ones(n), e, Ma, gamma_g*1.0)  # GM=gamma
    R = orient_matrix(n, rng)
    sx, sy, _ = to_sky(x, y, R)
    ux, uy, _ = to_sky(vx, vy, R)
    s2 = np.hypot(sx, sy)
    v2 = np.hypot(ux, uy)
    return v2*np.sqrt(s2)                        # / sqrt(G_N M / s), G_N M=1

# ---------------- triples hiérarchiques (PS25 §3.3) ----------------
def triples(n, rng, cutset="ps25", d_max=300.0, gamma_g=1.0):
    """Génère n triples AVANT cuts ; renvoie dict avec r_p [kau], vtilde,
    masque de survie aux cuts, et diagnostics."""
    # masses + distance + sélection en magnitude
    M1 = kroupa_pdmf(n, rng); M2 = kroupa_pdmf(n, rng)
    q  = rng.uniform(0.02, 1.0, n); M3 = q*M2
    d  = d_max*rng.uniform(0, 1, n)**(1/3)
    ok = (mag_G(M1, d) < 17) & (mag_G(M2, d) < 17)
    idx = np.where(ok)[0]
    M1, M2, M3, q, d = M1[idx], M2[idx], M3[idx], q[idx], d[idx]
    m = len(M1)

    # orbite externe : etoile 1 vs barycentre(2+3)
    a_out = 10**rng.uniform(np.log10(300), np.log10(150000), m)   # au
    e_out = sample_e("super", m, rng, 1.3)
    Ma_o  = rng.uniform(0, 2*np.pi, m)
    xo, yo, vxo, vyo, _ = orbit_plane(a_out, e_out, Ma_o, gamma_g*G4*(M1+M2+M3))
    Ro = orient_matrix(m, rng)
    sxo, syo, _ = to_sky(xo, yo, Ro)
    uxo, uyo, _ = to_sky(vxo, vyo, Ro)
    r_p = np.hypot(sxo, syo)                                       # au

    # orbite interne (2-3) : lognormale Offner + stabilite Tokovinin
    a_in = 10**rng.normal(np.log10(40), 1.5, m)
    a_max = a_out*np.maximum(0.342*(1-e_out)**2, 0.01)
    for _ in range(12):                                            # re-tirage
        bad = a_in > a_max
        if not bad.any(): break
        a_in[bad] = 10**rng.normal(np.log10(40), 1.5, bad.sum())
    a_in = np.minimum(a_in, a_max)
    e_in = sample_e("tokovinin", m, rng)
    P_in = np.sqrt(a_in**3/(M2+M3))                                # yr
    Ma0  = rng.uniform(0, 2*np.pi, m)
    Ri   = orient_matrix(m, rng)
    # positions internes aux deux epoques -> vitesse moyennee (PS25 : 34 mois)
    x0, y0, _, _, _ = orbit_plane(a_in, e_in, Ma0, G4*(M2+M3))
    x1, y1, _, _, _ = orbit_plane(a_in, e_in, Ma0 + 2*np.pi*T_DR3/P_in, G4*(M2+M3))
    sx0, sy0, _ = to_sky(x0, y0, Ri)
    sx1, sy1, _ = to_sky(x1, y1, Ri)
    vinx, viny = (sx1-sx0)/T_DR3, (sy1-sy0)/T_DR3                  # au/yr (ciel)
    theta_in = np.hypot(sx0, sy0)/d                                # arcsec

    # photocentre / centre observable (PS25 eq 17) + masse apparente
    L2, L3 = M2**4, M3**4
    fpb_unres = M3/(M2+M3) - L3/(L2+L3)
    fpb = np.where(theta_in < 1.0, fpb_unres, M3/(M2+M3))
    M2est = np.where(theta_in < 1.0, (L2+L3)**0.25, M2)

    vx_obs = uxo - fpb*vinx
    vy_obs = uyo - fpb*viny
    v2d = np.hypot(vx_obs, vy_obs)*KMS                             # km/s
    vc  = 29.784*np.sqrt((M1+M2est)/r_p)
    vt  = v2d/vc

    # ---- cuts simules ----
    G2app = mag_G(M2, d)
    z = 10**(0.4*(np.maximum(G2app, 14) - 15))
    u = np.sqrt(-1.631 + 680.766*z + 32.732*z**2)
    sigAL = (100 + 7.75*u)/3.0                                     # micro-arcsec
    sig_cen = np.abs(fpb_unres)*a_in/(2*d)*1e6                     # micro-arcsec
    n_orb = np.minimum(1.0, T_DR3/P_in)
    ruwe_sim = np.sqrt(n_orb**4*sig_cen**2 + sigAL**2)/sigAL
    dmag = 2.5*np.log10(L2/L3)

    ruwe_lim = 1.2 if cutset == "ps25" else 1.4
    cut_ruwe = ruwe_sim < ruwe_lim
    cut_ipd  = ~((theta_in > 0.3) & (dmag < 4))
    cut_lob  = q < 0.8
    if cutset == "ps25":
        survive = cut_ruwe & cut_ipd & cut_lob
    else:            # variante "notre echantillon" (El-Badry/Chae) : ruwe seul
        survive = cut_ruwe
    return dict(r_p_kau=r_p/1000, vt=vt, survive=survive, Mapp=M1+M2est,
                cut_ruwe=cut_ruwe, cut_ipd=cut_ipd, cut_lob=cut_lob)

# ================= VALIDATION =================
if __name__ == "__main__":
    rng = np.random.default_rng(11)
    print("="*70)
    print("V-F1/V-F2/V-F3 — BINAIRES PURES (oracles PS25 §3.1) :")
    print(f"  {'f(e)':10s} {'med':>7s} {'P90':>7s} {'frac>=0.8':>10s}   oracle")
    oracles = {"flat": ("0.94", "23.2%"), "thermal": ("0.94", "21.0%"),
               "super": ("0.94", "20.6%")}
    for ed in ["flat", "thermal", "super"]:
        vt = binaries_vtilde(600000, ed, rng)
        print(f"  {ed:10s} {np.median(vt):7.4f} {np.quantile(vt,0.90):7.4f} "
              f"{np.mean(vt>=0.8)*100:9.1f}%   P90~{oracles[ed][0]}, f>=0.8~{oracles[ed][1]}")
    vt_th = binaries_vtilde(600000, "thermal", rng)
    print(f"  V-F1 mediane thermique = {np.median(vt_th):.4f} (attendu 0.549)")

    print("\nV-F4/V-F5 — TRIPLES (oracle Table 2 + §3.3.1 de PS25) :")
    tr = triples(1200000, rng, cutset="ps25")
    sv = tr["survive"]
    print(f"  survie totale : {sv.mean()*100:.1f}%  (oracle ~40.5%)")
    print(f"  rejet ruwe seul : {(~tr['cut_ruwe']).mean()*100:.1f}%  (oracle ~37%)")
    print(f"  rejet lobster : {(tr['cut_ruwe'] & ~tr['cut_lob']).mean()*100:.1f}% suppl.  (oracle ~11%)")
    T2 = {(1.25,1.77): (0.87,1.96), (1.77,2.5): (0.94,2.22), (2.5,3.5): (1.00,2.51),
          (3.5,5.0): (1.10,2.90), (5.0,7.1): (1.20,3.37), (7.1,10.0): (1.32,3.86),
          (10.0,14.1): (1.45,4.43), (14.1,20.0): (1.62,5.16)}
    print(f"  {'bin r_p [kau]':>14s} {'N':>7s} {'med':>6s} {'P90':>6s}   Table2: med  P90")
    for (lo,hi),(mT,pT) in T2.items():
        mbin = sv & (tr["r_p_kau"]>=lo) & (tr["r_p_kau"]<hi)
        v = tr["vt"][mbin]
        print(f"  [{lo:5.2f},{hi:5.2f}) {mbin.sum():7d} {np.median(v):6.2f} "
              f"{np.quantile(v,0.90):6.2f}   Table2: {mT:.2f}  {pT:.2f}")
    print("\nFIN VALIDATION PHASE F")
