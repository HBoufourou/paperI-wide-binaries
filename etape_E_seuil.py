#!/usr/bin/env python3
"""etape_E_seuil.py — scan de sensibilité du seuil de troncature de E2 (référé, point 1).
Réutilise le pipeline H (données, perspective, gabarits) et fait varier Tcut.
Exécuter depuis la racine du dépôt. Résultat cité dans l'article, Sect. 5."""
import sys, numpy as np; sys.path.insert(0, "code")
src = open("code/phase_H.py").read()
exec(src[:src.index('print("\\nH-2')])          # charge tout jusqu'aux définitions
mV, mT, mD = zone_masks(fid)
rows = []
for em in ["emp", "thermal"]:
    tV = ZoneTemplate(mV, em, seed=3); tZ = ZoneTemplate(mT, em, seed=4)
    vV, vZ = vt_corr[mV], vt_corr[mT]
    for T in [1.2, 1.4, SQ2, 1.6, 1.8]:
        denom = np.median(tV.observed(Tcut=T))
        pred  = np.array([np.median(tZ.observed(g, Tcut=T)) for g in GGRID])/denom
        g_hat = np.interp(np.median(vZ[vZ<T])/np.median(vV[vV<T]), pred, GGRID)
        gb = []
        for _ in range(200):
            bV = rng.choice(vV, len(vV)); bZ = rng.choice(vZ, len(vZ))
            gb.append(np.interp(np.median(bZ[bZ<T])/np.median(bV[bV<T]), pred, GGRID))
        lo, hi = np.quantile(gb, [0.16, 0.84])
        rows.append((em, round(T,3), round(g_hat,3), round(lo,3), round(hi,3)))
        print(f"f(e)={em:8s} Tcut={T:.3f} : gamma={g_hat:.3f} 68%[{lo:.3f},{hi:.3f}]")
import csv
with open("data/etape_E_seuil_resultats.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["fe","Tcut","gamma","lo68","hi68"]); w.writerows(rows)
