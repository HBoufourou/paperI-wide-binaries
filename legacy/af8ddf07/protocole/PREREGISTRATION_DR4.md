# PROTOCOLE PRÉ-ENREGISTRÉ — TEST DE GRAVITÉ DES BINAIRES LARGES AVEC GAIA DR4
Version 1.0 — rédigé le 15/07/2026, AVANT la publication de Gaia DR4.
À geler : dépôt horodaté (Zenodo) du présent document + code (hash SHA-256) avant DR4.

## 1. Engagement

Le pipeline, les estimateurs, les coupures, les zones et les critères de décision
ci-dessous sont définis avant tout accès aux données DR4 et ne seront pas modifiés.
Toute déviation sera signalée comme telle dans la publication. Le résultat sera
publié quel qu'il soit (critère C4).

## 2. Question et hypothèses

γ = G_eff/G_N dans les binaires à basse accélération interne (g_N < a₀ = 1,2×10⁻¹⁰ m/s²).
  H0 (Newton/RG) : γ = 1,00.
  H1 (MOND-AQUAL + effet de champ externe, g_ext ≈ 1,5-1,9 a₀) : γ ≈ 1,35-1,40.
Observable : ṽ = Δv_⊥ / √(G_N·M_tot/s_2D), Δv_⊥ = composante transverse de la
vitesse relative 3D (perspective corrigée via les vitesses radiales).

## 3. Échantillon (DR4)

- Catalogue de paires type El-Badry+ mis à jour DR4 (ou reconstruction équivalente
  documentée) ; d < 300 pc ; G < 18 ; |b| > 15°.
- Coupures : RUWE < 1,2 (les deux composantes) ; ipd_frac_multi_peak ≤ 2 ;
  R_chance < 0,01 ; cut HRD/lobster (corps corrélé : |Δlobster| ≤ 0,25 et
  sur-luminosité < 0,4 mag) ; σ_ṽ < 0,10 (propagée des erreurs de μ).
- Toute coupure de sélection du catalogue amont (ex. troncature Δv) est APPLIQUÉE
  À L'IDENTIQUE aux gabarits simulés (règle N2), y compris le bruit de mesure,
  resamplé PAR ZONE (leçon du biais de Rice, phase G 2026).

## 4. Zones (pré-enregistrées)

  Validation : s_2D ∈ [0,2 ; 2] kau (régime newtonien — calibration obligatoire).
  Test global : s_2D ∈ [2 ; 30] kau.
  Bin profond : g_N/a₀ ∈ [0,03 ; 0,3] (régime MOND saturé — le juge de paix).

## 5. Estimateurs gelés (validés par injection-récupération, phase G 2026)

E2 (primaire) — rapport des médianes TRONQUÉES à ṽ < √2 (troncature identique
   données/gabarits), inversé en γ par gabarits par zone : signal scale-free
   képlérien × f(e), masses/σ/troncature resamplées de la zone. Biais mesuré sur
   vérité connue : ≤ +0,04 (Newton + 30 % triples) ; récupère 1,31-1,36 pour 1,4 injecté.
E3 (secondaire) — vraisemblance Poisson jointe (γ, f_trip) sur l'histogramme de ṽ
   de la zone test [0 ; 2,4], modèle = (1−f)·binaires_γ + f·gabarit de triples
   génératives (Offner/Tokovinin/photocentre/moyennage ΔT_DR4/biais de masse,
   coupures simulées). Biais mesuré : 0,00 (gabarit exact) ; dégénérescence γ↔f
   notée à f→0.
Bande f(e) systématique OBLIGATOIRE : {empirique par zone, thermique, superthermique
   α=1,3} — publiée comme fourchette, pas moyennée.

## 6. Critères de décision (gelés)

  D1. Contrôle : médiane ṽ de la zone de validation reproduite par le gabarit
      newtonien à < 2 % — sinon STOP et publication de l'échec du pipeline.
  D2. Verdict Newton : γ(E2, bin profond) et γ(E3) tous deux < 1,15 avec
      1,35 exclu à > 5σ (stat) et hors de la bande systématique complète.
  D3. Verdict anomalie : γ(E2, bin profond) et γ(E3) tous deux > 1,20 avec
      1,00 exclu à > 5σ et hors bande systématique.
  D4. Zone grise : tout autre résultat — publié comme non concluant avec la
      carte d'injection-récupération refaite sur les conditions DR4.

## 7. Prévision DR4 (pour dimensionner D2/D3)

Précision μ ×2-4 vs DR3 (loi t^{-3/2}) ; échantillon étendu G<18, d<300-400 pc :
N_test attendu ≳ 3-5× DR3 (≳ 25-40 k paires en zone test après coupures).
σ_stat(γ, E3) attendue ≲ 0,005 ; les systématiques (f(e), modèle de triples,
perspective) domineront — d'où la bande obligatoire du §5 et la carte du §6-D4.
Étalon DR3 (2026, ce travail) : γ_test = 1,04-1,05 ± 0,01 (stat), bande f(e)
[0,955 ; 1,042] ; bin profond 0,97 ± 0,13 ; f_trip ajusté 0,14 ;
Δln L(γ=1,4) = 138 (≈16σ).

## 8. Livraison

Publication des résultats bruts (histogrammes par zone), du code, des gabarits,
de la grille d'injection-récupération DR4, et du présent document inchangé.
