# Potentiel d'interaction et interpolant conservateur — protocole avant calcul

8 septembre 2026. Nouvelle étape ; aucun résultat antérieur n'est modifié. Un seul processus de calcul lourd est autorisé, avec OMP/BLAS limités à un fil.

## Domaine physique et origine

Deux sources de Plummer prescrites, G=Mtotal=a0=1, q∈{0.1,0.3,1}, mi=(1,q)/(1+q), bi=0.00125√mi. Champ newtonien extérieur (0,0,1), fonction ν simple inchangée. Le vecteur relatif a longueur r et cosinus **signé** μ avec l'axe extérieur. Aucun pliage μ→|μ| n'est appliqué, y compris lors de la construction de la grille. L'interaction relative est Ψ=Uinteraction/(m1m2), avec origine Ψ(∞)=0. Le champ physique extérieur uniforme n'est pas ajouté au mouvement relatif.

L'énergie d'action utilise l'inclusion/exclusion de la source totale, de chaque source seule dans le même champ extérieur et du fond extérieur. La source totale est non linéaire ; on ne superpose pas deux solutions QUMOND. L'énergie est séparée en une interaction newtonienne entre Plummer, intégrée indépendamment en 1D, et une correction Qph=Q−|u|². Les expressions stables de Qph et sa dérivée sont contrôlées. Dans la zone |u1|+|u2|<0.01, une quadrature mixte Gauss2×2 du Hessien de Qph évite la soustraction de quatre valeurs presque égales.

## Quadrature spatiale et queue

La correction est intégrée avec deux quadratures sphériques pondérées par la partition multicentre déjà définie pour les forces. Niveau initial : largeur logarithmique radiale maximale0.6, Gauss radial6,24 nœuds en cosθ,48 azimuts. Niveau de référence :0.45,8,40,80. Frontière intérieure10^-12√mi ; frontière extérieure R=64max(1,r).

La queue spatiale fantôme ajoutée est −m1m2[νe(1+Ke/3)−1]/R. Elle corrige la convergence en1/R de l'énergie. Elle ne remplace pas l'asymptote **orbitale** anisotrope −m1m2νe[1+Ke(1−μ²)/2]/r. La frontière multicentre n'est pas exactement une sphère commune ; doubler R contrôle aussi cette approximation, sans la déclarer rigoureusement bornée.

## Pilotage numérique et grille

Avant la grille complète, mesurer le coût d'un premier cas q1,r1,μ=1/√2 au niveau initial ; calculer aussi le niveau de référence et R doublé. Ce pilote peut révéler une erreur d'implémentation, qui sera corrigée et tracée, mais ne servira pas à sélectionner un résultat gravitationnel.

Grille initiale :17 rayons uniformes en log r entre0.03 et30,9 valeurs uniformes de μ entre−1 et1, pour chacun des trois q, soit459 énergies. L'interpolation porte sur A=−rΨ et sur l'équivalent newtonien. Un spline tensoriel cubique avec conditions «not-a-knot» conserve les polynômes jusqu'au degré3. Ses dérivées définissent le champ ; potentiel et force ne sont pas interpolés indépendamment.

Avec l=ln r, le champ est a_rel=−(A−A_l)rhat/r²+A_μ(ez−μ rhat)/r². Cette expression reste régulière aux pôles. L'API Ψ(positions,q) et acceleration(positions,q) traite des lots, supporte exactement les trois q et refuse les rayons hors[0.03,30]. Les q intermédiaires ne sont pas interpolés sans validation.

## Contrôles de raffinement et validation finale séparée

Points de développement, pour chaque q : (r,μ)=(0.055,−0.83),(0.18,0.28),(0.55,−0.37),(1.4,0.68),(4.7,−0.55),(17,0.1). Les forces directes sont calculées avec le solveur deux-masses antérieur, niveau fin, sans le modifier. Ces18 comparaisons servent uniquement au choix de résolution : si l'écart vectoriel maximal dépasse0.5%, construire la grille raffinée33×17 (valeurs initiales réutilisées), puis réévaluer ces mêmes contrôles. Aucun troisième niveau automatique n'est autorisé : une insuffisance sera rapportée. Le choix de grille sera arrêté avant les validations finales.

Points finaux réservés, pour chaque q : (0.045,0.61),(0.12,−0.22),(0.41,0.87),(0.9,−0.71),(2.8,0.19),(11.5,−0.46),(24,0.79). Le critère final est **moins de1% d'écart vectoriel** aux forces directes dans chacune de ces21 configurations, pour QUMOND et NewtonPlummer. Les critères portent sur des points définis, pas un continuum. Un échec final est conservé et n'entraîne pas de nouveau réglage.

Autres contrôles : comparer les valeurs d'énergie au niveau de référence et avec R doublé aux six configurations q1 et.3, r∈{.05,1,20}, μ=−0.6. Repère0.1% relatif pour chaque différence. Comparer aux trois anciennes forces de référence àr1,μ=+1/√2 et champ extérieur1 (un cas par q ; les références isolées restent hors du présent modèle). Une quadrature d'énergie indépendante est développée séparément et sa portée sera rapportée.

Sur les21 points finaux, vérifier aussi le gradient du potentiel par différences finies (pas10^-5r, repère10^-5 relatif) et l'absence de composante de rotation du champ via symétrie du jacobien (pas10^-4r, repère10^-4). Ces tests sont des contrôles de l'interpolant, sans les assimiler à une validation physique indépendante.

## Borne du modèle interpolé et limites

Chaque patch bicubique de A est converti en base de Bernstein sur[0,1]². Le maximum de ses16 coefficients borne A sur ce patch. Fournir Cmax par q et modèle ainsi que le maximum global, avec arrondi vers le haut. Cette borne concerne exactement le modèle interpolé, pas le champ QUMOND continu inconnu entre les points. Fournir également la borne de A au rayon extérieur pour vérifier les bandes d'énergie admissibles.

Le présent travail fournit un potentiel utilisable par un échantillonneur ; il ne génère aucune population ni aucune décision. La stationnarité, les coupes fondées sur les invariants, les sélections instantanées, la contamination et les observations relèvent du couplage suivant. Une réussite ne démontre ni une population astrophysique correcte ni un test décisif sur Gaia.

Sources primaires : [Milgrom2010, action QUMOND](https://arxiv.org/html/0911.5464), [Pflamm-Altenburg2025, forces à deux masses](https://arxiv.org/html/2509.01493v1), [Banik & Zhao, limite extérieure dominante](https://arxiv.org/html/1509.08457v3).
