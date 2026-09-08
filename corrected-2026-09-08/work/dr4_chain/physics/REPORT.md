# Potentiel conservateur QUMOND à deux masses : résultat de cette étape

Statut des contrôles déclarés : **PASS**. Le modèle sélectionné est `model_refined.npz`. Ce statut porte sur les critères et points prédéfinis, sans borne uniforme de l'erreur physique entre les points.

Deux sources de Plummer prescrites ont été traitées ensemble dans l’action QUMOND non linéaire. Le potentiel d’interaction est calculé par inclusion/exclusion de la configuration complète, de chaque source seule et du champ extérieur. Le modèle newtonien utilise les mêmes deux distributions de masse. Le potentiel relatif est Ψ=U/(m1m2), avec Ψ(∞)=0. Le champ newtonien extérieur uniforme vaut 1, G=Mtotal=a0=1, q∈{0.1,0.3,1}, bi=0.00125√mi et r∈[0.03,30].

Le cosinus μ avec le champ extérieur reste signé. Le potentiel de masses inégales n’est pas artificiellement symétrisé. Cette étape fournit une implémentation pour la chaîne ; elle ne constitue pas une nouvelle théorie physique.

## Précision mesurée

Le maximum de l’écart de force sur les 18 configurations de développement de la grille initiale vaut 0.687324%. Le critère de raffinement fixé était 0.5 %. Après le raffinement prévu, le maximum de développement du modèle choisi vaut 0.0379856%. Le choix de grille a été enregistré avant de calculer les 21 configurations finales.

| q | Maximum final QUMOND | Maximum final Newton Plummer |
|---:|---:|---:|
| 0.1 | 0.0102524% | 0.00021042% |
| 0.3 | 0.00472587% | 0.000214992% |
| 1 | 0.00674426% | 0.000218289% |

Le seuil final est 1 % dans chaque cas. Les références sont les forces intégrées directement avec le solveur antérieur, non modifié. Sur les mêmes points, l’écart maximal entre gradient numérique du potentiel et accélération analytique de l’interpolant vaut 1.01766e-10, et l’asymétrie relative maximale du jacobien vaut 2.38035e-08. Ces deux derniers tests vérifient la conservation numérique du scalaire interpolé, pas l’exactitude indépendante de QUMOND.

| Contrôle d’énergie sur six configurations | Maximum relatif | Repère |
|---|---:|---:|
| Raffinement de la quadrature | 4.02016e-05% | 0.1% |
| Doublement de la frontière spatiale | 5.36958e-05% | 0.1% |

Aux trois anciens cas de référence q∈{0.1,0.3,1}, r=1 et μ=1/√2, l’écart maximal entre l’accélération interpolée et la force de référence plus fine est 0.0305924%.

La correction analytique de queue spatiale est proportionnelle à [νe(1+Ke/3)−1]/R. Elle ne remplace pas l’asymptote orbitale anisotrope dépendant de μ. Doubler la frontière mesure la sensibilité résiduelle ; ce test ne fournit pas à lui seul une borne analytique de tous les termes omis.

## Bornes et usage orbital

| Modèle | q | Borne inférieure de −rΨ | Cmax | Borne extérieure de −rΨ |
|---|---:|---:|---:|---:|
| newton | 0.1 | 0.9991305343 | 0.9999999991 | 0.9999999991 |
| newton | 0.3 | 0.9991279255 | 0.9999999991 | 0.9999999991 |
| newton | 1 | 0.9991260162 | 0.9999999991 | 0.9999999991 |
| qumond | 0.1 | 1.030485179 | 1.63417324 | 1.623536074 |
| qumond | 0.3 | 1.028739983 | 1.6254339 | 1.621267044 |
| qumond | 1 | 1.027849736 | 1.617791899 | 1.617791899 |

Le potentiel et l’accélération dérivent du même spline scalaire bicubique de A=−rΨ en ln r et μ. Les bornes sont obtenues en base Bernstein avec arrondis vers l’extérieur. Elles encadrent le polynôme effectivement enregistré, et peuvent donc entrer dans une sélection fondée sur les invariants E,Lz de ce modèle. Elles ne certifient pas l’erreur physique de QUMOND sur tout le domaine.

L’API refuse les rayons et rapports de masses hors domaine. Son chargement normal exige une validation réussie et les empreintes concordantes du modèle et de son évaluateur. Les q intermédiaires ne sont pas interpolés.

## Revue indépendante et limites restantes

Une autre implémentation, sans import du producteur, calcule l’énergie dans des coordonnées prolates et traite la queue jusqu’à l’infini. Ses gradients aux configurations q=1 et 0.3, r=1, θ=45° concordent avec les forces antérieures. Les résultats détaillés et leur domaine sont enregistrés séparément dans `../review/energy_validation.json`. La reconstruction indépendante du modèle sélectionné a le statut PASS : 49 152 intervalles Bernstein ont été vérifiés en arithmétique rationnelle exacte depuis les nombres enregistrés. Les bornes concernent le polynôme encodé ; elles ne deviennent pas une borne uniforme de l’erreur QUMOND.

Cette étape ne produit aucune population orbitale, aucune distribution astrophysique validée et aucun verdict Newton/MOND sur Gaia. Elle prépare un potentiel conservateur réellement à deux masses dans un champ extérieur prescrit. La construction de populations stationnaires doit employer ce potentiel avec ses invariants, puis déclarer les sélections instantanées et les contaminants. L’ancienne population de trajectoires centrales ne peut pas être réutilisée par simple substitution de cette force anisotrope.

Les sources restent des Plummer finies, les masses sont limitées à trois rapports, le champ extérieur à une valeur et la fonction ν à une forme. Les contrôles portent sur un ensemble fini de points. Une extrapolation à un continuum de masses, au champ galactique variable ou à une validation sur le ciel demanderait des preuves supplémentaires.

## Traçabilité

Modèle SHA256 : `d15c7a58dfda747bc92cf4c67e6a52b4a42bb5cba38438598307f78f5ce2f97f`.
Évaluateur SHA256 : `3551de4294b95d32c1c0b90d2d29b463a6432e821948ad8c07f57cdcbc20cef7`.

Le protocole, les valeurs directes, le choix de grille avant validation finale et les écarts sont conservés. Aucun ancien résultat n’a été modifié.

Sources primaires : [Milgrom2010](https://arxiv.org/html/0911.5464), [Pflamm-Altenburg2025](https://arxiv.org/html/2509.01493v1), [Banik & Zhao](https://arxiv.org/html/1509.08457v3).
