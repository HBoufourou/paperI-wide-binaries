# Portée du contrôle statistique indépendant

Le vérificateur `verify_statistics.py` reconstruit les résultats sans importer `stats/engine.py`, `stats/features.py` ni un module physique ou astrométrique de production. Les observations enregistrées constituent ses entrées. Il ne régénère pas les propositions physiques rejetées ni les observations élémentaires ; le contrôle des banques est une preuve distincte.

Les deux représentations complètes sont recalculées depuis les mesures : positions projetées, masses photométriques, distances, mouvements propres, covariances et diagnostics. Pour Full34, le programme ne reçoit que les quatre champs d’observation de 34 mois ; aucune mesure de 66 mois ni différence entre produits n’est disponible à cette branche. L’inverse des covariances 2×2 est évalué par une expression algébrique indépendante. Les projections utilisent un regroupement explicite des identifiants de cases plutôt que les axes de remodelage de production.

Le programme reconstruit les poids, les probabilités lissées et mélangées des 39 points de nuisance par hypothèse. Le score Conditional66 est calculé directement avec les probabilités après mélange `p(y|z)` et l’entropie au sein des strates. Il ne soustrait ni deux minima ni les surfaces du moteur de production. Les cinq surfaces, minima, différences directionnelles, points optimaux et nombres de solutions proches sont comparés.

Les indices de chacun des 234 cas de calibration et 180 cas de test sont régénérés selon les graines déclarées. Cette partie réutilise volontairement l’algorithme NumPy de tirage annoncé : elle vérifie l’identité du flux, sans prétendre réimplémenter un générateur pseudo-aléatoire. Les histogrammes sont reconstruits indépendamment par regroupement vectorisé et leurs empreintes canoniques sont comparées. Les mêmes indices servent aux trois scénarios, aux cinq variantes et aux deux règles de calibration.

Les rangs supérieurs comptent les égalités. La maximisation sur les 39 nuisances précède la conversion des comptes en probabilités ; les deux statistiques sont combinées au sein de chaque bibliothèque, puis l’enveloppe maximise les trois valeurs combinées. Tous les `p`, rejets, décisions et raisons d’indétermination doivent être exactement égaux aux fichiers produits. Les déviances et probabilités des modèles admettent seulement des tolérances numériques documentées dans le JSON.

Les 5 400 lignes de métriques, les intervalles de Wilson marginaux, les six critères de passage et les 7 020 comparaisons appariées sont recalculés. Les cas purs restent uniques par gravité, réalisation et β. Les empreintes des entrées et sorties lues sont vérifiées avant et après le contrôle. Le nombre de contrôles est un compte de vérifications groupées ; les nombres de valeurs et d’indices comparés sont indiqués séparément.

Le premier essai logiciel complet, à deux catalogues par cas sur des banques artificielles de 32 lignes, a passé 22 896 contrôles ; tous ses `p` et décisions sont identiques. Son pouvoir scientifique est nul : la résolution des rangs impose ici l’indétermination. La première écriture du rapport a rencontré une conversion `int64` vers JSON après les comparaisons ; seul le vérificateur a été corrigé, puis l’essai complet répété. Les deux journaux sont conservés.

La version finale étend ce contrôle aux onze tableaux descriptifs : son essai logiciel complet passe 36 925 contrôles. Les journaux intermédiaires préservent aussi une correction de l’étiquette de portée attendue (`final_methodological_strengthening`), l’intégration des colonnes descriptives ajoutées au writer et une optimisation des produits matriciels sans changement de formule. Ces ajustements concernent uniquement le vérificateur ; aucun paramètre, seuil, source ou résultat scientifique n’a été modifié. Le contrôle scientifique final passe 28 842 contrôles : son nombre de lignes d’échecs descriptifs est plus faible que celui de l’essai artificiel, d’où un total de vérifications groupées plus petit.

Le fichier `statistics_verification.json` est la preuve du contrôle scientifique : seul `complete: true` associé à `passed: true` établit que son parcours exhaustif s’est terminé sans divergence. En l’absence du manifeste final, `--wait` attend explicitement ; il ne produit aucun faux bilan scientifique achevé.

Commande à partir de la racine du dossier reproductible, avec un interpréteur Python disposant de NumPy :

```text
python work/final_paper/review/verify_statistics.py --root work/final_paper --out verification/statistics_verification.json
```

Le choix de `--out` permet une vérification après déplacement sans écrire dans les résultats de production. `--results` et `--manifest` servent notamment à l’essai logiciel isolé. Les chemins historiques absolus inscrits à titre informatif dans le manifeste d’exécution ne sont pas utilisés pour localiser les banques : le contrôleur emploie le manifeste fourni et ses chemins relatifs.

Un PASS numérique ne valide pas les hypothèses astrophysiques, l’échangeabilité entre nouvelles bibliothèques finies, une population continue de nuisances, ni un verdict sur Gaia. La revue est une implémentation séparée avec assistance IA, pas une expertise humaine extérieure. Les résultats restent conditionnels aux bibliothèques de modèles et au processus d’observation déclarés.
