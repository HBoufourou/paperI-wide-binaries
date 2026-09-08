"""Report saved results without altering grids or validation evidence."""
from pathlib import Path
import json,hashlib
HERE=Path(__file__).resolve().parent


def report():
    v=json.loads((HERE/'validation.json').read_text());choice=json.loads((HERE/'results'/'grid_choice_before_final.json').read_text())
    rows=v['rows'];dev=json.loads((HERE/'results'/'development_model_initial.json').read_text())
    selected_dev=json.loads((HERE/'results'/('development_'+Path(choice['selected_model']).stem+'.json')).read_text())
    review_name='interpolant_refined_verification.json' if 'refined' in choice['selected_model'] else 'interpolant_initial_verification.json'
    review_path=HERE.parent/'review'/review_name
    independent_text='La reconstruction indépendante de l’interpolant et des bornes constitue une preuve distincte ; son statut doit être lu dans le dossier de revue.'
    if review_path.exists():
        independent=json.loads(review_path.read_text())
        if independent.get('model_sha256')==v['model_sha256']:
            count=sum(row['exact_bernstein_intervals_checked'] for row in independent['models'])
            independent_text=f"La reconstruction indépendante du modèle sélectionné a le statut {'PASS' if independent['passed'] else 'FAIL'} : {count:,} intervalles Bernstein ont été vérifiés en arithmétique rationnelle exacte depuis les nombres enregistrés. Les bornes concernent le polynôme encodé ; elles ne deviennent pas une borne uniforme de l’erreur QUMOND."
    text=['# Potentiel conservateur QUMOND à deux masses : résultat de cette étape','',
      f"Statut des contrôles déclarés : **{'PASS' if v['passed'] else 'FAIL'}**. Le modèle sélectionné est `{choice['selected_model']}`. Ce statut porte sur les critères et points prédéfinis, sans borne uniforme de l'erreur physique entre les points.",'',
      'Deux sources de Plummer prescrites ont été traitées ensemble dans l’action QUMOND non linéaire. Le potentiel d’interaction est calculé par inclusion/exclusion de la configuration complète, de chaque source seule et du champ extérieur. Le modèle newtonien utilise les mêmes deux distributions de masse. Le potentiel relatif est Ψ=U/(m1m2), avec Ψ(∞)=0. Le champ newtonien extérieur uniforme vaut 1, G=Mtotal=a0=1, q∈{0.1,0.3,1}, bi=0.00125√mi et r∈[0.03,30].','',
      'Le cosinus μ avec le champ extérieur reste signé. Le potentiel de masses inégales n’est pas artificiellement symétrisé. Cette étape fournit une implémentation pour la chaîne ; elle ne constitue pas une nouvelle théorie physique.','',
      '## Précision mesurée','',
      f"Le maximum de l’écart de force sur les 18 configurations de développement de la grille initiale vaut {100*dev['max_force_error']:.6g}%. Le critère de raffinement fixé était 0.5 %. Après le raffinement prévu, le maximum de développement du modèle choisi vaut {100*selected_dev['max_force_error']:.6g}%. Le choix de grille a été enregistré avant de calculer les 21 configurations finales.",'',
      '| q | Maximum final QUMOND | Maximum final Newton Plummer |','|---:|---:|---:|']
    for q in [.1,.3,1.]:
        errors={g:max(r['relative_force_error'] for r in rows if r['q']==q and r['gravity']==g) for g in ['qumond','newton']}
        text.append(f"| {q:g} | {100*errors['qumond']:.6g}% | {100*errors['newton']:.6g}% |")
    text.extend(['',f"Le seuil final est 1 % dans chaque cas. Les références sont les forces intégrées directement avec le solveur antérieur, non modifié. Sur les mêmes points, l’écart maximal entre gradient numérique du potentiel et accélération analytique de l’interpolant vaut {max(r['finite_difference_gradient_error'] for r in rows):.6g}, et l’asymétrie relative maximale du jacobien vaut {max(r['jacobian_asymmetry'] for r in rows):.6g}. Ces deux derniers tests vérifient la conservation numérique du scalaire interpolé, pas l’exactitude indépendante de QUMOND.",'',
       '| Contrôle d’énergie sur six configurations | Maximum relatif | Repère |','|---|---:|---:|',
       f"| Raffinement de la quadrature | {100*max(r['mesh_relative'] for r in v['energy_controls']):.6g}% | 0.1% |",
       f"| Doublement de la frontière spatiale | {100*max(r['outer_relative'] for r in v['energy_controls']):.6g}% | 0.1% |",'',
       f"Aux trois anciens cas de référence q∈{{0.1,0.3,1}}, r=1 et μ=1/√2, l’écart maximal entre l’accélération interpolée et la force de référence plus fine est {100*max(r['relative_error'] for r in v['old_reference_comparisons']):.6g}%.",'',
       'La correction analytique de queue spatiale est proportionnelle à [νe(1+Ke/3)−1]/R. Elle ne remplace pas l’asymptote orbitale anisotrope dépendant de μ. Doubler la frontière mesure la sensibilité résiduelle ; ce test ne fournit pas à lui seul une borne analytique de tous les termes omis.','',
       '## Bornes et usage orbital','',
       '| Modèle | q | Borne inférieure de −rΨ | Cmax | Borne extérieure de −rΨ |','|---|---:|---:|---:|---:|'])
    for g in ['newton','qumond']:
        for i,q in enumerate([.1,.3,1.]):
            text.append(f"| {g} | {q:g} | {v['cmin'][g][i]:.10g} | {v['cmax'][g][i]:.10g} | {v['outer_amax'][g][i]:.10g} |")
    text.extend(['','Le potentiel et l’accélération dérivent du même spline scalaire bicubique de A=−rΨ en ln r et μ. Les bornes sont obtenues en base Bernstein avec arrondis vers l’extérieur. Elles encadrent le polynôme effectivement enregistré, et peuvent donc entrer dans une sélection fondée sur les invariants E,Lz de ce modèle. Elles ne certifient pas l’erreur physique de QUMOND sur tout le domaine.','',
      'L’API refuse les rayons et rapports de masses hors domaine. Son chargement normal exige une validation réussie et les empreintes concordantes du modèle et de son évaluateur. Les q intermédiaires ne sont pas interpolés.','',
      '## Revue indépendante et limites restantes','',
      'Une autre implémentation, sans import du producteur, calcule l’énergie dans des coordonnées prolates et traite la queue jusqu’à l’infini. Ses gradients aux configurations q=1 et 0.3, r=1, θ=45° concordent avec les forces antérieures. Les résultats détaillés et leur domaine sont enregistrés séparément dans `../review/energy_validation.json`. '+independent_text,'',
      'Cette étape ne produit aucune population orbitale, aucune distribution astrophysique validée et aucun verdict Newton/MOND sur Gaia. Elle prépare un potentiel conservateur réellement à deux masses dans un champ extérieur prescrit. La construction de populations stationnaires doit employer ce potentiel avec ses invariants, puis déclarer les sélections instantanées et les contaminants. L’ancienne population de trajectoires centrales ne peut pas être réutilisée par simple substitution de cette force anisotrope.','',
      'Les sources restent des Plummer finies, les masses sont limitées à trois rapports, le champ extérieur à une valeur et la fonction ν à une forme. Les contrôles portent sur un ensemble fini de points. Une extrapolation à un continuum de masses, au champ galactique variable ou à une validation sur le ciel demanderait des preuves supplémentaires.','',
      '## Traçabilité','',
      f"Modèle SHA256 : `{v['model_sha256']}`.",f"Évaluateur SHA256 : `{v['evaluator_sha256']}`.",'',
      'Le protocole, les valeurs directes, le choix de grille avant validation finale et les écarts sont conservés. Aucun ancien résultat n’a été modifié.','',
      'Sources primaires : [Milgrom2010](https://arxiv.org/html/0911.5464), [Pflamm-Altenburg2025](https://arxiv.org/html/2509.01493v1), [Banik & Zhao](https://arxiv.org/html/1509.08457v3).',''])
    (HERE/'REPORT.md').write_text('\n'.join(text).replace('49,152','49 152'),encoding='utf-8')

if __name__=='__main__':report()
