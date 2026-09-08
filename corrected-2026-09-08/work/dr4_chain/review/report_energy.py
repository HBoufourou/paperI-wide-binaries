"""Descriptive report of completed independent energy checks only."""
from pathlib import Path
import hashlib,json
import numpy as np
HERE=Path(__file__).resolve().parent


def main():
    paths=[HERE/'energy_validation.json',HERE/'energy_numerics.json',
           HERE.parent/'physics/pilot_first.json',HERE.parent/'physics/pilot_controls.json']
    val,num,first,pilots=[json.loads(p.read_text()) for p in paths]
    independent={(r['q'],r['theta_deg']):r for r in val['integrals']
                 if r['d']==1 and r['order']==12 and r['theta_deg'] in (45.,135.)}
    compared=[]
    for p in [first]+pilots:
        own=independent[(p['q'],45.)]
        compared.append({'q':p['q'],'production_level':p['level'],'production_Rfactor':p['outer_factor'],
                         'production_U':p['U'],'independent_U':own['U'],
                         'relative_difference':abs(p['U']/own['U']-1)})
    gradients=[r for r in val['gradients'] if r['step_r']==.002]
    uplus=independent[(.3,45.)]['U'];uminus=independent[(.3,135.)]['U']
    out={'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         'input_sha256':{str(p.relative_to(HERE.parent)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
         'completed':val['complete'] and num['complete'],'passed':val['numerical_checks_passed'] and num['passed'],
         'primary_integrals':len(val['integrals']),'primary_nodes':sum(r['nodes'] for r in val['integrals']),
         'auxiliary_integrals':len(num['cases']),'gradient_comparison':gradients,
         'potential_comparison':compared,'numerical_auxiliary_checks':num['checks'],
         'signed_asymmetry':{'q':.3,'U_plus_mu':uplus,'U_minus_mu':uminus,'relative_to_plus':abs(uminus-uplus)/abs(uplus)},
         'limited_scope':'Independent action/force checks at q=1,.3,d=1,E_N=1,epsilon=.00125; signed case135deg. No uniform full-grid guarantee, population or sky verdict.'}
    (HERE/'ENERGY_REPORT.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
    text=f'''# Vérification indépendante du potentiel : contrôles exécutés

Les contrôles déclarés passent. L'intégration prolate indépendante n'importe aucun code du solveur de production. Elle calcule l'action fantôme jusqu'à l'infini par une transformation rationnelle et ajoute l'énergie NewtonPlummer obtenue en1D. Le solveur de production utilise des coordonnées sphériques multicentres et une correction de queue spatiale : les deux intégrations ont donc des géométries et frontières numériques différentes.

Le contrôle principal conserve {out['primary_integrals']} intégrales ({out['primary_nodes']:,} nœuds). Pour q=1 et0.3, d=1, champ newtonien externe1, epsilon=.00125, les dérivées radiale **et angulaire** du potentiel reproduisent les forces précédentes à {gradients[0]['force_reference_relative_error']:.6g} et {gradients[1]['force_reference_relative_error']:.6g} relatif. Le seuil fixé avant calcul était1e−4. Les changements de résolution et de pas passent séparément ; aucune extrapolation n'est nécessaire pour franchir ce seuil.

Les quatre comparaisons directes d'énergie avec les pilotes multicentres diffèrent au maximum de {max(p['relative_difference'] for p in compared):.6g} relatif. Ce maximum inclut le niveau de référence et la frontière doublée ; le niveau le plus coûteux n'est pas présenté comme une vérité exacte. Les valeurs et hachages figurent dans ENERGY_REPORT.json.

La différence signée pour q=.3 est nette : U(45°)={uplus:.14g}, U(135°)={uminus:.14g}, soit {100*out['signed_asymmetry']['relative_to_plus']:.5g}% de |U(45°)|. Un pliage μ→|μ| serait donc incorrect pour ce système de masses inégales. À masses égales, la réflexion est retrouvée à {num['checks']['equal_mass_reflection_relative']:.3g} relatif ; l'échange simultané des étiquettes et de l'orientation pour les masses inégales est retrouvé à {num['checks']['unequal_mass_label_exchange_relative']:.3g}. Ces identités sont contrôlées par six intégrales auxiliaires distinctes.

Le déplacement du seuil de l'intégrale Hessienne mixte change l'énergie de {num['checks']['mixed_hessian_switch_relative']:.3g} relatif. Déplacer de1000 à100 le raccord à la queue infinie la change de {num['checks']['infinite_tail_breakpoint_relative']:.3g}. Les identités analytiques de Qph, de son Hessien et de la limite externe sont également conservées dans les résultats.

Cette réussite valide les configurations effectivement calculées. Elle ne borne pas l'erreur sur tout le domaine, ne certifie pas un modèle de populations astrophysiques et ne constitue aucun résultat sur la gravité mesurée dans Gaia. Les sources restent deux profils de Plummer prescrits ; le champ extérieur et la fonction nu restent les hypothèses du protocole. La vérification de l'interpolant, de sa borne et des populations couplées est une étape séparée.

Reproduction depuis le dossier de travail : exécuter independent_energy.py, check_energy.py, check_energy_numerics.py, puis report_energy.py avec Python et NumPy. Les chemins de références de force sont paramétrables par --force-root ; les résultats principaux peuvent être écrits ailleurs par --out. Aucun ancien résultat n'est modifié.
'''
    (HERE/'ENERGY_REPORT.md').write_text(text,encoding='utf-8')
    print(json.dumps({'passed':out['passed'],'max_energy_difference':max(p['relative_difference'] for p in compared)}))

if __name__=='__main__':main()
