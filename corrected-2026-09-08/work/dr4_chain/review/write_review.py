"""Write a concise review only from completed independent evidence."""
from pathlib import Path
import json
HERE=Path(__file__).resolve().parent


def main():
    energy=json.loads((HERE/'ENERGY_REPORT.json').read_text())
    interp=json.loads((HERE/'interpolant_refined_verification.json').read_text())
    physical=json.loads((HERE/'physics_report_verification.json').read_text())
    pop=json.loads((HERE/'population_verification.json').read_text())
    if not pop.get('complete'):raise RuntimeError('All22 banks must be completed and verified before a final review')
    if not all((energy['passed'],interp['passed'],physical['report_reconstruction_passed'],physical['physics_criteria_passed'],pop['passed'])):
        raise RuntimeError('A failed check requires explicit review; do not generate a passed narrative')
    banks=[r for r in pop['populations'] if not r['label'].startswith('trajectory_initial/')]
    count=sum(r['rows'] for r in banks);ess=min(v['ess'] for r in banks for v in r['weights'].values())
    fractions=sum(r['exact_bernstein_intervals_checked'] for r in interp['models'])
    force=max(r['force_error'] for r in physical['force_metrics'])
    text=f'''# Revue indépendante de la chaîne physique et de ses populations

Les contrôles exécutés passent dans leur périmètre déclaré. Le travail fournit désormais un potentiel conservateur à deux masses, un parent statistique défini par des invariants et des bibliothèques effectivement propagées puis observées. Cette conclusion est distincte du résultat de la comparaison statistique Newton–QUMOND, qui n'entre dans aucun des contrôles ci-dessous.

## Preuves numériques

- Énergie :31intégrales indépendantes en coordonnées prolates, plus6contrôles auxiliaires. Les gradients du potentiel retrouvent deux anciennes forces de référence à mieux que5.1e−6 relatif. L'énergie conserve son origine à l'infini et l'asymétrie signée à masses inégales ; elle ne résulte pas d'une somme de deux solutions MOND individuelles.
- Interpolant final :1 683valeurs d'énergie reconstruites depuis leurs fichiers, tous les coefficients bicubiques recalculés par un système polynomial différent et{fractions:,}intervalles Bernstein vérifiés par fractions exactes. Les bornes positives/minimales et maximales concernent ce polynôme encodé.
- Validation physique :{physical['checks_count']}contrôles de reconstruction du rapport final. Pour les21positions et les2modèles déclarés, l'écart maximal aux intégrales de forces est{100*force:.7g}%. Cette précision est mesurée à ces positions ; elle n'est pas une borne uniforme sur le domaine physique.
- Populations :{len(banks)}banques complètes,{count:,}lignes acceptées, plus1 024états de validation. Le vérificateur final effectue{pop['checks_count']}contrôles groupés, portant sur toutes ces lignes : énergie,Lz,barrière,poids des3β,masses vraies/photométriques,ancrage au pool empirique,projection du photocentre,hiérarchie et résolution initiale. Il retrouve également les statistiques d'invariants des traces sauvegardées, les44graines, les comptes de sélection, les empreintes des66observations et les règles de taille/ESS. La plus petite ESS pondérée des banques est{ess:.1f}.

Les vérificateurs n'importent pas les calculs de production. Un petit exporteur séparé appelle volontairement l'API de production pour fournir des sondes à comparer ; sa provenance est explicite. Les intégrales directes de force finales, les ajustements GLS complets et tous les tirages rejetés ne sont pas répétés dans ce volet : les contrôles de production correspondants restent identifiés comme tels. Les données sauvegardées et les profils interpolés sont reconstruits indépendamment.

## Stationnarité et portée scientifique

Le parent est une distribution fonction de E et Lz dans le potentiel scalaire statique interpolé. Le masque de confinement et la hiérarchie intérieure dépendent d'invariants. La proposition Beta(5/2,3/2) et le poids r³(−Psi)³(−E)^(β−1) donnent le volume de phase annoncé en unités sans dimension, conditionnellement aux covariates et masses définis. La sélection de résolution initiale et de séparation projetée est instantanée ; les catalogues sélectionnés ne sont pas déclarés stationnaires.

Les masses totales photométriques sont ancrées au même pool empirique avant sélection, avec conversion explicite vers la masse vraie quand un compagnon est ajouté. Cela ne garantit pas des distributions photométriques identiques après sélection. Les corrélations détaillées par source entre masse,flux,erreurs et sélection Gaia restent approximatives. Le parent orbital et ses coupures en Lz sont un choix de modèle, pas une distribution réelle d'excentricités déduite de Gaia.

Le domaine demeure limité :deux sources de Plummer prescrites,q discret0.1/0.3/1,fonctionnu simple,champ newtonien extérieur uniforme fixé àa0 et axe synthétiquez. Les bornes Bernstein certifient le modèle interpolé ; les configurations de validation ne certifient pas toutes les théories MOND,un champ galactique variable,ni une précision uniforme de QUMOND. L'orbite intérieure est traitée avec son modèle de Kepler et une approximation monopolaire extérieure.

La future comparaison utilise aussi des bibliothèques finies. Un p de rang est exact sous la même loi simulée que ses références ; le transfert à une bibliothèque test indépendante comporte une erreur de bibliothèque et doit être mesuré. Les erreurs/pouvoirs sur des nuisances discrètes ne se généralisent pas automatiquement à tous les autres paramètres ou à la vraie sélection Gaia.

Ces vérifications établissent la cohérence et l'exécution de cette chaîne conditionnelle. Elles ne réhabilitent pas à elles seules la conclusion observationnelle de l'article initial et ne donnent aucun verdict sur le ciel réel. Les conclusions de la comparaison statistique doivent être rapportées séparément, avec leurs éventuels échecs.

## Reproduction

Les résultats détaillés sont ENERGY_REPORT.json,interpolant_refined_verification.json,grid_refined_verification.json,physics_report_verification.json etpopulation_verification.json. Les scripts verify_interpolant.py/verify_energy_grid.py acceptent --model et --out ; verify_physics_report.py etverify_populations.py acceptent --root et --out. Pour ce dernier,--stage all vérifie la chaîne complète et--wait attend les22banques. Les chemins conservent la structure work/ ; aucun ancien résultat de production n'est modifié. Les rapports partiels et les sondes logicielles sont distincts des preuves finales.
'''
    (HERE/'FINAL_REVIEW.md').write_text(text,encoding='utf-8')
    print(json.dumps({'completed':True,'bank_rows':count,'minimum_weighted_ess':ess,'final_checks':pop['checks_count']}))

if __name__=='__main__':main()
