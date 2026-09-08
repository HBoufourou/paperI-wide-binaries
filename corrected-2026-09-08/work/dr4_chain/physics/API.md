# Contrat du potentiel d'interaction

`interaction_model.InteractionModel.load()` charge `results/interaction_model.npz` seulement si `validation.json` déclare la réussite et correspond aux empreintes du modèle et de l'évaluateur. L'argument interne `require_validation=False` sert au développement ; son emploi ne constitue pas une validation.

Méthodes : `potential(positions,q,gravity='qumond')`, alias `psi`, et `acceleration(...)`. Les positions ont une forme `(...,3)` ; q est scalaire ou diffusé sur les positions. Les résultats ont respectivement les formes `(...)` et `(...,3)`. L'axe z est celui du champ newtonien extérieur positif. Le vecteur position est R2−R1 avec m1=1/(1+q), m2=q/(1+q). Le signe de z est physiquement conservé.

Les modèles sont `'qumond'` et `'newton'`. Tous deux utilisent les mêmes sources de Plummer bi=0.00125√mi. Unités G=Mtotal=a0=1. Le potentiel relatif est Ψ=U/(m1m2) ; l'énergie relative spécifique vaut E=v²/2+Ψ. L'accélération retournée est a2−a1=−∇Ψ. L'origine est Ψ(∞)=0 ; le champ uniforme commun ne figure pas dans l'accélération relative.

Seuls q∈{0.1,0.3,1} et 0.03≤r≤30 sont supportés. Une tolérance numérique de10^-12 aux limites absorbe les arrondis ; elle ne permet pas une extrapolation scientifique. Les valeurs de q intermédiaires et les rayons hors domaine sont refusés. Aucune interpolation en q et aucun pliage μ→|μ| ne sont utilisés.

`cmax(q,gravity)`, `cmin(...)` et `outer_amax(...)` renvoient des bornes de A=−rΨ sur l'interpolant encodé. Le minimum et le maximum proviennent d'une conversion bicubique en base Bernstein avec arrondis vers l'extérieur. La borne extérieure s'applique uniquement à r=30. Ces bornes peuvent servir aux barrières en E,Lz du modèle interpolé ; elles ne certifient pas l'erreur uniforme du potentiel QUMOND physique.

L'archive NPZ expose les nœuds, valeurs, coefficients et bornes : `log_r`, `mu_nodes`, `q_nodes`, `gravity_names`, `A_values[gravity,q,r,mu]`, `coefficients[gravity,q,patch_r,patch_mu,power_r,power_mu]`, `bernstein_coefficients`, `bernstein_lower`, `bernstein_upper`, `cmax`, `cmin`, `outer_amax`, `metadata`.

Les puissances locales ne sont pas normalisées : (lnr−lnri)^p(μ−μj)^k. La base Bernstein, elle, utilise les coordonnées locales normalisées sur[0,1]². Les deux axes sont interpolés par des splines cubiques «not-a-knot». Potentiel et accélération dérivent d'un même scalaire ; il n'existe pas de grille de force interpolée séparément.

La validation des forces porte sur les positions déclarées dans le protocole. Le chargement réussi ne signifie ni précision uniforme garantie entre ces positions, ni modèle astrophysique complet, ni verdict sur Gaia.
