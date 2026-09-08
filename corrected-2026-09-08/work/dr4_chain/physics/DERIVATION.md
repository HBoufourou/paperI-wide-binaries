# Scalaire utilisé par la chaîne

Les conventions sont G=Mtotal=a0=1, m1=1/(1+q), m2=q/(1+q), R1=−m2 r et R2=m1 r. Le champ newtonien extérieur uniforme vaut e_z ; son gradient de potentiel vaut u_e=−e_z. Pour la source i,

u_i(x)=m_i(x−R_i)/(|x−R_i|²+b_i²)^(3/2), avec b_i=0.00125√m_i.

La théorie QUMOND considérée possède ν(y)=(1+√(1+4/y))/2 et Q'(z)=ν(√z). On isole Q_ph(y)=Q(y²)−y², dont la dérivée est dQ_ph/dy=2y(ν−1). L'énergie d'interaction s'écrit

U=U_N−(1/8π)∫[Q_ph(|u_e+u_1+u_2|)−Q_ph(|u_e+u_1|)−Q_ph(|u_e+u_2|)+Q_ph(|u_e|)]d³x.

Le terme U_N est l'interaction exacte entre les deux distributions de Plummer, intégrée par une quadrature unidimensionnelle indépendante. Une somme de deux solutions QUMOND à une masse ne donnerait pas cette expression : la somme des champs newtoniens intervient **avant** l'application de Q_ph.

Pour éviter la perte de chiffres significatifs dans la partie extérieure, la différence mixte de Q_ph est évaluée à faible champ interne par

ΔQ_ph=2∫_0^1∫_0^1 u_1ᵀ D[(ν−1)u](u_e+s u_1+t u_2)u_2 ds dt.

Le Hessien vaut (ν−1)I+(ν'/|u|)uuᵀ. Une quadrature Gauss2×2 est utilisée uniquement lorsque |u1|+|u2|<0.01. Ailleurs l'inclusion/exclusion est évaluée directement avec une expression stable de Q_ph.

L'intégrale spatiale décroît seulement comme1/R. La queue fantôme principale est −m1m2[ν_e(1+K_e/3)−1]/R, K_e=dlnν/dlny à y=1. Le doublement de R teste les contributions résiduelles, notamment la frontière multicentre non exactement sphérique. L'origine physique Ψ(∞)=0 est définie grâce au champ extérieur non nul. Le cas isolé, dont l'énergie comporte une divergence logarithmique à l'infini, ne fait pas partie de cette API.

La masse réduite est μ_red=m1m2, et Ψ=U/μ_red. La dérivée de l'énergie à centre de masse fixé donne ∇_r U=m2F1−m1F2, donc −∇Ψ=F2/m2−F1/m1. Cette relation permet de comparer directement le gradient de l'interpolant aux anciennes intégrales de forces indépendantes sur chacune des sources, sans imposer numériquement leur réciprocité.

L'axisymétrie autour de e_z implique les invariants E=v²/2+Ψ et Lz=(r×v)_z. Elle n'implique pas une symétrie r_z→−r_z quand q≠1. Une population f(E,Lz) peut être stationnaire dans ce potentiel statique, mais une sélection instantanée en séparation projetée ne conserve généralement pas cette stationnarité. Le choix astrophysique de la fonction de distribution demeure une hypothèse séparée.

Références primaires : [action QUMOND, Milgrom2010](https://arxiv.org/html/0911.5464), [traitement complet des deux masses, Pflamm-Altenburg2025](https://arxiv.org/html/2509.01493v1), [limite extérieure dominante, Banik & Zhao](https://arxiv.org/html/1509.08457v3). Les sources physiques sont prescrites, étendues et rigides ; ces équations ne modélisent pas leur déformation ni un champ galactique variant le long de la trajectoire.
