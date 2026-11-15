# Modifications de la V1 à la V1.1

Auteur : Benjamin Albouy — novembre 2026

La V1.1 est une révision corrective. Les équations, les paramètres des scénarios et les résultats publiés sont ceux de la V1 : les dix CSV et le journal d'exécution produits par le moteur V1.1 sont identiques octet pour octet à ceux de la V1.

## Moteur `eh_reproductible.py`

| Objet | V1 | V1.1 |
|---|---|---|
| Valeurs non finies | non contrôlées | tout paramètre et le pas Δt sont des nombres finis ; erreur explicite sinon |
| Horizon | non contrôlé | T fini et strictement positif |
| Grille | N = round(T/Δt), τ_steps = round(τ/Δt) | T et τ multiples entiers de Δt (tolérance 1e-9 sur le quotient) ; N et τ_steps entiers exacts ; N ≥ 1 ; τ_steps ≥ 1 pour τ > 0 |
| Flux terminal | flux_rel = 0 lorsque q(T) = 0 | critère relatif pour q(T) > 0 ; flux nul constant et extinction rapportés par le champ `flux_regime` |
| Point fixe bouclé | sans changement de signe sur la grille : point de grille minimisant \|g\| | racine en un point de grille, ou encadrée sur [0, 1e-4] et obtenue par bissection ; résidu au plus égal à 1e-12, erreur explicite sinon |
| Détecteur de retournements | pas significatif si \|m_k\| > θ | pas significatif si \|m_k\| ≥ θ, comme dans l'annexe A.3 de l'article |
| Tableau T8 | proportion dans K de A et B écrite en constante | proportion lue dans les résultats de la simulation |
| Figure 6, panneau A | encart placé sur le point initial | encart placé dans une zone libre du panneau ; trajectoire et point initial entièrement visibles |
| En-tête | documents de référence datés d'octobre 2026 | novembre 2026 ; liste des modules (`fractions`) ; origine de l'instant de franchissement continu du scénario C |
| Empreinte SHA-256 | `7a81d614ba4efcb6bffafffd37e490fa2475f1a9811577676109619560d5d15b` | `ffc2339bcb2b61adb5a0c71f6d30f063636d1015a3db547e36e7bd62f89bd9b5` |

Les configurations publiées satisfont toutes les nouveaux contrôles.

## Sorties

- `csv/` (T1 à T10) et `execution.log` : identiques octet pour octet à la V1.
- Figures 1 à 5, S1 et S2 (PNG et PDF) : identiques octet pour octet à la V1, dans l'environnement de référence.
- Figure 6 : régénérée. Les données tracées sont les mêmes ; l'encart du panneau A est déplacé, et la mise en page automatique décale légèrement les panneaux C et F.

## Fichiers ajoutés

- `controle/controle_ordre_superieur.py`, `controle/requirements.txt`, `controle/sortie_attendue.txt` : contrôle par l'intégrateur DOP853 et par la méthode des pas, hors du moteur de référence ; paramètres retranscrits de l'article ; événements de la zone cible complète (contrainte biophysique et fenêtre monétaire).
- `tests/tests_moteur.py`, `tests/test_non_regression.py` et leurs sorties attendues.
- `README.md`, `CHANGELOG_V1.1.md`.

## Article

- Version et note de version sur la page de titre ; empreinte du moteur au §13.
- §4 et annexe A.6 : vecteurs écrits avec le point-virgule séparateur ; notation e_base.
- §5, §10, §11, §12, §13 : le contrôle par intégrateur d'ordre supérieur et par la méthode des pas est décrit comme un contrôle complémentaire réalisé hors du moteur de référence ; le §13 renvoie au script par son nom.
- §5 : tolérance des contrôles de domaine et d'admissibilité ; définition du flux stabilisé.
- Annexes A.2, A.5, A.6 : recherche du point fixe, domaine du critère relatif de flux, contrôles de grille et leur tolérance.
- §2.1, §5, §8 G, chapeau de l'annexe A : précisions de rédaction.
- §7 et §8 A : bornes inférieures de l'IED écrites 0,274 et 1,360.
- §8 C : l'instant de franchissement continu est attribué au script de contrôle.
- §8 E, annexe A.3 et légende de la figure S1 : l'état initial est nommé dans la portée de l'observation ; le critère du détecteur est précisé.
- Annexes renommées A, B, C ; signes « % » ; valeur de IRNR₀ dans l'annexe d'ancrage empirique (0,54).
- Figure 6 remplacée par la figure régénérée.
