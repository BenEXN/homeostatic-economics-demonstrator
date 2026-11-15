# Paquet reproductible du démonstrateur EH — V1.1

Auteur : Benjamin Albouy
Rattachement : programme de recherche « Économie homéostatique »
Version : V1.1 — novembre 2026

Ce paquet accompagne l'article « Mise en calcul d'une architecture biophysique–monétaire sous contraintes : un démonstrateur numérique exploratoire » (V1.1). Il contient le moteur de référence, ses sorties, un script de contrôle complémentaire et des tests. Il se suffit à lui-même : aucun autre fichier n'est nécessaire pour reproduire les tableaux et les figures de l'article.

## Contenu

| Élément | Rôle |
|---|---|
| `eh_reproductible.py` | Moteur de référence. Une exécution régénère les dix tableaux, les CSV et les figures de l'article. |
| `csv/` | Dix tableaux de référence, T1 à T10. |
| `figures/` | Huit figures de l'article, en PNG (300 dpi) et en PDF vectoriel. |
| `execution.log` | Sortie standard de l'exécution de référence. |
| `EXECUTION_REFERENCE.txt` | Environnement, commande et empreinte de l'exécution de référence ; historique des révisions. |
| `controle/` | Script de contrôle par intégrateurs d'ordre supérieur, ses dépendances et sa sortie attendue. |
| `tests/` | Tests du moteur et test de non-régression du paquet, avec leurs sorties attendues. |
| `CHANGELOG_V1.1.md` | Modifications de la V1 à la V1.1. |
| `MANIFEST_SHA256.txt` | Empreintes SHA-256 de tous les fichiers du paquet. |

Empreinte SHA-256 du moteur, identique à celle imprimée au §13 de l'article :

    ffc2339bcb2b61adb5a0c71f6d30f063636d1015a3db547e36e7bd62f89bd9b5

## Dépendances

Obligatoire : Python 3.8 ou ultérieur. Le moteur, ses tableaux et ses CSV reposent sur la bibliothèque standard (`math`, `csv`, `os`, `fractions`). Les tests reposent eux aussi sur la bibliothèque standard.

Facultatives :

- `matplotlib`, pour les figures. En son absence, le moteur produit les tableaux et les CSV et signale que les figures ne sont pas générées.
- `numpy` et `scipy`, pour le script de contrôle `controle/controle_ordre_superieur.py`.

Installation des dépendances facultatives :

    python3 -m pip install matplotlib
    python3 -m pip install -r controle/requirements.txt

Environnement de l'exécution de référence : Python 3.12.3, matplotlib 3.10.8 ; script de contrôle : numpy 2.5.4, scipy 1.18.1.

## Commandes et résultats attendus

Le moteur écrit ses sorties dans `csv/` et `figures/`, à côté de lui-même. Pour conserver les fichiers de référence, copier le moteur dans un répertoire vide avant de l'exécuter ; le test de non-régression le fait automatiquement.

**1. Vérifier l'intégrité du paquet**

    sha256sum -c MANIFEST_SHA256.txt

Attendu : `OK` pour chaque fichier.

**2. Régénérer les tableaux et les figures**

    mkdir essai && cp eh_reproductible.py essai/ && cd essai
    python3 eh_reproductible.py > execution.log

Attendu : code de retour 0 ; dix CSV dans `csv/` et le journal `execution.log` identiques octet pour octet à ceux du paquet ; seize fichiers dans `figures/`. Durée : quelques secondes.

**3. Test de non-régression (étapes 2 et comparaison, automatisées)**

    python3 tests/test_non_regression.py

Attendu : `BILAN : dix CSV et journal identiques octet pour octet`. Les figures sont identiques octet pour octet dans l'environnement de référence ; dans un autre environnement, leur rendu peut différer sans que les données changent, et le test le signale à titre d'information.

**4. Tests du moteur**

    python3 tests/tests_moteur.py

Attendu : `BILAN : 111 tests conformes, 0 non conformes` (sortie complète dans `tests/tests_moteur_sortie_attendue.txt`). Les tests couvrent la formule fermée du schéma d'Euler dans le scénario A, les métriques de la zone cible, le détecteur de retournements, les contrôles préalables (paramètres inadmissibles, valeurs non finies, horizon, divisibilité de T et de τ par Δt), l'indexation du retard, le certificat de l'instance B en arithmétique rationnelle exacte, les paramètres du scénario E à τ = 10 pour deux états initiaux, le régime du flux terminal, les contrôles de grille (au moins un intervalle, au moins un pas de retard) et le point fixe bouclé.

**5. Contrôle par intégrateurs d'ordre supérieur**

    python3 controle/controle_ordre_superieur.py

Attendu : `BILAN : tous les contrôles sont CONFORMES` (sortie complète dans `controle/sortie_attendue.txt`). Ce script intègre les mêmes équations hors du moteur de référence, par l'intégrateur DOP853 (rtol = 1e-12) et, pour τ > 0, par la méthode des pas. Ses paramètres sont retranscrits du §6 de l'article et comparés à ceux du moteur. Pour le scénario F, il contrôle les entrées et les sorties de la zone cible complète, contrainte biophysique et fenêtre monétaire vérifiées conjointement. Il fournit l'instant de franchissement IED = k_env de la solution continue du scénario C, cité au §8 C de l'article : t ≈ 40,048363. Durée : une dizaine de secondes.

## Reproductibilité

Le moteur ne comporte aucune source d'aléa. Deux exécutions successives dans un même environnement produisent des fichiers identiques octet pour octet. Entre environnements différents, les CSV et le journal dépendent de l'arithmétique flottante en double précision et le rendu des figures dépend de la version de matplotlib et des polices disponibles ; les données tracées sont les mêmes.

Le moteur exige que tous les paramètres soient des nombres finis et que l'horizon T et le retard τ soient des multiples entiers du pas Δt (quotients contrôlés à 1e-9 près), avec au moins un intervalle pour l'horizon et au moins un pas pour tout retard strictement positif ; il signale par une erreur explicite toute configuration qui ne satisfait pas ces conditions.

Le critère relatif de stabilisation du flux est défini pour q(T) > 0 ; le moteur rapporte séparément un flux nul constant (champ `flux_regime`). Le point fixe bouclé est renvoyé avec un résidu au plus égal à 1e-12, ou refusé par une erreur explicite.

## Usage d'outils d'intelligence artificielle

Claude Fable 5.2 (Anthropic) a été utilisé pour la préparation de la version V1.1 : contrôle de reproductibilité de la V1, application des corrections à l'article et au moteur, scripts de contrôle et de test, constitution du paquet de reproductibilité. L'auteur a vérifié les contenus produits et reste seul responsable du contenu.
