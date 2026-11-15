#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Auteur : Benjamin Albouy
# Rattachement : programme de recherche « Économie homéostatique »
"""
================================================================================
Démonstrateur EH V1.1 — TEST DE NON-RÉGRESSION DU PAQUET
================================================================================
Copie le moteur eh_reproductible.py dans un répertoire temporaire, l'exécute, et
compare ses sorties aux fichiers de référence du paquet :

  - les dix CSV et le journal execution.log : identité octet pour octet attendue ;
  - les figures : identité octet pour octet dans l'environnement de référence
    (Python 3.12.3, matplotlib 3.10.8) ; dans un autre environnement, le rendu
    graphique peut varier sans que les données soient modifiées, et l'écart est
    alors rapporté à titre d'information, sans échec du test.

Le paquet lui-même n'est pas modifié.

DÉPENDANCES : bibliothèque standard ; matplotlib facultatif (figures).
USAGE : python3 test_non_regression.py
Code de retour : 0 si les dix CSV et le journal sont identiques, 1 sinon.
================================================================================
"""
import sys, os, shutil, subprocess, tempfile, filecmp

here = os.path.dirname(os.path.abspath(__file__)); paquet = os.path.dirname(here)
tmp = tempfile.mkdtemp(prefix="eh_v11_")
try:
    shutil.copy(os.path.join(paquet, "eh_reproductible.py"), tmp)
    with open(os.path.join(tmp, "execution.log"), "w", encoding="utf-8") as log:
        rc = subprocess.run([sys.executable, "eh_reproductible.py"], cwd=tmp, stdout=log,
                            env=dict(os.environ, PYTHONIOENCODING="utf-8")).returncode
    print(f"Exécution du moteur : code de retour {rc}")
    ok = (rc == 0)
    noms = sorted(os.listdir(os.path.join(paquet, "csv")))
    for n in ["csv/" + x for x in noms] + ["execution.log"]:
        a, b = os.path.join(tmp, n), os.path.join(paquet, n)
        same = os.path.exists(a) and filecmp.cmp(a, b, shallow=False); ok &= same
        print(f"[{'IDENTIQUE' if same else 'DIFFÉRENT'}] {n}")
    figdir = os.path.join(tmp, "figures")
    if os.path.isdir(figdir) and os.listdir(figdir):
        for n in sorted(os.listdir(os.path.join(paquet, "figures"))):
            a, b = os.path.join(figdir, n), os.path.join(paquet, "figures", n)
            same = os.path.exists(a) and filecmp.cmp(a, b, shallow=False)
            print(f"[{'IDENTIQUE' if same else 'RENDU DIFFÉRENT (information)'}] figures/{n}")
    else:
        print("Figures non générées (matplotlib absent) : comparaison graphique omise.")
    print("BILAN : " + ("dix CSV et journal identiques octet pour octet" if ok else "ÉCART sur les CSV ou le journal"))
finally:
    shutil.rmtree(tmp, ignore_errors=True)
sys.exit(0 if ok else 1)
