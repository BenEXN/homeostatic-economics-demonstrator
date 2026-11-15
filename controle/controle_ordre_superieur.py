#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Auteur : Benjamin Albouy
# Rattachement : programme de recherche « Économie homéostatique »
"""
================================================================================
Démonstrateur EH V1.1 — CONTRÔLE COMPLÉMENTAIRE PAR INTÉGRATEURS D'ORDRE SUPÉRIEUR
================================================================================
Ce script est extérieur au moteur de référence (eh_reproductible.py). Il intègre
les mêmes équations par une méthode indépendante du schéma d'Euler :

  - τ = 0 : intégrateur DOP853 (Runge–Kutta explicite d'ordre 8, pas adaptatif,
    scipy.integrate.solve_ivp), rtol = 1e-12, atol = 1e-14 ;
  - τ > 0 : méthode des pas — intégration DOP853 sur les intervalles successifs
    [nτ, (n+1)τ], le terme retardé IED(t−τ) étant lu dans la solution dense de
    l'intervalle précédent (histoire constante IED_0 pour t−τ < 0).

Il produit :
  1. la comparaison, scénario par scénario, des trajectoires du moteur (Euler,
     Δt = 0,05 et 0,025) avec la solution d'ordre supérieur ;
  2. l'instant de franchissement IED = k_env de la solution continue du
     scénario C, cité au §8 C de l'article (valeur attendue : t ≈ 40,048363) ;
  3. pour le scénario E, l'écart maximal à la trajectoire sans retard ;
  4. pour chaque valeur de la grille β du scénario F (tableau T7), les entrées
     et les sorties de la zone cible complète Z_cible, contrainte biophysique et
     fenêtre monétaire vérifiées conjointement à l'instant de l'événement.

Les paramètres des scénarios sont retranscrits du §6 de l'article ; le script
contrôle qu'ils sont identiques à ceux du moteur avant de les employer.

Équations (article, §4) : dS_i/dt = −α_i S_i + β_i u + γ_i ; dM/dt = u − λM ;
u(t) = min(IED(t−τ), 1) ; IED = (S_1 S_2 S_3)^(1/3) ;
α⃗ = (α ; α ; 1,3α), β⃗ = (β ; β ; 0,5β), γ⃗ = (0,02 ; 0,018 ; 0,012).

DÉPENDANCES : numpy, scipy (voir requirements.txt). Le moteur de référence est
importé pour la comparaison ; il n'est ni modifié ni exécuté comme programme.
Les intégrations de ce script n'appellent aucune fonction du moteur.

USAGE : python3 controle_ordre_superieur.py [chemin/vers/eh_reproductible.py]
        (par défaut : ../eh_reproductible.py). Sortie attendue : sortie_attendue.txt
Code de retour : 0 si tous les contrôles sont conformes, 1 sinon.
================================================================================
"""
import sys, os, math, importlib.util
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

RTOL, ATOL = 1e-12, 1e-14
T_FINAL = 70.0
GAMMA = np.array([0.02, 0.018, 0.012])
TOL_EVT = 0.15                   # tolérance sur les instants d'événement de Z_cible (trois pas de grille)

# Paramètres retranscrits du §6 de l'article (tableaux des scénarios), indépendamment du moteur.
def _p(alpha, beta, lam, tau, s0, m0, k_env, k_social, z_cible, m_lo, m_hi):
    return {'alpha': alpha, 'beta': beta, 'lambda': lam, 'tau': tau, 'ibd0': s0[0], 'iee0': s0[1], 'irnr0': s0[2],
            'm0': m0, 'k_env': k_env, 'k_social': k_social, 'z_cible': z_cible, 'm_cible_lo': m_lo, 'm_cible_hi': m_hi}
PARAMS = {
    'A': _p(0.10, 0.18, 0.12, 0,  (1.5, 1.4, 1.2),   1.2, 0.40, 0.30, 0.90, 0.8, 1.5),
    'B': _p(0.07, 0.03, 0.15, 0,  (0.7, 0.7, 0.5),   0.8, 0.30, 0.20, 0.90, 0.7, 1.4),
    'C': _p(0.07, 0.02, 0.15, 0,  (0.7, 0.7, 0.5),   0.8, 0.30, 0.20, 0.90, 0.7, 1.4),
    'E': _p(0.22, 0.28, 0.15, 10, (0.45, 0.45, 0.4), 0.8, 0.20, 0.15, 1.20, 0.7, 1.4),
    'F': _p(0.12, 0.16, 0.82, 0,  (0.6, 0.55, 0.5),  1.0, 0.30, 0.25, 0.90, 0.7, 1.5),
}
BETA_GRILLE = [i / 100 for i in range(10, 23)]

T_C_ATTENDU = 40.048363          # franchissement continu du scénario C (§8 C)
TOL_T_C = 5e-6

here = os.path.dirname(os.path.abspath(__file__))
engine_path = sys.argv[1] if len(sys.argv) > 1 else os.path.join(here, "..", "eh_reproductible.py")
spec = importlib.util.spec_from_file_location("eh_reproductible", engine_path)
eh = importlib.util.module_from_spec(spec); spec.loader.exec_module(eh)


def vecteurs(p):
    a, b = p['alpha'], p['beta']
    return np.array([a, a, 1.3 * a]), np.array([b, b, 0.5 * b])

def ied(S):
    return float(np.prod(S)) ** (1.0 / 3.0)

def integrer_sans_retard(p, T=T_FINAL):
    """DOP853 sur [0, T], τ = 0. Renvoie une fonction t -> (S1, S2, S3, M)."""
    AL, BE = vecteurs(p); lam = p['lambda']
    def f(t, y):
        u = min(ied(y[:3]), 1.0)
        return np.concatenate((-AL * y[:3] + BE * u + GAMMA, [u - lam * y[3]]))
    y0 = [p['ibd0'], p['iee0'], p['irnr0'], p['m0']]
    sol = solve_ivp(f, (0.0, T), y0, method="DOP853", rtol=RTOL, atol=ATOL, dense_output=True)
    assert sol.success, sol.message
    return sol.sol

def integrer_avec_retard(p, T=T_FINAL):
    """Méthode des pas : DOP853 sur [nτ, (n+1)τ], terme retardé lu dans la
    solution dense de l'intervalle précédent. Renvoie une fonction t -> état."""
    tau = float(p['tau']); AL, BE = vecteurs(p); lam = p['lambda']
    y = np.array([p['ibd0'], p['iee0'], p['irnr0'], p['m0']], dtype=float)
    ied0 = ied(y[:3]); morceaux = []          # liste de (t_debut, t_fin, solution dense)
    t0 = 0.0
    while t0 < T - 1e-12:
        t1 = min(t0 + tau, T)
        precedent = morceaux[-1][2] if morceaux else None
        def f(t, z, precedent=precedent):
            ied_ret = ied0 if precedent is None else ied(precedent(t - tau)[:3])
            u = min(ied_ret, 1.0)
            return np.concatenate((-AL * z[:3] + BE * u + GAMMA, [u - lam * z[3]]))
        sol = solve_ivp(f, (t0, t1), y, method="DOP853", rtol=RTOL, atol=ATOL, dense_output=True)
        assert sol.success, sol.message
        morceaux.append((t0, t1, sol.sol)); y = sol.y[:, -1]; t0 = t1
    def etat(t):
        for a, b, s in morceaux:
            if a - 1e-12 <= t <= b + 1e-12:
                return s(min(max(t, a), b))
        raise ValueError(t)
    return etat

def solution(p):
    return integrer_avec_retard(p) if p['tau'] > 0 else integrer_sans_retard(p)

def comparer(nom, p):
    """Écart maximal sur la grille entre Euler (moteur) et la solution d'ordre supérieur."""
    ref = solution(p); lignes = []
    for dt in (0.05, 0.025):
        r = eh.run(p, dt)
        e_ied = max(abs(w['IED'] - ied(ref(w['t'])[:3])) for w in r['rows'])
        e_m = max(abs(w['M'] - ref(w['t'])[3]) for w in r['rows'])
        lignes.append((dt, r['ied_final'], e_ied, e_m))
    yT = ref(T_FINAL)
    print(f"  {nom:<14} ordre supérieur : IED(70) = {ied(yT[:3]):.10f} ; M(70) = {yT[3]:.8f}")
    for dt, iedf, e_ied, e_m in lignes:
        print(f"      Euler Δt={dt:<6} IED(70) = {iedf:.10f} ; écart max sur la grille : IED {e_ied:.2e} ; M {e_m:.2e}")
    rapport = lignes[0][2] / lignes[1][2]
    print(f"      rapport des écarts max sur l'IED entre Δt = 0,05 et 0,025 : {rapport:.2f}")
    return ref, rapport


def evenements_zcible(ref, p, T=T_FINAL, h=0.005):
    """Entrées et sorties de Z_cible pour une solution dense ref : t -> (S1, S2, S3, M).
    Renvoie (occupation initiale, [(instant d'entrée, contrainte active)], [(instant de sortie, contrainte active)]).
    Une excursion plus courte que h entre deux instants d'évaluation n'est pas détectée."""
    def marges(t):
        y = ref(t); return (ied(y[:3]) - p['z_cible'], y[3] - p['m_cible_lo'], p['m_cible_hi'] - y[3])
    c = lambda t: min(marges(t))
    noms = ("IED = z_cible", "M = m_lo", "M = m_hi")
    n = int(round(T / h)); ts = [min(i * h, T) for i in range(n + 1)]
    dedans = [c(t) >= 0 for t in ts]; entrees, sorties = [], []
    for i in range(1, n + 1):
        if dedans[i] != dedans[i - 1]:
            t0 = brentq(c, ts[i - 1], ts[i], xtol=1e-12, rtol=1e-14)
            m = marges(t0); active = noms[min(range(3), key=lambda k: abs(m[k]))]
            (entrees if dedans[i] else sorties).append((t0, active))
    return dedans[0], entrees, sorties

def controler_zcible(nom, p):
    """Compare les événements de Z_cible de la solution d'ordre supérieur à ceux du moteur."""
    ref = solution(p); init, entrees, sorties = evenements_zcible(ref, p)
    r = eh.run(p)
    sortie_finale = sorties[-1][0] if (sorties and (not entrees or sorties[-1][0] > entrees[-1][0])) else None
    derniere_sortie = sorties[-1][0] if sorties else None
    accord = (init == r['init_in_z']) and (len(entrees) == r['z_n_entries'])
    if entrees:
        accord &= r['z_first_entry_t'] is not None and abs(r['z_first_entry_t'] - entrees[0][0]) <= TOL_EVT
    accord &= (derniere_sortie is None) == (r['z_last_exit_t'] is None)
    if derniere_sortie is not None and r['z_last_exit_t'] is not None:
        accord &= abs(r['z_last_exit_t'] - derniere_sortie) <= TOL_EVT
    fmt = lambda ev: " ; ".join(f"t = {t:.4f} ({a})" for t, a in ev) if ev else "aucune"
    print(f"  {nom}")
    print(f"      ordre supérieur : occupation initiale = {init} ; entrées : {fmt(entrees)} ; sorties : {fmt(sorties)}")
    print(f"      moteur          : occupation initiale = {r['init_in_z']} ; entrées : {r['z_n_entries']}"
          f" (première : {r['z_first_entry_t']}) ; dernière sortie : {r['z_last_exit_t']} -> {'CONFORME' if accord else 'NON CONFORME'}")
    return accord


if __name__ == "__main__":
    import scipy
    ok = True
    print("=" * 96)
    print("CONTRÔLE COMPLÉMENTAIRE PAR INTÉGRATEURS D'ORDRE SUPÉRIEUR — démonstrateur EH V1.1")
    print(f"DOP853, rtol = {RTOL:g}, atol = {ATOL:g} ; méthode des pas pour τ > 0")
    print("=" * 96)

    SC = PARAMS
    print("\n0. Paramètres retranscrits du §6 de l'article, comparés à ceux du moteur")
    ident = all(PARAMS[k[0]] == dict(v) for k, v in eh.SCN.items()) and len(eh.SCN) == len(PARAMS) \
            and BETA_GRILLE == list(eh.BETA_GRID) and tuple(GAMMA) == tuple(eh.GAMMA)
    print(f"  cinq scénarios, grille β et vecteur γ identiques -> {'CONFORME' if ident else 'NON CONFORME'}")
    ok &= ident
    print("\n1. Moteur de référence (Euler) contre solution d'ordre supérieur, t ∈ [0, 70]")
    refs = {}
    for nom, p in (("A", SC['A']), ("B", SC['B']), ("C", SC['C']), ("F (β=0,16)", SC['F'])):
        refs[nom], rapport = comparer(nom, p)
        ok &= 1.8 < rapport < 2.2
    for tau in (0, 3, 6, 10, 15):
        refs[f"E{tau}"], rapport = comparer(f"E (τ={tau})", dict(SC['E'], tau=tau))
        ok &= 1.8 < rapport < 2.2
    print("  Contrôle : les écarts sont divisés par deux environ lorsque le pas est divisé par deux"
          f" -> {'CONFORME' if ok else 'NON CONFORME'}")

    print("\n2. Scénario C — franchissement IED = k_env de la solution continue")
    pC = SC['C']; refC = refs["C"]
    g = lambda t: ied(refC(t)[:3]) - pC['k_env']
    tC = brentq(g, 30.0, 50.0, xtol=1e-12, rtol=1e-14)
    rC = eh.run(pC)
    print(f"  solution d'ordre supérieur (DOP853) : t = {tC:.6f}")
    print(f"  moteur de référence, premier état de la grille hors de K (Δt = 0,05) : t = {rC['exit_t']}")
    c2 = abs(tC - T_C_ATTENDU) < TOL_T_C
    print(f"  valeur attendue : t ≈ {T_C_ATTENDU:.6f} (tolérance {TOL_T_C:g}) -> {'CONFORME' if c2 else 'NON CONFORME'}")
    ok &= c2

    print("\n3. Scénario E — écart maximal à la trajectoire sans retard, sur la grille t_k = 0,05·k")
    grille = [round(0.05 * k, 6) for k in range(1401)]
    base = [ied(refs["E0"](t)[:3]) for t in grille]
    gaps = eh.delay_gaps(0.05)
    for tau in (3, 6, 10, 15):
        ecart = max(abs(ied(refs[f"E{tau}"](t)[:3]) - b) for t, b in zip(grille, base))
        print(f"  τ = {tau:<3} ordre supérieur : {ecart:.5f} ; moteur (Euler, Δt = 0,05, tableau T6) : {gaps[tau][0]:.5f}")
        ok &= abs(ecart - gaps[tau][0]) < 5e-4

    print("\n4. Zone cible complète Z_cible = {IED ≥ z_cible} ∩ {m_lo ≤ M ≤ m_hi} : entrées et sorties, t ∈ [0, 70]")
    print("   Événements de la solution d'ordre supérieur : appartenance évaluée tous les 0,005, instants affinés")
    print("   par recherche de racine sur c(t) = min(IED − z_cible, M − m_lo, m_hi − M) ; contrainte active indiquée.")
    print("   Accord exigé avec le moteur (grille, Δt = 0,05) : même occupation initiale, même nombre d'entrées,")
    print(f"   même présence d'une sortie finale, instants à moins de {TOL_EVT:g} près.")
    pF = PARAMS['F']
    cas = [(f"F, β = {beta:.2f}", dict(pF, beta=beta)) for beta in BETA_GRILLE]
    cas.append(("A (occupation initiale, puis sortie)", PARAMS['A']))
    cas.append(("contrôle : F, β = 0,16, m_hi = 1,15 (entrée puis sortie)", dict(pF, m_cible_hi=1.15)))
    for nom, p in cas:
        ok &= controler_zcible(nom, p)
    print("   Les deux derniers cas vérifient que le contrôle détecte une sortie par la fenêtre monétaire ;")
    print("   le dernier est une configuration de test, hors des scénarios publiés.")

    print("\nEnvironnement : Python %s ; numpy %s ; scipy %s" % (sys.version.split()[0], np.__version__, scipy.__version__))
    print("BILAN : " + ("tous les contrôles sont CONFORMES" if ok else "au moins un contrôle est NON CONFORME"))
    sys.exit(0 if ok else 1)
