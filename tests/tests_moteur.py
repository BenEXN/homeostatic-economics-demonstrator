#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Auteur : Benjamin Albouy
# Rattachement : programme de recherche « Économie homéostatique »
"""
================================================================================
Démonstrateur EH V1.1 — TESTS DU MOTEUR DE RÉFÉRENCE
================================================================================
Tests unitaires du moteur eh_reproductible.py, importé comme module (son bloc
principal n'est pas exécuté, aucun fichier n'est écrit). Chaque test compare un
résultat OBTENU à un résultat ATTENDU connu indépendamment du moteur (formule
fermée, suite construite à la main, arithmétique rationnelle exacte).

  1. Scénario A : formule fermée du schéma d'Euler en régime saturé (u = 1).
  2. Métriques de Z_cible sur des suites construites à la main.
  3. Détecteur de retournements sur des trajectoires analytiques.
  4. Contrôles préalables : paramètres inadmissibles, valeurs non finies,
     horizon, divisibilité de T et de τ par Δt.
  5. Indexation du retard : u_k = min(IED_{k−τ_steps}, 1), histoire IED_0.
  6. Certificat de l'instance B en arithmétique rationnelle exacte.
  7. Paramètres du scénario E, τ = 10 : état initial de l'article et état
     initial (0,1 ; 1,9 ; 1,9).
  8. Régime du flux terminal : critère relatif, flux nul constant, extinction.
  9. Grille : au moins un intervalle pour l'horizon, au moins un pas pour tout
     retard strictement positif.
 10. Point fixe bouclé : racine inférieure à la borne basse du balayage,
     contrôle du résidu, refus d'un point fixe non encadré.

DÉPENDANCES : bibliothèque standard uniquement.
USAGE : python3 tests_moteur.py [chemin/vers/eh_reproductible.py]
        (par défaut : ../eh_reproductible.py). Sortie attendue : tests_moteur_sortie_attendue.txt
Code de retour : 0 si tous les tests sont conformes, 1 sinon.
================================================================================
"""
import sys, os, math, importlib.util
from fractions import Fraction as Fr

here = os.path.dirname(os.path.abspath(__file__))
engine_path = sys.argv[1] if len(sys.argv) > 1 else os.path.join(here, "..", "eh_reproductible.py")
spec = importlib.util.spec_from_file_location("eh_reproductible", engine_path)
eh = importlib.util.module_from_spec(spec); spec.loader.exec_module(eh)
SC = {k[0]: dict(v) for k, v in eh.SCN.items()}

BILAN = {"CONFORME": 0, "NON CONFORME": 0}
def test(ok, libelle, detail=""):
    v = "CONFORME" if ok else "NON CONFORME"; BILAN[v] += 1
    print(f"[{v}] {libelle}" + (f"\n      {detail}" if detail else ""))
def refus(fn):
    """Renvoie le message d'erreur si fn() est refusé par le moteur, None sinon."""
    try: fn(); return None
    except (ValueError, RuntimeError) as e: return f"{type(e).__name__} : {str(e)[:80]}"

# ───────── 1. Scénario A : formule fermée du schéma d'Euler ─────────
print("== 1. Scénario A — formule fermée du schéma d'Euler (u = 1 sur tout l'horizon)")
pA = SC['A']; a, b, lam, dt = pA['alpha'], pA['beta'], pA['lambda'], 0.05
AL = [a, a, 1.3 * a]; BE = [b, b, 0.5 * b]; G = eh.GAMMA
S0 = [pA['ibd0'], pA['iee0'], pA['irnr0']]; Seq = [(BE[i] + G[i]) / AL[i] for i in range(3)]
r = eh.run(pA)
test(all(w['u'] == 1.0 for w in r['rows']), "A : coefficient de régulation saturé à chaque état", f"{len(r['rows'])} états")
e_ied = e_m = 0.0
for k, w in enumerate(r['rows']):
    S = [Seq[i] + (S0[i] - Seq[i]) * (1 - AL[i] * dt) ** k for i in range(3)]
    M = 1 / lam + (pA['m0'] - 1 / lam) * (1 - lam * dt) ** k
    e_ied = max(e_ied, abs(w['IED'] - (S[0] * S[1] * S[2]) ** (1 / 3))); e_m = max(e_m, abs(w['M'] - M))
test(e_ied < 1e-12 and e_m < 1e-12, "A : moteur contre S_k = S* + (S_0 − S*)(1 − αΔt)^k et M_k = 1/λ + (M_0 − 1/λ)(1 − λΔt)^k",
     f"écart maximal sur 1 401 états : IED {e_ied:.1e} ; M {e_m:.1e} (tolérance 1e-12)")
test(len(r['rows']) == 1401 and r['rows'][0]['t'] == 0.0 and r['rows'][-1]['t'] == 70.0, "N + 1 = 1 401 états, t_0 = 0, t_N = 70")

# ───────── 2. Métriques de Z_cible ─────────
print("== 2. Métriques de Z_cible sur des suites construites à la main (Δt = 0,5)")
def zm(seq): return eh.zcible_metrics([bool(x) for x in seq], [0.5 * k for k in range(len(seq))], 0.5)
for lab, seq, att in [
    ("jamais dans la zone", [0, 0, 0, 0], dict(n_entries=0, first_entry_t=None, last_exit_t=None, longest_stay=0.0)),
    ("occupation initiale sans entrée dynamique", [1, 1, 1, 1], dict(n_entries=0, first_entry_t=None, last_exit_t=None, longest_stay=1.5)),
    ("occupation initiale puis sortie", [1, 1, 0, 0], dict(n_entries=0, first_entry_t=None, last_exit_t=1.0, longest_stay=0.5)),
    ("une entrée, pas de sortie", [0, 0, 1, 1, 1], dict(n_entries=1, first_entry_t=1.0, last_exit_t=None, longest_stay=1.0)),
    ("deux séjours, le second plus long", [0, 1, 0, 1, 1, 1, 0], dict(n_entries=2, first_entry_t=0.5, last_exit_t=3.0, longest_stay=1.0)),
    ("état isolé dans la zone", [0, 0, 1, 0], dict(n_entries=1, first_entry_t=1.0, last_exit_t=1.5, longest_stay=0.0)),
    ("entrée au dernier état", [0, 0, 0, 1], dict(n_entries=1, first_entry_t=1.5, last_exit_t=None, longest_stay=0.0))]:
    obt = zm(seq); test(obt == att, f"Z_cible : {lab}", f"attendu {att} ; obtenu {obt}")

# ───────── 3. Détecteur de retournements ─────────
print("== 3. Détecteur de retournements sur des trajectoires analytiques (Δt = 0,05, T = 70)")
def osc(f, dt=0.05, T=70.0): return eh.detect_oscillations([f(k * dt) for k in range(round(T / dt) + 1)], dt)
for lab, f, natt, oatt in [
    ("constante", lambda t: 1.0, 0, False),
    ("rampe croissante", lambda t: 0.5 + 0.005 * t, 0, False),
    ("relaxation exponentielle", lambda t: 1 - 0.5 * math.exp(-0.1 * t), 0, False),
    ("un seul maximum (parabole)", lambda t: 1 - 0.3 * ((t - 35) / 35) ** 2, 1, False),
    ("sinus d'amplitude 0,1 et de période 20 (7 extremums)", lambda t: 1 + 0.1 * math.sin(2 * math.pi * t / 20), 7, True)]:
    n, o = osc(f); test(n == natt and o == oatt, f"détecteur : {lab}", f"attendu {natt} retournement(s), oscillation = {oatt} ; obtenu {n}, {o}")
th = eh.OSC_THETA_RATE; n, o = eh.detect_oscillations([0.0, th * 0.5, 0.0, th * 0.5, 0.0], 0.5)
test(n == 3 and o, "détecteur : pentes de valeur absolue exactement égale à θ (pas significatif si |m_k| ≥ θ)", f"attendu 3 retournements ; obtenu {n}")

# ───────── 4. Contrôles préalables ─────────
print("== 4. Contrôles préalables (attendu : refus avec message explicite)")
def avec(base='A', **kw):
    p = dict(SC[base]); p.update(kw); return p
nan, inf = float('nan'), float('inf')
cas_refus = [
    ("α = 0", lambda: eh.run(avec(alpha=0.0))), ("α < 0", lambda: eh.run(avec(alpha=-0.1))), ("β < 0", lambda: eh.run(avec(beta=-0.01))),
    ("λ = 0", lambda: eh.run(avec(**{'lambda': 0.0}))), ("τ < 0", lambda: eh.run(avec(tau=-1))),
    ("β + γ_1 > 2α (β = 0,19 pour α = 0,10)", lambda: eh.run(avec(beta=0.19))),
    ("S_0 hors de [0,2]³ (IBD_0 = 2,1)", lambda: eh.run(avec(ibd0=2.1))), ("S_0 négatif", lambda: eh.run(avec(iee0=-0.1))), ("M_0 < 0", lambda: eh.run(avec(m0=-0.5))),
    ("α_3 Δt > 1 (Δt = 8)", lambda: eh.run(avec(), 8.0)), ("λ Δt > 1 (λ = 0,82 ; Δt = 1,3)", lambda: eh.run(avec('F'), 1.3)),
    ("Δt = 0", lambda: eh.run(avec(), 0.0)), ("Δt < 0", lambda: eh.run(avec(), -0.05)), ("Δt non fini", lambda: eh.run(avec(), nan)),
    ("T = 0", lambda: eh.run(avec(), 0.05, 0.0)), ("T < 0", lambda: eh.run(avec(), 0.05, -5.0)), ("T non fini", lambda: eh.run(avec(), 0.05, inf)),
    ("T/Δt non entier (Δt = 0,03 ; T = 70)", lambda: eh.run(avec(), 0.03)),
    ("τ/Δt non entier (τ = 0,125 ; Δt = 0,05)", lambda: eh.run(avec('E', tau=0.125))),
    ("τ/Δt non entier (τ = 10 ; Δt = 0,07)", lambda: eh.run(avec('E'), 0.07))]
for nom in ('alpha', 'beta', 'lambda', 'tau', 'ibd0', 'iee0', 'irnr0', 'm0', 'k_env', 'k_social', 'z_cible', 'm_cible_lo', 'm_cible_hi'):
    cas_refus.append((f"{nom} = NaN", lambda nom=nom: eh.run(avec(**{nom: nan}))))
    cas_refus.append((f"{nom} = +∞", lambda nom=nom: eh.run(avec(**{nom: inf}))))
for lab, fn in cas_refus:
    e = refus(fn); test(e is not None, f"refus : {lab}", e or "ACCEPTÉ sans erreur")
print("   (attendu : acceptation)")
for lab, fn in [("τ = 10, Δt = 0,025", lambda: eh.run(avec('E'), 0.025)), ("T = 3000, Δt = 0,05", lambda: eh.run(avec(), 0.05, 3000.0)),
                ("égalité d'invariance β + γ_1 = 2α (F, β = 0,22)", lambda: eh.run(avec('F', beta=0.22))),
                ("S_0 = (2 ; 2 ; 2)", lambda: eh.run(avec(ibd0=2.0, iee0=2.0, irnr0=2.0))), ("S_0 = (0 ; 0 ; 0)", lambda: eh.run(avec(ibd0=0.0, iee0=0.0, irnr0=0.0)))]:
    e = refus(fn); test(e is None, f"acceptation : {lab}", e or "exécution sans erreur")

# ───────── 5. Indexation du retard ─────────
print("== 5. Retard : u_k = min(IED_{k−τ_steps}, 1), histoire constante IED_0")
for tau in (3, 6, 10, 15):
    for dt in (0.05, 0.025):
        r = eh.run(avec('E', tau=tau), dt); ts = round(tau / dt); rows = r['rows']
        ok = all(w['u'] == min(rows[k - ts]['IED'] if k >= ts else r['ied0'], 1.0) for k, w in enumerate(rows))
        test(ok, f"τ = {tau}, Δt = {dt} : égalité exacte sur {len(rows)} états (τ_steps = {ts})")

# ───────── 6. Certificat de l'instance B ─────────
print("== 6. Certificat de l'instance B — arithmétique rationnelle exacte")
pB = SC['B']; al = [Fr(7, 100), Fr(7, 100), Fr(91, 1000)]; be = [Fr(3, 100), Fr(3, 100), Fr(15, 1000)]; ga = [Fr(2, 100), Fr(18, 1000), Fr(12, 1000)]
test((pB['alpha'], pB['beta'], pB['lambda']) == (0.07, 0.03, 0.15), "paramètres du scénario B : α = 0,07 ; β = 0,03 ; λ = 0,15")
def mul(p, q):
    out = [Fr(0)] * (len(p) + len(q) - 1)
    for i, x in enumerate(p):
        for j, y in enumerate(q): out[i + j] += x * y
    return out
prod = [Fr(1)]
for i in range(3): prod = mul(prod, [ga[i] / al[i], be[i] / al[i]])      # coefficients par degré croissant
prod[3] -= 1; Q = [c * 22295 for c in prod]
test(Q == [216, 954, 1395, -21620], "∏(β_i z + γ_i)/α_i − z³ = Q(z)/22295, Q(z) = −21620z³ + 1395z² + 954z + 216", f"coefficients obtenus (degré croissant) : {[int(c) for c in Q]}")
Qv = lambda z: sum(c * z ** k for k, c in enumerate(Q))
q1, q2 = Qv(Fr(31, 100)), Qv(Fr(311, 1000))
test(q1 == Fr(5369, 3125) and q2 == Fr(-135739961, 50000000), "Q(31/100) = 5369/3125 > 0 et Q(311/1000) = −135739961/50000000 < 0", f"obtenu {q1} et {q2}")
ie, res, Ss = eh.solve_fixed_point(pB['alpha'], pB['beta'])
test(0.310 < ie < 0.311 and res < 1e-12 and all(0 < s < 2 for s in Ss), "point fixe du moteur dans (0,310 ; 0,311), S* dans (0, 2)³", f"IED* = {ie:.10f} ; résidu {res:.1e}")
Ms = ie / pB['lambda']
test(31 / 15 < Ms < 311 / 150 and ie > pB['k_env'] and Ms > pB['k_social'], "M* = z*/λ dans (31/15 ; 311/150) ; (S*, M*) dans K", f"M* = {Ms:.6f}")
pe = dict(pB, ibd0=Ss[0], iee0=Ss[1], irnr0=Ss[2], m0=Ms); r = eh.run(pe)
test(r['pct_K'] == 100.0 and r['exit_t'] is None and max(r['ied_traj']) - min(r['ied_traj']) < 1e-12, "trajectoire issue de (S*, M*) : constante et dans K sur [0, 70]",
     f"dispersion de l'IED : {max(r['ied_traj']) - min(r['ied_traj']):.1e}")

# ───────── 7. Scénario E, τ = 10 : deux états initiaux ─────────
print("== 7. Paramètres du scénario E (α = 0,22 ; β = 0,28 ; λ = 0,15), τ = 10")
for dt in (0.05, 0.025):
    r = eh.run(avec('E'), dt)
    test(r['n_reversals'] == 0 and not r['oscillation'], f"état initial de l'article (0,45 ; 0,45 ; 0,4), Δt = {dt} : 0 retournement")
p2 = avec('E', ibd0=0.1, iee0=1.9, irnr0=1.9)
for dt in (0.05, 0.025):
    r = eh.run(p2, dt); x = r['ied_traj']
    sl = [(x[k + 1] - x[k]) / dt for k in range(len(x) - 1)]
    sg = [(k, 1 if m >= th else -1) for k, m in enumerate(sl) if abs(m) >= th]
    ext = [x[sg[i + 1][0]] for i in range(len(sg) - 1) if sg[i][1] != sg[i + 1][1]]       # valeurs de l'IED aux retournements
    amp = [abs(ext[i + 1] - ext[i]) for i in range(len(ext) - 1)]
    test(r['n_reversals'] == 10 and r['oscillation'], f"état initial (0,1 ; 1,9 ; 1,9), Δt = {dt} : 10 retournements relevés par le détecteur", f"obtenu {r['n_reversals']}")
    test(all(amp[i + 2] < amp[i] for i in range(len(amp) - 2)) and amp[-1] < 0.01 * amp[0],
         f"   amortissement : les alternances de même sens ont une amplitude strictement décroissante (Δt = {dt})",
         "amplitudes : " + " ; ".join(f"{v:.4f}" for v in amp))

# ───────── 8. Régime du flux terminal ─────────
print("== 8. Régime du flux terminal (critère relatif défini pour q(T) > 0)")
for lab, qT, q90, att in [
    ("q(T) = 1, q(0,9T) = 1 : critère relatif satisfait", 1.0, 1.0, ('relatif', 0.0, True)),
    ("q(T) = 0,5, q(0,9T) = 0,4995 : variation relative 1e-3, non strictement inférieure au seuil", 0.5, 0.4995, ('relatif', None, False)),
    ("q(T) = 0,5, q(0,9T) = 0,4 : critère relatif non satisfait", 0.5, 0.4, ('relatif', None, False)),
    ("q(T) = 0, q(0,9T) = 0 : flux nul constant, critère relatif non défini", 0.0, 0.0, ('nul_constant', None, False)),
    ("q(T) = 0, q(0,9T) = 0,3 : extinction", 0.0, 0.3, ('extinction', float('inf'), False))]:
    reg, rel, stab = eh.flux_terminal(qT, q90)
    ok = (reg == att[0] and stab == att[2] and (att[1] is None or rel == att[1]) and (reg != 'nul_constant' or rel is None))
    test(ok, f"flux_terminal : {lab}", f"obtenu régime = {reg}, flux_rel = {rel}, flux_stab = {stab}")
r = eh.run(avec('A', ibd0=0.0, iee0=0.0, irnr0=0.0, tau=100))
test(r['q_terminal'] == 0.0 and r['flux_regime'] == 'nul_constant' and r['flux_rel'] is None and r['flux_stab'] is False,
     "run : S_0 = (0 ; 0 ; 0), τ = 100 > T : u = 0 sur tout l'horizon, flux nul constant",
     f"q(T) = {r['q_terminal']} ; régime = {r['flux_regime']} ; flux_rel = {r['flux_rel']} ; flux_stab = {r['flux_stab']}")
for k in ('A', 'B', 'C', 'E', 'F'):
    r = eh.run(avec(k)); N = len(r['rows']) - 1; qT = r['rows'][-1]['u']; q90 = r['rows'][int(0.9 * N)]['u']
    test(r['flux_regime'] == 'relatif' and qT > 0 and r['flux_rel'] == abs(qT - q90) / qT and r['flux_stab'] == (r['flux_rel'] < 1e-3),
         f"scénario {k} : régime relatif, flux_rel = |q(T) − q(0,9T)|/q(T)", f"flux_rel = {r['flux_rel']:.3e} ; flux_stab = {r['flux_stab']}")

# ───────── 9. Grille : N >= 1 et τ_steps >= 1 ─────────
print("== 9. Grille : au moins un intervalle, au moins un pas de retard pour τ > 0")
for lab, fn in [("T = 1e-12, Δt = 0,05 (N = 0)", lambda: eh.run(avec(), 0.05, 1e-12)),
                ("T = 4e-11, Δt = 0,05 (N = 0)", lambda: eh.run(avec(), 0.05, 4e-11)),
                ("τ = 1e-12, Δt = 0,05 (τ_steps = 0)", lambda: eh.run(avec('E', tau=1e-12))),
                ("τ = 4e-11, Δt = 0,05 (τ_steps = 0)", lambda: eh.run(avec('E', tau=4e-11)))]:
    e = refus(fn); test(e is not None, f"refus : {lab}", e or "ACCEPTÉ sans erreur")
r = eh.run(avec(), 0.05, 0.05)
test(len(r['rows']) == 2 and r['rows'][-1]['t'] == 0.05, "acceptation : T = Δt = 0,05 : un intervalle, deux états, t_N = T")
r = eh.run(avec('E', tau=0.05)); r0 = eh.run(avec('E', tau=0))
test(r['ied_traj'] != r0['ied_traj'] and r['rows'][1]['u'] == min(r['ied0'], 1.0), "acceptation : τ = Δt = 0,05 : un pas de retard effectif")

# ───────── 10. Point fixe bouclé ─────────
print("== 10. Point fixe bouclé : racine sous la borne basse du balayage, résidu, encadrement")
for a in (1000.0, 300.0, 150.0):
    exact = (0.02 * 0.018 * 0.012 / (a * a * 1.3 * a)) ** (1.0 / 3.0)      # β = 0 : G_fp constante
    ie, res, Ss = eh.solve_fixed_point(a, 0.0)
    test(abs(ie - exact) <= 1e-9 * exact and res <= 1e-12, f"α = {a:g}, β = 0 : racine exacte (∏ γ_i/α_i)^(1/3) = {exact:.6e}",
         f"obtenu {ie:.6e} ; écart relatif {abs(ie - exact) / exact:.1e} ; résidu {res:.1e}")
import random
rnd = random.Random(20261113); pire_res = 0.0; pire_ecart = 0.0; n_cas = 0
for _ in range(2000):
    a = 10 ** rnd.uniform(-2, 3.2); b = rnd.uniform(0.0, max(0.0, 2 * a - 0.02))      # β + γ_1 ≤ 2α
    if b + 0.02 > 2 * a: continue
    ie, res, Ss = eh.solve_fixed_point(a, b); n_cas += 1
    AL = [a, a, 1.3 * a]; BE = [b, b, 0.5 * b]
    h = lambda z: ((BE[0] * min(z, 1) + G[0]) / AL[0] * (BE[1] * min(z, 1) + G[1]) / AL[1] * (BE[2] * min(z, 1) + G[2]) / AL[2]) ** (1 / 3) / z - 1
    lo_, hi_ = 1e-300, 2.5                                   # bissection indépendante sur H(z)/z − 1
    for _i in range(300):
        mid_ = 0.5 * (lo_ + hi_)
        if h(mid_) > 0: lo_ = mid_
        else: hi_ = mid_
    z_ = 0.5 * (lo_ + hi_)
    pire_res = max(pire_res, res); pire_ecart = max(pire_ecart, abs(ie - z_) / z_)
test(n_cas >= 1500 and pire_res <= 1e-12 and pire_ecart <= 1e-9, f"{n_cas} couples (α, β) admissibles tirés au hasard, α de 0,01 à 1 585 : accord avec une bissection indépendante",
     f"résidu maximal {pire_res:.1e} ; écart relatif maximal {pire_ecart:.1e}")
e = refus(lambda: eh.solve_fixed_point(0.1, 1.0))
test(e is not None, "refus : α = 0,1, β = 1 (paramètres inadmissibles, point fixe hors de [0, 2])", e or "valeur renvoyée sans erreur")
e = refus(lambda: eh.solve_fixed_point(float('nan'), 0.1))
test(e is not None, "refus : α non fini", e or "valeur renvoyée sans erreur")
ieB, resB, _ = eh.solve_fixed_point(0.07, 0.03); ieC, resC, _ = eh.solve_fixed_point(0.07, 0.02)
test(f"{ieB:.4f}" == "0.3104" and f"{ieC:.4f}" == "0.2695" and resB <= 1e-12 and resC <= 1e-12, "points fixes de B et C : 0,3104 et 0,2695 (tableau T3)", f"IED*_B = {ieB:.10f} ; IED*_C = {ieC:.10f}")

print(); print(f"BILAN : {BILAN['CONFORME']} tests conformes, {BILAN['NON CONFORME']} non conformes")
sys.exit(0 if BILAN['NON CONFORME'] == 0 else 1)
