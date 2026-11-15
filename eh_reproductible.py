#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Auteur : Benjamin Albouy
# Rattachement : programme de recherche « Économie homéostatique »
"""
================================================================================
Démonstrateur EH — MOTEUR DE RÉFÉRENCE ET PROTOCOLE REPRODUCTIBLE
Version V1.1 — novembre 2026 — concordance avec le Working Paper EH V1.0
(novembre 2026) et la note de cadrage P-F1 V1 (novembre 2026)
================================================================================
STATUT — dispositif de mise en calcul d'une instance particulière. Il ne calcule
ni ne caractérise le noyau de viabilité et ne résout pas P-F1. Vocabulaire des
documents de référence : C1 est une QUESTION de caractérisation (WP §5.5) ;
C2 et C4 (C4a–C4c) sont des CONJECTURES ; C3 est une propriété analytique de
l'IED. Les scénarios sont illustratifs PAR CONSTRUCTION (critères fonctionnels
posés a priori, cf. problème S2 du WP §9.5) : ils ne testent aucune conjecture ;
le scénario C joue le rôle de cas adverse (sortie de K).

Une seule exécution produit tous les tableaux numériques de l'article, écrit les
CSV correspondants dans ./csv/ et génère les figures de l'article dans
./figures/. Les nombres du manuscrit proviennent de ce fichier, à l'exception de
l'instant de franchissement continu du scénario C, produit par le script de
contrôle joint au paquet. Aucune figure ne repose sur une seconde implémentation
des équations.

DÉPENDANCES
 - calculs et CSV  : bibliothèque standard uniquement (math, csv, os, fractions).
 - génération des figures : dépendance OPTIONNELLE matplotlib.
   Si matplotlib est absent, les tableaux et les CSV sont produits normalement ;
   seule la génération graphique est ignorée, avec un message explicite. L'absence
   de matplotlib ne provoque jamais l'échec du protocole.
Compatible Python 3.8+.

USAGE D'OUTILS D'IA
Claude Fable 5.2 (Anthropic) a été utilisé pour la préparation de la version
V1.1 : contrôle de reproductibilité de la V1, application des corrections à
l'article et au moteur, scripts de contrôle et de test, constitution du paquet
de reproductibilité. L'auteur a vérifié les contenus produits et reste seul
responsable du contenu.

FIGURES — toutes les figures dynamiques sont construites à partir des sorties de
run() et de solve_fixed_point() (temps rows[k]["t"], IED rows[k]["IED"],
M rows[k]["M"], coefficient de régulation rows[k]["u"]) et des métriques d'événement déjà
calculées par le moteur ; aucun résultat n'est recopié en dur dans le code
graphique. Sortie en PNG 300 dpi et PDF vectoriel de même nom de base.

MOTEUR — conventions
 - horizon T, pas Δt, N = T/Δt intervalles et N+1 états ; l'état initial (t=0)
   est enregistré ; l'état obtenu après un pas d'Euler est associé à t+Δt ;
 - chaque enregistrement est temporellement cohérent : (t_k, S_k, M_k, IED_k, u_k)
   se rapportent au même instant t_k, u_k étant le coefficient de régulation appliqué en t_k ;
 - la masse monétaire M est conservée sans troncature, dans les données comme dans les
   indicateurs ;
 - min_{0<=t<=T} IED(t) est calculé sur tous les états ;
 - pour Z_cible, l'occupation initiale (état 0) est distinguée de la première
   entrée dynamique depuis l'extérieur.
Retard indexé en temps : u_k = phi_c(IED_{k-τsteps}), IED_{j<0}=IED_0,
avec phi_c(z)=min(z,1) sur [0,2].
Spécialisation monétaire : e_base(S) ≡ 1 (constante E_BASE du code), R(t) ≡ 0,
D(M;Λ)=λM ; ainsi dM/dt = u(t)·e_base − λM.
Système triangulaire : M n'entre pas dans dS ; IED* ne dépend que de (α,β,γ).
La classe affine retenue est restreinte aux paramètres rendant [0,2]^3 invariant ;
aucune projection post-Euler de S n'est utilisée.

POINT FIXE BIOPHYSIQUE BOUCLÉ : g(x)=G_fp(x)-x=0,
G_fp(x)=[Π_i ((β_i·phi_c(x)+γ_i)/α_i)]^(1/3).
Méthode : balayage de [1e-4, 2] en 4001 points à la recherche d'un changement
de signe de g, puis DEUX branches selon le résultat du balayage.
 - Changement de signe détecté : bissection sur l'encadrement, arrêtée lorsque
   la LARGEUR de l'encadrement devient inférieure à 1e-15 (au plus
   200 itérations) ; IED* est le milieu de l'encadrement final.
 - Aucun changement de signe sur la grille : si g s'annule exactement en un
   point de la grille, ce point est la racine ; sinon, si g(0) > 0 et
   g(1e-4) < 0, la racine est encadrée sur [0, 1e-4] et obtenue par la même
   bissection. Tout autre cas (point fixe non encadré sur [0, 2]) déclenche
   une erreur explicite.
Dans les deux cas le résidu |g(IED*)| est calculé a posteriori et rapporté ;
il n'est jamais utilisé comme critère d'arrêt. Un résidu supérieur à 1e-12
déclenche une erreur explicite : aucune valeur non convergée n'est renvoyée.
Les cinq scénarios retenus empruntent tous la première branche.
Simulation T=3000 = contrôle complémentaire uniquement.

DÉTECTEUR D'OSCILLATIONS
 - pente en IED PAR UNITÉ DE TEMPS : m_k = (IED_{k+1}-IED_k)/Δt ;
 - tolérance θ = 1e-4, exprimée en unités d'IED par unité de temps ;
   |m_k| < θ  -> pas plat (ignoré pour la direction) ; |m_k| >= θ -> significatif ;
 - retournement = changement de signe entre deux pentes significatives consécutives ;
 - 1 retournement = « extremum isolé ou dépassement transitoire » (NON compté
   comme oscillation) ; oscillation ssi nb de retournements >= R_min = 2 ;
 - fenêtre analysée = [0, T]. Répété à Δt=0.05 et 0.025.
 Portée : ne détecte que les retournements résolus par le pas ; ne prouve pas
 l'absence d'oscillation à toute échelle.

CIBLE MONÉTAIRE QUASI-STATIONNAIRE : M*(T) = [e(T)+R(T)]/λ = q(T)/λ, avec
q(T)=u(T)·e_base et R≡0 — c'est la cible quasi-stationnaire M*(t) du WP (§5.5)
évaluée au temps final ; ce n'est pas un équilibre du système complet. Flux
« stabilisé » ssi |q(T)-q(0.9T)|/|q(T)| < ε_flux = 1e-3 (test empirique de la
séparation d'échelles du WP, temps d'ajustement ~1/λ). Ce critère relatif est
défini pour q(T) > 0 (flux_regime = 'relatif'). Pour q(T) = q(0.9T) = 0, le flux
est rapporté comme nul constant (flux_regime = 'nul_constant', flux_rel = None,
flux_stab = False) ; pour q(T) = 0 et q(0.9T) > 0, comme éteint
(flux_regime = 'extinction', flux_rel = inf, flux_stab = False). Lorsque u≡1 sur la fin
de l'horizon, la stabilisation résulte de la SATURATION de phi_c et non d'un
équilibre biophysique : le moteur le signale (colonne origine_stabilisation).

CONTRÔLES DE GRILLE : T et τ sont des multiples entiers de Δt, les quotients
T/Δt et τ/Δt étant contrôlés avec une tolérance absolue de 1e-9 (distincte de
la tolérance 1e-12 des contrôles de domaine et d'admissibilité). L'horizon
comporte au moins un intervalle (N >= 1) et tout retard τ > 0 au moins un pas
(τ_steps >= 1).

BORNITUDE MONÉTAIRE ET PROPRIÉTÉ (P) DE LA NOTE P-F1 : avec u∈[0,1], e_base ≡ 1,
R≡0 et λ>0, on a 0 <= M(t) <= max(M0, 1/λ) pour tout t (et pour le schéma
d'Euler dès que λΔt<=1). Aucune explosion en temps fini n'est possible : (P) est
vérifiée pour cette instance. Le moteur contrôle cette borne le long de chaque
trajectoire (T2, T9).

ADMISSIBILITÉ AU SENS DE LA NOTE P-F1 (tableau T9) : (F1) continuité de F ;
(F3′) viabilité globale existentielle du cube, ici impliquée par l'invariance
forte (β_i+γ_i <= 2α_i, γ_i >= 0) ; (F4) non-trivialité (β>0) ; (L) F affine ;
phi_c ∈ 𝒫_Φ ; D(M;λ)=λM vérifie (D1)–(D4) ; e_base ≡ 1 et r ≡ 0 continues et
lipschitziennes. Pour τ>0 (scénario E) l'instance sort du périmètre EDO de P-F1
(extension exploratoire rattachée à C2/P-F2) : seule la partie τ=0 est qualifiée.

CERTIFICAT D'ÉQUILIBRE — INSTANCE B (tableau T10) : le scénario B du démonstrateur
coïncide avec l'instance B de la note P-F1 (§4). Le moteur vérifie en
arithmétique rationnelle exacte que le polynôme Π_i(β_i z+γ_i)/α_i − z³ est
proportionnel à Q(z) = −21620z³+1395z²+954z+216 donné par la note, localise sa
racine positive z*, vérifie 0.310 < z* < 0.311, z* < 1 (branche non saturée),
S* ∈ (0,2)^3, M* = z*/λ ∈ (31/15, 311/150) et (S*,M*) ∈ K strictement. Il ne
s'agit que de la condition suffisante élémentaire de P-F1 (Viab(K) ≠ ∅ pour
CETTE instance) ; ni le noyau ni l'ensemble des triplets admissibles ne sont
caractérisés.

RÉGULARITÉ : l'IED n'est pas lipschitzien en certains points du bord nul du cube
mais l'est sur tout domaine IED >= δ > 0. Le moteur rapporte min IED le long de
chaque trajectoire (T9) : c'est la valeur de δ pour laquelle le champ bouclé est
lipschitzien sur la région effectivement parcourue.

Z_CIBLE — définition retenue ici : Z_cible = {IED >= z_cible} ∩ {M ∈ [m_cible_lo, m_cible_hi]}.
RELATION AVEC LA ZONE DE SÉCURITÉ Z_sat DU WORKING PAPER. Le WP définit Z_sat
comme le resserrement des DEUX contraintes de K par des marges strictement
positives (éq. 6) : g_env(S) >= ε_env et g_social(M) >= ε_social, soit deux
bornes INFÉRIEURES, sans borne supérieure sur M. La zone cible Z_cible ajoute
une borne monétaire SUPÉRIEURE m_cible_hi, choix propre à l'instance calculée
(cible bornée dans le plan de projection (IED, M)). Conséquences : (i)
l'inclusion Z_cible ⊂ K est vérifiée dans les scénarios retenus, car
z_cible > k_env et m_cible_lo > k_social ; (ii) M(70) dépasse m_cible_hi dans
les scénarios A, B, C, E (tableau T2) ; pour A, où IED(70) >= z_cible, cette
borne est la seule cause d'exclusion de Z_cible et n'a pas d'équivalent dans la
zone de sécurité du WP ; pour B, C et E, IED(70) < z_cible exclut aussi l'état.
Z_CIBLE — métriques par trajectoire : nombre d'entrées dynamiques ; premier
instant d'entrée ; dernier instant de sortie ; durée du plus long séjour continu,
définie comme (L-1)·Δt pour L états consécutifs dans Z_cible.
K — événements par trajectoire : init_in_k indique si l'état initial appartient
à K ; k_first_entry_t est le premier instant d'entrée DYNAMIQUE dans K (passage
de l'extérieur vers K, k >= 1 ; l'occupation initiale n'est pas une entrée) ;
exit_t est le premier instant d'une transition de K vers son extérieur
(ink[k-1] vrai, ink[k] faux), quel que soit l'état initial, et cause décrit le
franchissement constaté à cet instant. Une trajectoire issue de l'extérieur de
K qui n'y entre jamais n'a donc pas de sortie : exit_t reste None et
init_in_k=False permet de ne pas l'interpréter comme un maintien dans K.

TRANSITION β : β_entry^Δt(T;p_{−β}) — transition NUMÉRIQUE dépendante du
protocole (T, x0, λ, Z_cible, Δt), non un seuil intrinsèque. Encadrement généré
automatiquement à partir de last_no / first_yes.

ROBUSTESSE NUMÉRIQUE : « comportement qualitativement inchangé lors du
raffinement du pas de 0.05 à 0.025 » ssi concordent : présence ET cause de
sortie ; signe de la marge au seuil environnemental ; présence d'une entrée dynamique ; état
d'oscillation ; instants d'événement (dans une tolérance liée au pas).
================================================================================
"""
import math, csv, os

# Dépendance graphique optionnelle : son absence ne doit pas empêcher les
# tableaux ni les CSV d'être produits.
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle, FancyArrowPatch
    HAS_MPL = True
    MPL_ERR = None
except Exception as _e:            # ImportError ou toute erreur de backend
    HAS_MPL = False
    MPL_ERR = repr(_e)

GAMMA = [0.02, 0.018, 0.012]
E_BASE = 1.0                    # e_base(S) ≡ 1 dans cette instance ; R(t) ≡ 0
STATE_TOL = 1e-12               # tolérance de contrôle des bornes numériques
EPS_FLUX = 1e-3
FP_RES_TOL = 1e-12        # résidu maximal admis pour le point fixe bouclé
OSC_THETA_RATE = 1e-4          # tolérance en IED par unité de temps
OSC_RMIN = 2
EVENT_TOL = 3 * 0.05          # tolérance d'instant liée au pas grossier (quelques pas)
CSV_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "csv")
FIG_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figures")

def cbrt(x): return math.copysign(abs(x)**(1/3), x) if x else 0.0
def sgn(x): return (x > 0) - (x < 0)

def control_phi(z):
    """Rétroaction prescrite du démonstrateur : phi_c(z)=min(z,1) sur [0,2].

    Statut au regard des conditions posées sur Φ par le WP V1.0 (§5.5) et par la
    note P-F1 V1 (§3) : continue sur [0,2], croissante au sens large (condition
    retenue par ces deux documents), phi_c(0)=0, phi_c(1)=1, à valeurs dans
    U=[0,u_max] avec u_max=1 (valeur minimale admise, u_max >= 1). phi_c appartient
    donc à la classe 𝒫_Φ ; sa constance sur [1,2] (saturation) est admissible.

    Conséquence de u_max=1 : le régime d'expansion u>1 décrit par le WP est
    structurellement exclu de cette instance ; seuls les régimes de contraction
    (u<1) et neutre (u=1) sont représentés. Le candidat exploratoire du WP
    Φ₀(z)=z, qui exige u_max >= 2, n'est pas utilisé ici. La part de l'horizon à
    coefficient de régulation non saturé est mesurée par l'indicateur u_unsat renvoyé par run().
    """
    if z < -STATE_TOL or z > 2.0 + STATE_TOL:
        raise ValueError(f"IED hors du domaine [0,2] : {z}")
    return min(max(z, 0.0), 1.0)

def validate_parameters(p, dt):
    """Vérifie que l'instance affine respecte le domaine d'état annoncé.

    Pour dS_i/dt = -alpha_i S_i + beta_i u + gamma_i, u∈[0,1],
    l'invariance de [0,2]^3 est garantie si le champ est rentrant aux bornes.
    La borne inférieure est automatique pour beta_i,gamma_i>=0 ; la borne
    supérieure impose beta_i + gamma_i <= 2 alpha_i pour chaque composante.
    """
    a, b, lam = p['alpha'], p['beta'], p['lambda']
    for name in ('alpha', 'beta', 'lambda', 'tau', 'ibd0', 'iee0', 'irnr0', 'm0',
                 'k_env', 'k_social', 'z_cible', 'm_cible_lo', 'm_cible_hi'):
        if not math.isfinite(p[name]):
            raise ValueError(f"Paramètre non fini : {name}={p[name]}")
    if not math.isfinite(dt):
        raise ValueError(f"Pas de temps non fini : {dt}")
    if a <= 0 or b < 0 or lam <= 0 or p['tau'] < 0:
        raise ValueError(f"Paramètres invalides : alpha={a}, beta={b}, lambda={lam}, tau={p['tau']}")
    if dt <= 0:
        raise ValueError(f"Pas de temps non positif : {dt}")
    AL = [a, a, a*1.3]; BE = [b, b, b*0.5]
    for j in range(3):
        if AL[j]*dt > 1.0 + STATE_TOL:
            raise ValueError(
                f"Pas trop grand pour préserver [0,2] par Euler : alpha_i*dt={AL[j]*dt} "
                f"(composante {j+1})."
            )
        outward_upper = -2.0*AL[j] + BE[j] + GAMMA[j]
        if outward_upper > STATE_TOL:
            raise ValueError(
                f"Paramètres incompatibles avec l'invariance de [0,2]^3 : "
                f"composante {j+1}, F_i(2,u=1)={outward_upper:.6g} > 0 "
                f"(alpha={a}, beta={b})."
            )
    initials = (p['ibd0'], p['iee0'], p['irnr0'])
    if any(v < -STATE_TOL or v > 2.0 + STATE_TOL for v in initials):
        raise ValueError(f"Condition initiale S0 hors [0,2]^3 : {initials}")
    if p['m0'] < -STATE_TOL:
        raise ValueError(f"Condition initiale M0 négative : {p['m0']}")
    # Pour l'Euler explicite monétaire avec entrée non négative, lambda*dt<=1
    # évite une négativité numérique créée par le seul pas de discrétisation.
    if lam*dt > 1.0 + STATE_TOL:
        raise ValueError(f"Pas trop grand pour préserver M>=0 par Euler : lambda*dt={lam*dt}")

def detect_oscillations(traj, dt, theta_rate=OSC_THETA_RATE, rmin=OSC_RMIN):
    """Pentes m_k=(IED_{k+1}-IED_k)/Δt (IED/temps). Retournements = changements de
    signe entre pentes significatives (|m|>=θ). 1 retournement = extremum isolé ou
    dépassement transitoire ; oscillation ssi retournements>=R_min."""
    slope = [(traj[k+1]-traj[k])/dt for k in range(len(traj)-1)]
    signs = [(1 if m >= theta_rate else (-1 if m <= -theta_rate else 0)) for m in slope]
    nz = [s for s in signs if s != 0]
    reversals = sum(1 for k in range(len(nz)-1) if nz[k] != nz[k+1])
    return reversals, (reversals >= rmin)

def zcible_metrics(inz, times, dt):
    """Nombre d'entrées dynamiques (False->True, k>=1), 1er instant d'entrée,
    dernier instant de sortie (True->False), durée du plus long séjour continu."""
    entries = [k for k in range(1, len(inz)) if inz[k] and not inz[k-1]]
    exits   = [k for k in range(1, len(inz)) if not inz[k] and inz[k-1]]
    longest = cur = 0
    for v in inz:
        cur = cur + 1 if v else 0
        longest = max(longest, cur)
    return dict(n_entries=len(entries),
                first_entry_t=(times[entries[0]] if entries else None),
                last_exit_t=(times[exits[-1]] if exits else None),
                longest_stay=(round((longest-1)*dt, 4) if longest > 0 else 0.0))

def flux_terminal(q_T, q_90):
    """Régime du flux terminal, à partir de q(T) et de q(0.9T), tous deux >= 0.
    Renvoie (flux_regime, flux_rel, flux_stab) :
      - q(T) > 0            : 'relatif', flux_rel = |q(T)-q(0.9T)|/q(T), flux_stab = (flux_rel < EPS_FLUX) ;
      - q(T) = q(0.9T) = 0  : 'nul_constant', flux_rel = None (critère relatif non défini), flux_stab = False ;
      - q(T) = 0 < q(0.9T)  : 'extinction', flux_rel = inf, flux_stab = False."""
    if q_T > 0:
        rel = abs(q_T - q_90)/q_T
        return 'relatif', rel, (rel < EPS_FLUX)
    if q_90 == 0:
        return 'nul_constant', None, False
    return 'extinction', float('inf'), False

# ─────────────────────────── MOTEUR ───────────────────────────
def run(p, dt=0.05, T=70.0):
    validate_parameters(p, dt)
    if not math.isfinite(T) or T <= 0:
        raise ValueError(f"Horizon non fini ou non positif : {T}")
    for name, val in (('T', T), ('tau', p['tau'])):
        ratio = val/dt
        if abs(ratio - round(ratio)) > 1e-9:
            raise ValueError(f"{name}/dt = {ratio} n'est pas entier : choisir un pas diviseur de {name} "
                             "(aucun arrondi silencieux de l'horizon ni du retard).")
    N = round(T/dt)
    if N < 1:
        raise ValueError(f"Horizon T={T} sans intervalle au pas dt={dt} (N = 0) : choisir T >= dt.")
    a, b, lam = p['alpha'], p['beta'], p['lambda']
    tau_steps = round(p['tau']/dt)
    if p['tau'] > 0 and tau_steps < 1:
        raise ValueError(f"Retard tau={p['tau']} sans pas entier au pas dt={dt} (tau_steps = 0) : "
                         "choisir tau = 0 ou tau >= dt.")
    AL = [a, a, a*1.3]; BE = [b, b, b*0.5]
    kenv, ksoc = p['k_env'], p['k_social']
    zs, mlo, mhi = p['z_cible'], p['m_cible_lo'], p['m_cible_hi']
    S = [[p['ibd0'], p['iee0'], p['irnr0']]]; M = [p['m0']]
    ied0 = cbrt(S[0][0]*S[0][1]*S[0][2]); hist = []
    rows = []
    for k in range(N+1):
        ied_k = cbrt(S[k][0]*S[k][1]*S[k][2]); hist.append(ied_k)
        idx = k - tau_steps
        u_k = control_phi(hist[idx] if idx >= 0 else ied0)
        in_k = (ied_k >= kenv) and (M[k] >= ksoc)
        # Z_cible du démonstrateur : seuil sur IED + FENÊTRE sur M (borne
        # supérieure mhi sans contrepartie dans l'éq. (6) du WP, cf. en-tête).
        in_z = (ied_k >= zs) and (mlo <= M[k] <= mhi)
        rows.append(dict(t=round(k*dt, 6), IED=ied_k, M=M[k], u=u_k, in_k=in_k, in_z=in_z))
        if k < N:
            dS = [-AL[j]*S[k][j] + BE[j]*u_k + GAMMA[j] for j in range(3)]
            S_next = [S[k][j] + dS[j]*dt for j in range(3)]
            if any(v < -STATE_TOL or v > 2.0 + STATE_TOL for v in S_next):
                raise RuntimeError(
                    f"Le pas d'Euler quitte [0,2]^3 à t={k*dt:.6g} : S_next={S_next}. "
                    "Réduire Δt ou revoir les paramètres ; aucune projection silencieuse n'est appliquée."
                )
            S.append(S_next)
            M_next = M[k] + (u_k*E_BASE - lam*M[k])*dt
            if M_next < -STATE_TOL:
                raise RuntimeError(
                    f"Le pas d'Euler produit M<0 à t={k*dt:.6g} : M_next={M_next}. "
                    "Réduire Δt ; aucune projection silencieuse n'est appliquée."
                )
            M.append(M_next)
    ied = [r['IED'] for r in rows]; ink = [r['in_k'] for r in rows]
    inz = [r['in_z'] for r in rows]; times = [r['t'] for r in rows]
    zm = zcible_metrics(inz, times, dt)
    dyn_entry = zm['first_entry_t']
    # Sortie de K : première transition K -> extérieur, quel que soit l'état
    # initial (pour un état initial dans K, c'est le premier état hors de K).
    exit_t = cause = None
    for k in range(1, len(ink)):
        if ink[k-1] and not ink[k]:
            exit_t = rows[k]['t']
            cause = ('env' if rows[k]['IED'] < kenv else '') + ('mon' if rows[k]['M'] < ksoc else '')
            break
    # Entrée dynamique dans K : première transition extérieur -> K (k >= 1).
    k_entries = [k for k in range(1, len(ink)) if ink[k] and not ink[k-1]]
    k_first_entry_t = times[k_entries[0]] if k_entries else None
    q_T = rows[-1]['u'] * E_BASE
    q_90 = rows[int(0.9*N)]['u'] * E_BASE
    flux_regime, flux_rel, flux_stab = flux_terminal(q_T, q_90)
    nrev, osc = detect_oscillations(ied, dt)
    return dict(rows=rows, ied_traj=ied, ied0=ied0, ied_final=ied[-1], ied_min=min(ied),
                M0=p['m0'], M_final=rows[-1]['M'], M_max=max(r['M'] for r in rows),
                q_terminal=q_T, M_star_T=(q_T/lam if lam > 0 else float('inf')),
                M_bound=max(p['m0'], 1.0*E_BASE/lam), M_bound_ok=(max(r['M'] for r in rows) <= max(p['m0'], 1.0*E_BASE/lam) + 1e-12),
                u_sat_final=all(r['u'] >= 1.0 for r in rows[int(0.9*N):]),
                flux_rel=flux_rel, flux_stab=flux_stab, flux_regime=flux_regime,
                pct_K=round(sum(ink)/len(ink)*100, 1), init_in_z=inz[0], dyn_entry=dyn_entry,
                dwell_z=round(sum(inz)/len(inz)*100, 1), exit_t=exit_t, cause=cause,
                init_in_k=ink[0], k_first_entry_t=k_first_entry_t,
                last_in_k=ink[-1], u_unsat=round(sum(1 for r in rows if r['u'] < 1.0)/len(rows)*100, 1),
                n_reversals=nrev, oscillation=osc, margin_env_threshold=min(ied)-kenv, **{'z_'+k: v for k, v in zm.items()})

def fixed_point_map(x, a, b):
    """Carte scalaire induite par les équilibres biophysiques à IED supposé x.

    Ce n'est pas la fonction de couplage théorique Φ ; le coefficient de régulation est control_phi(x).
    """
    AL = [a, a, a*1.3]; BE = [b, b, b*0.5]
    u = control_phi(x)
    s = [(BE[j]*u + GAMMA[j])/AL[j] for j in range(3)]
    if any(v < -STATE_TOL or v > 2.0 + STATE_TOL for v in s):
        raise ValueError(f"Carte de point fixe hors [0,2]^3 : {s} (alpha={a}, beta={b})")
    return cbrt(s[0]*s[1]*s[2]), s

def solve_fixed_point(a, b):
    g = lambda x: fixed_point_map(x, a, b)[0] - x
    xs = [1e-4 + i*(2.0-1e-4)/4000 for i in range(4001)]
    lo = hi = None; gp = g(xs[0])
    for x in xs[1:]:
        gc = g(x)
        if gp*gc < 0: lo, hi = x-(2.0-1e-4)/4000, x; break
        gp = gc
    if lo is None:
        # Aucun changement de signe strict sur la grille : racine en un point de la
        # grille, ou racine sous la borne basse du balayage (g(0) > 0 > g(1e-4)).
        zeros = [x for x in xs if g(x) == 0.0]
        if zeros:
            lo = hi = zeros[0]
        elif g(0.0) > 0 and g(xs[0]) < 0:
            lo, hi = 0.0, xs[0]
        else:
            raise ValueError(f"Point fixe non encadré sur [0, 2] pour alpha={a}, beta={b}.")
    for _ in range(200):
        mid = 0.5*(lo+hi)
        if g(lo)*g(mid) <= 0: hi = mid
        else: lo = mid
        if hi-lo < 1e-15: break
    ied = 0.5*(lo+hi)
    _, Ss = fixed_point_map(ied, a, b)
    res = abs(g(ied))
    if not res <= FP_RES_TOL:
        raise RuntimeError(f"Point fixe non convergé pour alpha={a}, beta={b} : résidu {res:.3e} > {FP_RES_TOL:g}.")
    return ied, res, Ss

# ─────────────────────── SCÉNARIOS RETENUS ───────────────────────
# Retenus après vérification fonctionnelle et test minimal de robustesse numérique.
def P(**k): k['lambda'] = k.pop('lam'); return k
SCN = {
 "A_reference":  P(alpha=0.10, beta=0.18, lam=0.12, tau=0, ibd0=1.5, iee0=1.4, irnr0=1.2, m0=1.2,
                   k_env=0.40, k_social=0.30, z_cible=0.9, m_cible_lo=0.8, m_cible_hi=1.5),
 "B_approche_seuil":    P(alpha=0.07, beta=0.03, lam=0.15, tau=0, ibd0=0.7, iee0=0.7, irnr0=0.5, m0=0.8,
                   k_env=0.30, k_social=0.20, z_cible=0.9, m_cible_lo=0.7, m_cible_hi=1.4),
 "C_sortie_K":   P(alpha=0.07, beta=0.02, lam=0.15, tau=0, ibd0=0.7, iee0=0.7, irnr0=0.5, m0=0.8,
                   k_env=0.30, k_social=0.20, z_cible=0.9, m_cible_lo=0.7, m_cible_hi=1.4),
 "E_delai":      P(alpha=0.22, beta=0.28, lam=0.15, tau=10, ibd0=0.45, iee0=0.45, irnr0=0.4, m0=0.8,
                   k_env=0.20, k_social=0.15, z_cible=1.2, m_cible_lo=0.7, m_cible_hi=1.4),
 "F_acces_Zcible": P(alpha=0.12, beta=0.16, lam=0.82, tau=0, ibd0=0.6, iee0=0.55, irnr0=0.5, m0=1.0,
                   k_env=0.30, k_social=0.25, z_cible=0.9, m_cible_lo=0.7, m_cible_hi=1.5),
}
BETA_GRID = [i/100 for i in range(10, 23)]  # 0.10..0.22, classe affine invariant [0,2]^3 pour alpha=0.12

LABELS = {"A_reference":"A · Référence — maintien sur horizon fini",
          "B_approche_seuil":"B · Approche progressive du seuil environnemental",
          "C_sortie_K":"C · Sortie environnementale de K",
          "E_delai":"E · Délai actif sans oscillation observée",
          "F_acces_Zcible":"F · Accès conditionnel à Z_cible"}

def emit(fname, header, data):
    os.makedirs(CSV_DIR, exist_ok=True)
    with open(os.path.join(CSV_DIR, fname), "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f); w.writerow(header); w.writerows(data)
    print("  [CSV] csv/" + fname)

def qualitative_concordance(sc):
    r5, r25 = run(sc, 0.05), run(sc, 0.025)
    def close(a, b):
        if a is None and b is None: return True
        if a is None or b is None: return False
        return abs(a-b) <= EVENT_TOL
    chk = {
      'sortie_présence': (r5['exit_t'] is None) == (r25['exit_t'] is None),
      'sortie_cause':    r5['cause'] == r25['cause'],
      'signe_marge_seuil_env': sgn(r5['margin_env_threshold']) == sgn(r25['margin_env_threshold']),
      'entrée_présence': (r5['dyn_entry'] is None) == (r25['dyn_entry'] is None),
      'état_oscillation': r5['oscillation'] == r25['oscillation'],
      'instant_sortie':  close(r5['exit_t'], r25['exit_t']),
      'instant_entrée':  close(r5['dyn_entry'], r25['dyn_entry']),
    }
    return all(chk.values()), chk

def bar(): print("-"*100)

# ═════════════════════════════ TABLEAUX ═════════════════════════════
def t_robustesse():
    print("\n[T1] ROBUSTESSE — Δt=0.05 et 0.025, T=70 (test qualitatif renforcé)"); bar()
    hdr = ["scenario","dt","IED0","IED_min(0..70)","IED_final","pct_K","entree_Zcible","sortie_K","cause","oscillation"]
    rows = []
    print(f"{'scénario':38}{'dt':>6}{'IED0':>6}{'min':>7}{'fin':>7}{'%K':>6}{'entréeZ':>8}{'sortieK':>8}{'osc':>5}")
    for key, sc in SCN.items():
        for dt in (0.05, 0.025):
            r = run(sc, dt)
            rows.append([key, dt, f"{r['ied0']:.4f}", f"{r['ied_min']:.4f}", f"{r['ied_final']:.4f}",
                         r['pct_K'], r['dyn_entry'], r['exit_t'], r['cause'], r['oscillation']])
            print(f"{LABELS[key][:38]:38}{dt:>6}{r['ied0']:6.3f}{r['ied_min']:7.3f}{r['ied_final']:7.3f}"
                  f"{r['pct_K']:6.1f}{str(r['dyn_entry']):>8}{str(r['exit_t']):>8}{str(r['oscillation']):>5}")
    emit("T1_robustesse.csv", hdr, rows)
    print("  Verdict (test renforcé : sortie+cause, signe marge, entrée, oscillation, instants) :")
    for key, sc in SCN.items():
        ok, chk = qualitative_concordance(sc)
        detail = "" if ok else "  [écarts: " + ", ".join(k for k, v in chk.items() if not v) + "]"
        print(f"    {LABELS[key][:38]:38} comportement qualitativement "
              f"{'inchangé' if ok else 'CHANGÉ'} lors du raffinement du pas de 0.05 à 0.025{detail}")

def t_monetaire():
    print("\n[T2] MONÉTAIRE — cible quasi-stationnaire M*(T)=[e(T)+R(T)]/λ (WP §5.5) ; borne analytique 0<=M<=max(M0,1/λ)"); bar()
    hdr = ["scenario","lambda","q_T","M(70)","M_star_T","flux_rel(90-100%)","flux_stabilise(<1e-3)",
           "origine_stabilisation","u_unsat_pct","borne_max(M0,1/lambda)","borne_respectee"]
    rows = []
    print(f"{'scénario':38}{'λ':>6}{'q_T':>8}{'M(70)':>8}{'M*(T)':>8}{'flux_rel':>10} {'flux':>14} {'origine':>12}{'u_unsat%':>9}{'borne':>7}")
    for key, sc in SCN.items():
        r = run(sc)
        orig = ("saturation u=1" if (r['flux_stab'] and r['u_sat_final']) else
                ("autre" if r['flux_stab'] else "-"))
        rows.append([key, sc['lambda'], f"{r['q_terminal']:.4f}", f"{r['M_final']:.3f}",
                     f"{r['M_star_T']:.3f}", f"{r['flux_rel']:.3e}", r['flux_stab'], orig,
                     r['u_unsat'], f"{r['M_bound']:.3f}", r['M_bound_ok']])
        print(f"{LABELS[key][:38]:38}{sc['lambda']:6.2f}{r['q_terminal']:8.3f}{r['M_final']:8.2f}"
              f"{r['M_star_T']:8.2f}{r['flux_rel']:10.2e} {('stabilisé' if r['flux_stab'] else 'NON stabilisé'):>14} {orig:>12}"
              f"{r['u_unsat']:9.1f}{('ok' if r['M_bound_ok'] else 'VIOL'):>7}")
    emit("T2_monetaire.csv", hdr, rows)

def t_pointsfixes():
    print("\n[T3] POINTS FIXES — g(x)=G_fp(x)−x=0 (bissection) ; contrôle simulation T=3000"); bar()
    hdr = ["scenario","alpha","beta","IED_star","residu","S1","S2","S3","IED_dyn_T3000","k_env","marge"]
    rows = []
    for key in ("B_approche_seuil", "C_sortie_K"):
        sc = SCN[key]; ie, res, Ss = solve_fixed_point(sc['alpha'], sc['beta'])
        iedL = run(sc, dt=0.05, T=3000.0)['ied_traj'][-1]; marg = ie - sc['k_env']
        rows.append([key, sc['alpha'], sc['beta'], f"{ie:.6f}", f"{res:.2e}",
                     f"{Ss[0]:.4f}", f"{Ss[1]:.4f}", f"{Ss[2]:.4f}", f"{iedL:.6f}", sc['k_env'], f"{marg:+.4f}"])
        print(f"  {LABELS[key]} : IED*={ie:.6f} (résidu {res:.2e}) ; contrôle T=3000={iedL:.6f} ; "
              f"IED*-k_env={marg:+.4f} (strict>: {ie > sc['k_env']})")
    rB = run(SCN["B_approche_seuil"])
    print(f"  B : min_(0<=t<=70) IED(t) = {rB['ied_min']:.4f} -> marge minimale sur T=70 = {rB['margin_env_threshold']:.4f}")
    print("  Note : IED* = point fixe numérique ; aucune stabilité asymptotique démontrée ; maintien formulé sur T=70.")
    emit("T3_points_fixes.csv", hdr, rows)

def t_c31_statique():
    print("\n[T4] C3.1 — ILLUSTRATION STATIQUE (algèbre pure, hors moteur) : IED(r)=(1.5·1.5·r)^(1/3)"); bar()
    hdr = ["IRNR","IED_statique"]; rows = []
    print(f"{'IRNR':>8}{'IED(r)':>14}")
    for r in IRNR_STAT:
        v = ied_static(r); rows.append([r, f"{v:.4f}"]); print(f"{r:8.2f}{v:14.4f}")
    emit("T4_C31_statique.csv", hdr, rows)

IRNR_STAT = [1.5, 1.0, 0.5, 0.2, 0.1, 0.05, 0.01]
IRNR_DYN = [1.5, 1.0, 0.5, 0.2, 0.05]
def d_family(irnr):
    """Famille dynamique D (§6 de l'article) : seules les CI varient."""
    return P(alpha=0.08, beta=0.12, lam=0.15, tau=0, ibd0=1.5, iee0=1.5, irnr0=irnr, m0=1.0,
             k_env=0.3, k_social=0.25, z_cible=0.9, m_cible_lo=0.7, m_cible_hi=1.4)
def ied_static(r): return (1.5*1.5*r)**(1/3)
def delay_gaps(dt=0.05, taus=(3, 6, 10, 15)):
    base = SCN["E_delai"]; b0 = dict(base); b0['tau'] = 0; tj0 = run(b0, dt)['ied_traj']; out = {}
    for tau in taus:
        bt = dict(base); bt['tau'] = tau; rt = run(bt, dt); tj = rt['ied_traj']
        out[tau] = (max(abs(tj[k]-tj0[k]) for k in range(min(len(tj), len(tj0)))), rt)
    return out

def t_c31_dynamique():
    print("\n[T5] C3.1 — FAMILLE DYNAMIQUE DISTINCTE (3 composantes évoluent ; seules les CI fixées)"); bar()
    hdr = ["IRNR0","IED0","IED_min","IED_final","observation"]; rows = []
    print(f"{'IRNR0':>7}{'IED0':>7}{'IEDmin':>8}{'IEDfin':>8}")
    for irnr in IRNR_DYN:
        r = run(d_family(irnr)); obs = "IED0 abaissé puis remonte (pas de chute dynamique)" if r['ied_final'] > r['ied0']+1e-3 else "-"
        rows.append([irnr, f"{r['ied0']:.4f}", f"{r['ied_min']:.4f}", f"{r['ied_final']:.4f}", obs])
        print(f"{irnr:7.2f}{r['ied0']:7.3f}{r['ied_min']:8.3f}{r['ied_final']:8.3f}   {obs}")
    emit("T5_C31_dynamique.csv", hdr, rows)

def t_delai():
    print("\n[T6] SCÉNARIO E — délai : écart de trajectoire + détecteur d'oscillations (pente/Δt)"); bar()
    print(f"     pentes m_k=(IED_k+1−IED_k)/Δt ; θ={OSC_THETA_RATE} IED/temps ; R_min={OSC_RMIN}")
    print("     (1 retournement = extremum isolé ou dépassement transitoire ; ≥2 = oscillation)")
    hdr = ["tau","dt","max_abs_dIED_vs_tau0","IED_final","u_unsat_pct","nb_retournements","oscillation"]; rows = []
    print(f"{'τ':>4}{'dt':>7}{'max|ΔIED|':>12}{'IEDfin':>8}{'u_unsat%':>9}{'#retour.':>9}{'oscill.':>8}")
    base = SCN["E_delai"]
    for dt in (0.05, 0.025):
        b0 = dict(base); b0['tau'] = 0; tj0 = run(b0, dt)['ied_traj']
        for tau in [0, 3, 6, 10, 15]:
            bt = dict(base); bt['tau'] = tau; rt = run(bt, dt); tj = rt['ied_traj']
            gap = max(abs(tj[k]-tj0[k]) for k in range(min(len(tj), len(tj0)))) if tau else 0.0
            rows.append([tau, dt, f"{gap:.4f}", f"{tj[-1]:.4f}", rt['u_unsat'], rt['n_reversals'], rt['oscillation']])
            print(f"{tau:>4}{dt:>7}{gap:12.4f}{tj[-1]:8.3f}{rt['u_unsat']:9.1f}{rt['n_reversals']:9}{str(rt['oscillation']):>8}")
    emit("T6_delai_tau.csv", hdr, rows)
    print('  Conclusion : « Aucun retournement significatif et répété n\'a été détecté selon le critère')
    print('  défini, pour les paramètres, les délais et l\'horizon testés. »')

def t_beta():
    print("\n[T7] SCÉNARIO F — β_entry^Δt(70; p_{−β}) : grille fine, robustesse, métriques Z_cible"); bar()
    sc = SCN["F_acces_Zcible"]; ied0 = cbrt(sc['ibd0']*sc['iee0']*sc['irnr0'])
    print(f"     IED0 = {ied0:.3f} < z_cible = {sc['z_cible']} (départ strictement hors Z_cible)")
    hdr = ["beta","eq_IED","entree_dt.05","entree_dt.025","nb_entrees","1er_entree","dernier_sortie",
           "plus_long_sejour","sejour_pct","M_final","M_dans_fenetre","pct_K"]
    rows = []
    print(f"{'β':>5}{'eqIED':>7}{'ent.05':>8}{'ent.025':>9}{'nE':>3}{'1ère':>7}{'dSort':>7}{'+long':>7}"
          f"{'séj%':>7}{'Mfin':>6}{'fen':>5}{'%K':>6}")
    last_no = first_yes = None; entries = []; dwell = {}
    for b in BETA_GRID:
        p = dict(sc); p['beta'] = b
        r = run(p, 0.05); r2 = run(p, 0.025); e, _, _ = solve_fixed_point(sc['alpha'], b)
        inwin = (sc['m_cible_lo'] <= r['M_final'] <= sc['m_cible_hi']); dwell[b] = r['dwell_z']
        if r['dyn_entry'] is None: last_no = b
        elif first_yes is None: first_yes = b
        if r['dyn_entry'] is not None: entries.append((b, r['dyn_entry']))
        rows.append([f"{b:.2f}", f"{e:.4f}", r['dyn_entry'], r2['dyn_entry'], r['z_n_entries'],
                     r['z_first_entry_t'], r['z_last_exit_t'], r['z_longest_stay'], r['dwell_z'],
                     f"{r['M_final']:.3f}", inwin, r['pct_K']])
        print(f"{b:5.2f}{e:7.3f}{str(r['dyn_entry']):>8}{str(r2['dyn_entry']):>9}{r['z_n_entries']:>3}"
              f"{str(r['z_first_entry_t']):>7}{str(r['z_last_exit_t']):>7}{r['z_longest_stay']:>7.2f}"
              f"{r['dwell_z']:>7.1f}{r['M_final']:>6.2f}{('oui' if inwin else 'NON'):>5}{r['pct_K']:>6.1f}")
    emit("T7_beta_grid.csv", hdr, rows)
    dec = all(entries[i+1][1] <= entries[i][1]+1e-9 for i in range(len(entries)-1))
    peak = max(dwell, key=lambda k: dwell[k])
    dwell_vals = [dwell[b] for b in BETA_GRID]
    nonmono = any(dwell_vals[i+1] < dwell_vals[i] - 1e-9 for i in range(len(dwell_vals)-1))
    print(f"\n  Encadrement généré automatiquement : {last_no} < β_entry^Δt(70; p_(−β)) ≤ {first_yes}")
    print(f"  (transition numérique dépendante du protocole, non un seuil intrinsèque)")
    print(f"  Instant de 1re entrée décroissant sur la grille testée : {dec}")
    print(f"  Durée d'occupation non monotone sur la grille : {nonmono} ; maximum observé vers β={peak:.2f}.")
    # métriques Z_cible détaillées pour le scénario retenu (β=0.16)
    rF = run(sc)
    print(f"  Z_cible (scénario retenu β={sc['beta']}) : nb_entrées={rF['z_n_entries']} ; "
          f"1er instant d'entrée={rF['z_first_entry_t']} ; dernier instant de sortie={rF['z_last_exit_t']} ; "
          f"plus long séjour continu={rF['z_longest_stay']}")

def t_registres():
    print("\n[T8] TABLEAU À QUATRE REGISTRES"); bar()
    rA, rB, rC = run(SCN["A_reference"]), run(SCN["B_approche_seuil"]), run(SCN["C_sortie_K"])
    ieB, resB, _ = solve_fixed_point(SCN["B_approche_seuil"]['alpha'], SCN["B_approche_seuil"]['beta'])
    ieC, _, _ = solve_fixed_point(SCN["C_sortie_K"]['alpha'], SCN["C_sortie_K"]['beta'])
    rF = run(SCN["F_acces_Zcible"])
    beta_runs = []
    for b in BETA_GRID:
        p = dict(SCN["F_acces_Zcible"]); p['beta'] = b
        beta_runs.append((b, run(p)))
    no_entry = [b for b, r in beta_runs if r['dyn_entry'] is None]
    yes_entry = [b for b, r in beta_runs if r['dyn_entry'] is not None]
    last_no = max(no_entry) if no_entry else None
    first_yes = min(yes_entry) if yes_entry else None
    dwell_seq = [r['dwell_z'] for _, r in beta_runs]
    dwell_nonmonotone = any(dwell_seq[i+1] < dwell_seq[i] - 1e-9 for i in range(len(dwell_seq)-1))
    scen_runs = [run(sc) for sc in SCN.values()]
    m_min, m_max = min(r['M_final'] for r in scen_runs), max(r['M_final'] for r in scen_runs)
    stabilized = [key for key, sc in SCN.items() if run(sc)['flux_stab']]
    # position de M*(T) relativement à la fenêtre monétaire de Z_cible, calculée
    pos_win = {"en deçà": [], "à l'intérieur": [], "au-delà": []}
    for key, sc in SCN.items():
        ms = run(sc)['M_star_T']
        pos_win["en deçà" if ms < sc['m_cible_lo'] else ("au-delà" if ms > sc['m_cible_hi'] else "à l'intérieur")].append(key[0])
    pos_txt = " ; ".join(f"{k} : {','.join(v)}" for k, v in pos_win.items() if v)
    kB = SCN["B_approche_seuil"]['k_env']
    stat_r = [1.5, 0.5, 0.1, 0.01]
    dyn = [run(d_family(i)) for i in IRNR_DYN]
    gaps = delay_gaps()
    gmax = max(g for g, _ in gaps.values())
    nrev_all = [r['n_reversals'] for _, r in gaps.values()]
    uuns = min(r['u_unsat'] for _, r in gaps.values())
    hdr = ["theme", "observation_directe", "calcul_derive", "interpretation", "limite_de_portee"]
    data = [
     ["A", f"pct_K={rA['pct_K']:g} ; min_(0<=t<=70)IED={rA['ied_min']:.3f} ; M(70)={rA['M_final']:.2f}",
      f"M*(T)={rA['M_star_T']:.2f} ; flux_rel={rA['flux_rel']:.1e} ; u_unsat={rA['u_unsat']}% (stabilisation par saturation u=1)",
      "maintien de la trajectoire dans K sur l'horizon, sans proximité du bord",
      "interne à l'instance affine, ces paramètres, horizon T=70"],
     ["B", f"pct_K={rB['pct_K']:g} ; min_(0<=t<=70)IED={rB['ied_min']:.3f}",
      f"IED*={ieB:.4f} (résidu {resB:.0e}) ; IED*-k_env={ieB-kB:+.4f} ; marge min T=70={rB['margin_env_threshold']:.3f} ; certificat d'équilibre P-F1 : Viab(K)≠∅ pour cette instance (T10)",
      "approche progressive du seuil environnemental sans le franchir sur l'horizon",
      "non-vacuité limitée à l'instance B ; ni noyau ni triplets caractérisés ; maintien de la trajectoire simulée sur T=70"],
     ["C", f"sortie de K à t={rC['exit_t']} ; pct_K={rC['pct_K']} ; min_(0<=t<=70)IED={rC['ied_min']:.3f}",
      f"IED*={ieC:.4f} < k_env ; cause={rC['cause']}",
      "sortie du corridor imputable à la contrainte environnementale",
      "dépend des paramètres et de x0 ; une instance"],
     ["D_statique", "IED(1.5,1.5,r), r=" + "/".join(f"{x:g}" for x in stat_r),
      " / ".join(f"{ied_static(x):.4f}" for x in stat_r), "non-compensabilité au bord (propriété de valeur)",
      "algébrique ; ≠ non-substituabilité empirique ; ≠ C3.2 (sensibilité différentielle absolue)"],
     ["D_dynamique", f"IED0 de {dyn[0]['ied0']:.3f} à {dyn[-1]['ied0']:.3f} ; IED_final de {min(r['ied_final'] for r in dyn):.3f} à {max(r['ied_final'] for r in dyn):.3f}",
      "écart IED_final − IED0 > 0", "valeur initiale abaissée puis redressée (pas de chute dynamique)",
      "famille distincte de l'illustration statique"],
     ["E", f"nb_retournements={max(nrev_all)} (τ=3..15) ; max|ΔIED| jusqu'à {gmax:.3f} ; u_unsat={uuns}%",
      "pente/Δt, θ=1e-4 IED/temps ; écart de trajectoire vs τ=0",
      "délai actif sans oscillation observée ; τ>0 hors périmètre EDO de P-F1 (C2/P-F2)",
      "τ∈{0,3,6,10,15}, T=70, ces paramètres ; n'exclut pas l'oscillation en général"],
     ["F", f"départ hors Z_cible ; nb_entrées={rF['z_n_entries']} ; 1re entrée t={rF['z_first_entry_t']} ; "
      f"dernière sortie t={rF['z_last_exit_t']} ; plus long séjour={rF['z_longest_stay']}",
      f"transition entrée : {last_no} < β_entry^Δt(70;p_(−β)) ≤ {first_yes} ; occupation non monotone={dwell_nonmonotone}",
      "accès à Z_cible conditionné par β sur la grille explorée ; l'occupation est décrite sans inférence d'attractivité",
      "transition dépendante du protocole ; ni convergence ni attractivité"],
     ["G_monetaire", f"M(70) de {m_min:.2f} à {m_max:.2f} selon scénario",
      f"M*(T)=[e(T)+R(T)]/λ (cible quasi-stationnaire du WP) ; flux stabilisé ssi flux_rel<1e-3 ; scénarios stabilisés={stabilized} ; borne analytique 0<=M(t)<=max(M0,1/λ)",
      f"M borné pour tout t (borne analytique, propriété (P) vérifiée) ; position de M*(T) / fenêtre de Z_cible — {pos_txt} ; stabilisation de A et F due à la saturation u=1",
      "M*(T) n'est pas un équilibre du système complet ; valeur de référence tant que le flux n'est pas stabilisé"],
    ]
    for row in data:
        print(f"  · {row[0]}")
        for lab, val in zip(("observation", "calcul", "interprét.", "limite"), row[1:]):
            print(f"      {lab:11}: {val}")
    emit("T8_quatre_registres.csv", hdr, data)


def t_admissibilite():
    """T9 — admissibilité de chaque instance au sens de la note P-F1 V1 (§3)."""
    print("\n[T9] ADMISSIBILITÉ AU SENS DE LA NOTE P-F1 (F1, F3′, F4, (L), 𝒫_Φ, (D1)–(D4), (P), régularité)"); bar()
    grid = [i/1000 for i in range(2001)]
    phi_vals = [control_phi(z) for z in grid]
    phi_ok = (control_phi(0.0) == 0.0 and control_phi(1.0) == 1.0
              and all(phi_vals[i+1] >= phi_vals[i] for i in range(len(grid)-1))
              and all(0.0 <= v <= 1.0 for v in phi_vals))
    hdr = ["scenario","perimetre_PF1","F1_continuite","F3prime_via_invariance_forte","F4_beta_pos",
           "L_lipschitz_S","Phi_dans_P_Phi","D1_D4_lambdaM","ebase_r","P_borne_M","u_max","u_unsat_pct","delta=min_IED"]
    rows = []
    for key, sc in SCN.items():
        a, b = sc['alpha'], sc['beta']
        AL = [a, a, a*1.3]; BE = [b, b, b*0.5]
        inv = all(BE[j] + GAMMA[j] <= 2*AL[j] + STATE_TOL and GAMMA[j] >= 0 for j in range(3))
        r = run(sc)
        per = "τ=0 : oui" if sc['tau'] == 0 else "τ>0 : hors périmètre EDO (C2/P-F2)"
        row = [key, per, True, inv, b > 0, True, phi_ok, sc['lambda'] > 0, "e_base≡1 ; r≡0",
               r['M_bound_ok'], 1.0, r['u_unsat'], f"{r['ied_min']:.4f}"]
        rows.append(row)
        print(f"  {LABELS[key][:38]:38} {per:36} F3′={inv} F4={b>0} Φ∈𝒫_Φ={phi_ok} (P)={r['M_bound_ok']} δ={r['ied_min']:.3f}")
    print("  u_max=1 : le régime d'expansion u>1 du WP est exclu de toutes les instances.")
    emit("T9_admissibilite_PF1.csv", hdr, rows)

def t_certificat_B():
    """T10 — certificat d'équilibre de l'instance B, en concordance avec la note P-F1 V1 (§4)."""
    from fractions import Fraction as Fr
    print("\n[T10] CERTIFICAT D'ÉQUILIBRE — INSTANCE B (concordance avec la note P-F1, §4)"); bar()
    sc = SCN["B_approche_seuil"]
    fa, fb, lam = Fr(str(sc['alpha'])), Fr(str(sc['beta'])), Fr(str(sc['lambda']))
    AL = [fa, fa, fa*Fr(13, 10)]; BE = [fb, fb, fb/2]; GA = [Fr(str(g)) for g in GAMMA]
    # polynôme exact Π_i (β_i z + γ_i)/α_i − z³, coefficients [c0, c1, c2, c3]
    poly = [Fr(1)]
    for j in range(3):
        c0, c1 = GA[j]/AL[j], BE[j]/AL[j]
        new = [Fr(0)]*(len(poly)+1)
        for i, c in enumerate(poly):
            new[i] += c*c0; new[i+1] += c*c1
        poly = new
    poly[3] -= 1
    Q_note = [Fr(216), Fr(954), Fr(1395), Fr(-21620)]
    ratio = poly[3]/Q_note[3]
    proportional = all(poly[i] == ratio*Q_note[i] for i in range(4))
    Qf = lambda z: sum(float(c)*z**i for i, c in enumerate(Q_note))
    lo, hi = 0.310, 0.311
    bracket = Qf(lo)*Qf(hi) < 0
    for _ in range(200):
        mid = 0.5*(lo+hi)
        if Qf(lo)*Qf(mid) <= 0: hi = mid
        else: lo = mid
        if hi-lo < 1e-15: break
    z = 0.5*(lo+hi)
    ie, res, Ss = solve_fixed_point(sc['alpha'], sc['beta'])
    Sstar = [(float(BE[j])*z + GAMMA[j])/float(AL[j]) for j in range(3)]
    Mstar = z*E_BASE/sc['lambda']
    ied_S = cbrt(Sstar[0]*Sstar[1]*Sstar[2])
    dS_res = max(abs(-float(AL[j])*Sstar[j] + float(BE[j])*control_phi(ied_S) + GAMMA[j]) for j in range(3))
    dM_res = abs(control_phi(ied_S)*E_BASE - sc['lambda']*Mstar)
    checks = [
      ("polynôme du moteur proportionnel à Q(z) de la note (arithmétique exacte)", proportional),
      ("changement de signe de Q sur (0.310 ; 0.311)", bracket),
      ("z* < 1 (branche non saturée de Φ_sat)", z < 1.0),
      ("|z*(note) − IED*(moteur)| < 1e-12", abs(z - ie) < 1e-12),
      ("S* ∈ (0,2)^3", all(0.0 < s < 2.0 for s in Sstar)),
      ("IED(S*) = z* (|écart| < 1e-12)", abs(ied_S - z) < 1e-12),
      ("M* ∈ (31/15 ; 311/150)", 31/15 < Mstar < 311/150),
      ("IED(S*) > k_env strictement", ied_S > sc['k_env']),
      ("M* > k_social strictement", Mstar > sc['k_social']),
      ("résidus d'équilibre du champ bouclé < 1e-12", dS_res < 1e-12 and dM_res < 1e-12),
    ]
    rows = []
    for lab, ok in checks:
        rows.append([lab, ok]); print(f"  [{'ok' if ok else 'ÉCHEC'}] {lab}")
    rows += [["z*", f"{z:.12f}"], ["S*", " ; ".join(f"{s:.6f}" for s in Sstar)], ["M*", f"{Mstar:.6f}"],
             ["portée", "condition suffisante élémentaire de P-F1 : Viab(K)≠∅ pour l'instance B seulement ; ni noyau ni triplets caractérisés"]]
    print(f"  z*={z:.12f} ; S*=({Sstar[0]:.4f} ; {Sstar[1]:.4f} ; {Sstar[2]:.4f}) ; M*={Mstar:.6f}")
    print("  Portée : Viab(K)≠∅ pour l'instance B seulement ; ni le noyau ni les triplets (F,Φ,Λ) ne sont caractérisés.")
    emit("T10_certificat_B.csv", ["verification", "resultat"], rows)
    if not all(ok for _, ok in checks):
        raise RuntimeError("Certificat de l'instance B : au moins une vérification a échoué.")

# ═════════════════════════════ FIGURES ═════════════════════════════
# Dépendance optionnelle matplotlib. Toutes les figures dynamiques dérivent des
# sorties de run()/solve_fixed_point() ; aucune valeur numérique de résultat
# (instants d'événement, encadrements, valeurs terminales) n'est écrite en dur.

# Palette sobre, compatible avec les principales déficiences de vision des
# couleurs (Okabe–Ito), doublée de styles de lignes distincts pour rester
# interprétable en niveaux de gris.
CB = {"bleu": "#0072B2", "orange": "#E69F00", "vert": "#009E73",
      "vermillon": "#D55E00", "violet": "#CC79A7", "cyan": "#56B4E9",
      "jaune": "#F0E442", "gris": "#555555"}
SCEN_STYLE = {   # (couleur, style de ligne)
    "A_reference":  (CB["bleu"],      "-"),
    "B_approche_seuil":    (CB["orange"],    "--"),
    "C_sortie_K":   (CB["vermillon"], "-."),
    "E_delai":      (CB["vert"],      ":"),
    "F_acces_Zcible": (CB["violet"],    "-"),
}
FIG_FILES = []

def setup_plot_style():
    """Réglages graphiques centralisés : fond blanc, axes nets, polices homogènes."""
    plt.rcdefaults()
    plt.rcParams.update({
        "figure.facecolor": "white", "axes.facecolor": "white",
        "savefig.facecolor": "white", "savefig.bbox": "tight",
        "font.family": "sans-serif", "font.size": 9,
        "axes.titlesize": 10, "axes.labelsize": 9,
        "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 8,
        "axes.linewidth": 0.9, "lines.linewidth": 1.8,
        "axes.spines.top": False, "axes.spines.right": False,
        "grid.linewidth": 0.5, "grid.alpha": 0.35, "grid.color": "#B0B0B0",
        "legend.frameon": True, "legend.framealpha": 0.92,
        "legend.edgecolor": "#B0B0B0", "legend.fancybox": False,
        "mathtext.default": "regular",
    })

def save_figure(fig, basename):
    """Enregistre en PNG 300 dpi et PDF vectoriel de même nom de base, puis ferme."""
    os.makedirs(FIG_DIR, exist_ok=True)
    for ext in ("png", "pdf"):
        path = os.path.join(FIG_DIR, basename + "." + ext)
        if ext == "png":
            fig.savefig(path, dpi=300)
        else:
            # Écrit les métadonnées finales dès la génération : aucune retouche
            # après exécution, pour préserver l'identité octet pour octet.
            fig.savefig(path, metadata={
                "Author": "Benjamin Albouy", "Creator": None, "Producer": None,
                "Title": None, "Subject": None, "Keywords": None,
                "CreationDate": None, "ModDate": None,
            })
        FIG_FILES.append(path)
    plt.close(fig)

# ── Figure 1 — schéma conceptuel : architecture du programme de recherche ──
def fig_architecture_programme():
    fig, ax = plt.subplots(figsize=(9.4, 3.3))
    etapes = ["Working Paper", "Axiomes\n(A1–A5)", "Hypothèses\n(H1–H3)",
              "Question C1\nConjectures\nC2, C4\n+ propriété C3",
              "Note de\ncadrage P-F1",
              "Simulation\nexploratoire\n(démonstrateur)",
              "Investigation\nanalytique\nfuture"]
    n = len(etapes); w, h = 1.95, 0.98; gap = 0.20
    for i, txt in enumerate(etapes):
        x = i*(w+gap)
        courant = (i == 5)
        futur = (i == 6)
        ax.add_patch(Rectangle((x, 0), w, h, facecolor=("#DCE9F5" if courant else "white"),
                               edgecolor=(CB["bleu"] if courant else CB["gris"]),
                               linewidth=(1.6 if courant else 1.0),
                               linestyle=("--" if futur else "-"), zorder=2))
        ax.text(x+w/2, h/2, txt, ha="center", va="center", fontsize=6.4,
                fontweight=("bold" if courant else "normal"), zorder=3)
        if i < n-1:
            ax.add_patch(FancyArrowPatch((x+w, h/2), (x+w+gap, h/2),
                                         arrowstyle="-|>", mutation_scale=11,
                                         color=CB["gris"], linewidth=1.0, zorder=1))
    ax.text(5*(w+gap)+w/2, -0.24, "position du présent article", ha="center",
            va="top", fontsize=7.5, color=CB["bleu"], style="italic")
    ax.set_xlim(-0.3, n*(w+gap)); ax.set_ylim(-0.62, h+0.22)
    ax.axis("off")
    save_figure(fig, "figure_01_architecture_programme")

# ── Figure 2 — schéma conceptuel : structure triangulaire ──
def fig_architecture_triangulaire():
    fig, ax = plt.subplots(figsize=(6.2, 4.0))
    pos = {"S": (0.4, 2.15), "IED": (2.5, 2.15), "u": (4.6, 2.15), "M": (2.5, 0.35)}
    lab = {"S": "$S=(\\mathrm{IBD},\\mathrm{IEE},\\mathrm{IRNR})$\nétat biophysique",
           "IED": "$\\mathrm{IED}(S)$\nmoyenne géométrique",
           "u": "$u=\\min(\\mathrm{IED}(t-\\tau),1)$\ncoefficient\nde régulation",
           "M": "$M$\nmasse monétaire"}
    w, h = 1.85, 0.86
    for k, (x, y) in pos.items():
        ax.add_patch(Rectangle((x-w/2, y-h/2), w, h, facecolor="white",
                               edgecolor=CB["bleu"] if k != "M" else CB["orange"],
                               linewidth=1.3, zorder=2))
        ax.text(x, y, lab[k], ha="center", va="center", fontsize=7.5, zorder=3)
    def arrow(a, b, rad=0.0, color=CB["gris"], ls="-"):
        ax.add_patch(FancyArrowPatch(a, b, arrowstyle="-|>", mutation_scale=12,
                                     color=color, linewidth=1.3, linestyle=ls,
                                     connectionstyle="arc3,rad=%.2f" % rad, zorder=1))
    xS, yS = pos["S"]; xI, yI = pos["IED"]; xU, yU = pos["u"]; xM, yM = pos["M"]
    # chaîne S → IED → u (l'état détermine l'indicateur, l'indicateur le coefficient de régulation)
    arrow((xS+w/2, yS), (xI-w/2, yI))
    arrow((xI+w/2, yI), (xU-w/2, yU))
    # retour u → S : le coefficient de régulation agit sur la dynamique biophysique
    arrow((xU, yU+h/2), (xS, yS+h/2), rad=0.30, color=CB["bleu"])
    ax.text(2.5, 3.16, "$\\beta_i\\,u$ dans $\\dot S_i$ : le coefficient de régulation agit sur $S$",
            ha="center", fontsize=7.5, color=CB["bleu"])
    # alimentation monétaire : IED → u → M ; pas de seconde multiplication par IED
    arrow((xU, yU-h/2), (xM+w/2, yM), rad=0.28, color=CB["orange"])
    ax.text(3.98, 1.16, "$\\dot M=u\\,e_{base}-\\lambda M$ ; $e_{base}\\equiv 1$, $R=0$", fontsize=7.3,
            color=CB["orange"], ha="left")
    # absence de retour M → S : seule flèche manquante de l'architecture
    ax.add_patch(FancyArrowPatch((xM-w/2, yM), (xS, yS-h/2),
                                 arrowstyle="-|>", mutation_scale=12,
                                 color="#C4C4C4", linewidth=1.2, linestyle=(0, (3, 3)),
                                 connectionstyle="arc3,rad=0.30", zorder=1))
    ax.plot([0.72], [1.02], marker="X", markersize=12, color=CB["vermillon"], zorder=4)
    ax.text(0.46, 0.72, "absent : aucun retour\nde $M$ vers $\\dot S$",
            ha="center", va="top", fontsize=7.5, color=CB["vermillon"])
    ax.set_xlim(-1.0, 6.3); ax.set_ylim(-0.55, 3.45); ax.axis("off")
    save_figure(fig, "figure_02_architecture_triangulaire")

# ── Figure 3 — schéma : géométrie de K et de Z_cible (projection (IED, M)) ──
def fig_geometrie_k_zcible():
    sc = SCN["F_acces_Zcible"]          # seuils d'UN scénario, annoncés en légende
    kenv, ksoc = sc['k_env'], sc['k_social']
    zs, mlo, mhi = sc['z_cible'], sc['m_cible_lo'], sc['m_cible_hi']
    xmax, ymax = 1.6, 2.2
    fig, ax = plt.subplots(figsize=(5.4, 4.2))
    ax.add_patch(Rectangle((kenv, ksoc), xmax-kenv, ymax-ksoc, facecolor=CB["bleu"],
                           alpha=0.10, edgecolor="none", zorder=0))
    ax.add_patch(Rectangle((zs, mlo), xmax-zs, mhi-mlo, facecolor=CB["vert"],
                           alpha=0.28, edgecolor=CB["vert"], linewidth=1.3, zorder=1))
    ax.axvline(kenv, color=CB["vermillon"], linestyle="--", linewidth=1.6)
    ax.axhline(ksoc, color=CB["orange"], linestyle="-.", linewidth=1.6)
    ax.axvline(zs, color=CB["vert"], linestyle=":", linewidth=1.4)
    ax.text(kenv+0.012, ymax*0.985, "seuil environnemental\n$\\mathrm{IED}\\geq k_{env}$",
            color=CB["vermillon"], fontsize=7.5, va="top")
    ax.text(xmax*0.985, ksoc+0.03, "plancher social  $M\\geq k_{social}$",
            color=CB["orange"], fontsize=7.5, ha="right", va="bottom")
    ax.text((zs+xmax)/2, (mlo+mhi)/2, "$Z_{cible}$", color="#00543D", fontsize=11,
            ha="center", va="center", fontweight="bold")
    ax.text((kenv+zs)/2, ymax*0.80, "corridor $K$", color=CB["bleu"], fontsize=10,
            ha="center", alpha=0.85)
    ax.text(zs+0.015, mhi+0.03, "$z_{cible}$", color="#00543D", fontsize=8)
    ax.annotate("", xy=(zs-0.045, mlo), xytext=(zs-0.045, mhi),
                arrowprops=dict(arrowstyle="<->", color="#00543D", linewidth=1.1))
    ax.text(zs-0.055, (mlo+mhi)/2, "$[m_{cible,lo},m_{cible,hi}]$", color="#00543D",
            fontsize=7.5, ha="right", va="center", rotation=90)
    ax.set_xlim(0, xmax); ax.set_ylim(0, ymax)
    ax.set_xlabel("$\\mathrm{IED}$"); ax.set_ylabel("$M$ (non tronqué)")
    ax.grid(True, linestyle=":", linewidth=0.4)
    ax.set_title("Projection illustrative dans le plan $(\\mathrm{IED},M)$\n"
                 "seuils du scénario F", fontsize=8.5)
    save_figure(fig, "figure_03_geometrie_K_Zcible")

# ── Figure 4 — C3.1 : illustration algébrique statique ──
def fig_c31_statique():
    """IED(r)=(1.5·1.5·r)^(1/3) : grille continue + points exacts du tableau T4."""
    pts_T4 = [1.5, 1.0, 0.5, 0.2, 0.1, 0.05, 0.01]     # mêmes valeurs que T4
    f = lambda r: (1.5*1.5*r)**(1/3)
    grid = [1e-4 + i*(1.5-1e-4)/1500 for i in range(1501)]
    fig, ax = plt.subplots(figsize=(5.4, 3.8))
    ax.plot(grid, [f(r) for r in grid], color=CB["bleu"], linestyle="-",
            label="$\\mathrm{IED}(r)=(1{,}5\\cdot1{,}5\\cdot r)^{1/3}$")
    ax.plot(pts_T4, [f(r) for r in pts_T4], linestyle="none", marker="o",
            markersize=5, markerfacecolor="white", markeredgecolor=CB["vermillon"],
            markeredgewidth=1.4, label="points du tableau T4")
    for r in (1.5, 0.5, 0.1, 0.01):
        ax.annotate("%.4f" % f(r), xy=(r, f(r)), xytext=(4, -9),
                    textcoords="offset points", fontsize=7, color=CB["vermillon"])
    ax.set_xlim(0, 1.55); ax.set_ylim(0, 1.65)
    ax.set_xlabel("$r=\\mathrm{IRNR}$  ($\\mathrm{IBD}=\\mathrm{IEE}=1{,}5$)")
    ax.set_ylabel("$\\mathrm{IED}$")
    ax.grid(True, linestyle=":", linewidth=0.4)
    ax.legend(loc="lower right")
    save_figure(fig, "figure_04_C31_statique")

# ── Figure 5 — trajectoires IED, scénarios A, B, C, E ──
def fig_trajectoires_ied():
    keys = ["A_reference", "B_approche_seuil", "C_sortie_K", "E_delai"]
    # Deux scénarios peuvent partager la même valeur de k_env (B et C valent
    # 0,30). Les étiquettes seraient alors superposées et l'une masquerait
    # l'autre. On regroupe donc les seuils par valeur distincte et on porte une
    # annotation unique nommant tous les scénarios concernés.
    seuil_groupes = {}
    for _k in keys:
        seuil_groupes.setdefault(round(SCN[_k]['k_env'], 6), []).append(_k.split("_")[0])
    seuils_annotes = set()
    fig, ax = plt.subplots(figsize=(7.6, 4.2))
    for key in keys:
        sc = SCN[key]; r = run(sc, dt=0.05, T=70.0)
        t = [row["t"] for row in r["rows"]]; y = r["ied_traj"]
        col, ls = SCEN_STYLE[key]
        ax.plot(t, y, color=col, linestyle=ls,
                label="%s\n(seuil propre $k_{env}=%.2f$)" % (LABELS[key], sc['k_env']))
        # Seuil environnemental PROPRE à chaque scénario : repère interne au
        # domaine [0,70], tracé dans le style du scénario et annoté par son nom.
        # Jamais une ligne horizontale unique (qui suggérerait un seuil commun),
        # jamais un élément placé hors du domaine de définition des trajectoires.
        ax.plot([0.0, 70.0], [sc['k_env']]*2, color=col, linestyle=ls,
                linewidth=0.9, alpha=0.45, zorder=1)
        _kv = round(sc['k_env'], 6)
        if _kv not in seuils_annotes:
            seuils_annotes.add(_kv)
            ax.annotate("$k_{env}$ %s = %.2f" % (" = ".join(seuil_groupes[_kv]), sc['k_env']),
                        xy=(69.4, sc['k_env']), xytext=(0, 3), textcoords="offset points",
                        fontsize=6.8, color=col, ha="right", va="bottom", zorder=6,
                        bbox=dict(boxstyle="round,pad=0.14", facecolor="white",
                                  edgecolor="none", alpha=0.85))
        # marqueur discret de sortie de K : instant issu des métriques du moteur
        if r["exit_t"] is not None:
            k = min(range(len(t)), key=lambda i: abs(t[i]-r["exit_t"]))
            ax.plot([t[k]], [y[k]], marker="v", markersize=8, color=col, zorder=5)
            ax.annotate("sortie de $K$ (%s) : $t=%.2f$, cause « %s »"
                        % (key.split("_")[0], r["exit_t"], r["cause"]),
                        xy=(t[k], y[k]), xytext=(-16, -26), textcoords="offset points",
                        fontsize=7, color=col, ha="center",
                        arrowprops=dict(arrowstyle="-", color=col, linewidth=0.7))
    # domaine temporel strictement égal au domaine de définition des trajectoires
    ax.set_xlim(0, 70); ax.set_ylim(0, None)
    ax.set_xlabel("$t$ (unités de temps)"); ax.set_ylabel("$\\mathrm{IED}(t)$")
    ax.grid(True, linestyle=":", linewidth=0.4)
    ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0), fontsize=6.8,
              borderaxespad=0.0)
    fig.tight_layout()
    save_figure(fig, "figure_05_trajectoires_IED")

# ── Figure 6 — portraits projetés (IED, M) : A, C, F ──
def fig_portraits_ied_m():
    keys = ["A_reference", "C_sortie_K", "F_acces_Zcible"]
    fig, axes = plt.subplots(1, 3, figsize=(10.6, 3.9))
    for ax, key in zip(axes, keys):
        sc = SCN[key]; r = run(sc, dt=0.05, T=70.0)
        x = r["ied_traj"]; M = [row["M"] for row in r["rows"]]   # M non tronqué
        t = [row["t"] for row in r["rows"]]
        col, _ = SCEN_STYLE[key]
        ax.plot(x, M, color=col, linewidth=1.6, zorder=3)
        # cadrage : trajectoire ET seuils propres au scénario, afin qu'aucun
        # panneau ne masque la position relative de la trajectoire et du corridor.
        # Les échelles restent propres à chaque panneau (pas d'échelle commune).
        xlo = min(min(x), sc['k_env']); xhi = max(max(x), sc['z_cible'])
        ylo = min(min(M), sc['k_social']); yhi = max(max(M), sc['m_cible_hi'])
        px = 0.10*(xhi-xlo or 1); py = 0.10*(yhi-ylo or 1)
        X0, X1, Y0, Y1 = xlo-px, xhi+px, max(0.0, ylo-py), yhi+py
        # zone K et bornes propres AU SCÉNARIO
        ax.add_patch(Rectangle((sc['k_env'], sc['k_social']), max(0, X1-sc['k_env']),
                               max(0, Y1-sc['k_social']), facecolor=CB["bleu"],
                               alpha=0.07, edgecolor="none", zorder=0))
        ax.axvline(sc['k_env'], color=CB["vermillon"], linestyle="--", linewidth=1.2)
        ax.axhline(sc['k_social'], color=CB["orange"], linestyle="-.", linewidth=1.2)
        ax.axvline(sc['z_cible'], color=CB["vert"], linestyle=":", linewidth=1.1)
        ax.add_patch(Rectangle((sc['z_cible'], sc['m_cible_lo']), max(0, X1-sc['z_cible']),
                               sc['m_cible_hi']-sc['m_cible_lo'], facecolor=CB["vert"],
                               alpha=0.20, edgecolor=CB["vert"], linewidth=1.0, zorder=1))
        ax.plot([x[0]], [M[0]], marker="o", markersize=7, markerfacecolor="white",
                markeredgecolor=col, markeredgewidth=1.5, zorder=6, label="$x_0$")
        ax.plot([x[-1]], [M[-1]], marker="s", markersize=6, color=col, zorder=6,
                label="$x(T)$")
        if r["exit_t"] is not None:                 # instant issu du moteur
            k = min(range(len(t)), key=lambda i: abs(t[i]-r["exit_t"]))
            ax.plot([x[k]], [M[k]], marker="v", markersize=8, color=CB["vermillon"],
                    zorder=7, label="sortie de $K$, $t=%.2f$" % r["exit_t"])
        if r["z_first_entry_t"] is not None:        # instant issu du moteur
            k = min(range(len(t)), key=lambda i: abs(t[i]-r["z_first_entry_t"]))
            ax.plot([x[k]], [M[k]], marker="*", markersize=12, color="#00543D",
                    zorder=7, label="1re entrée $Z_{cible}$, $t=%.2f$" % r["z_first_entry_t"])
        ax.set_xlim(X0, X1); ax.set_ylim(Y0, Y1)     # échelles propres, non communes
        ax.set_xlabel("$\\mathrm{IED}$")
        ax.set_ylabel("$M$ (non tronqué)")
        ax.set_title(LABELS[key].replace("Z_cible", "$Z_{cible}$"), fontsize=8)
        ax.grid(True, linestyle=":", linewidth=0.4)
        ax.legend(loc="best", fontsize=6.8)
        # Encart agrandi lorsque la fenêtre monétaire de Z_cible occupe une part
        # trop faible de l'ordonnée pour être lisible. La trajectoire complète
        # reste tracée dans le panneau principal, à échelle propre et sans
        # troncature de M ; l'encart n'est qu'un agrandissement d'une portion
        # du même tracé, aux mêmes données.
        span = (Y1-Y0) or 1.0
        if (sc['m_cible_hi']-sc['m_cible_lo'])/span < 0.20:
            iax = ax.inset_axes([0.12, 0.40, 0.37, 0.34])
            iy0 = max(0.0, min(sc['k_social'], sc['m_cible_lo']) - 0.18)
            iy1 = sc['m_cible_hi'] + 0.18
            iax.add_patch(Rectangle((sc['k_env'], sc['k_social']),
                                    max(0, X1-sc['k_env']), max(0, iy1-sc['k_social']),
                                    facecolor=CB["bleu"], alpha=0.07, edgecolor="none",
                                    zorder=0))
            iax.add_patch(Rectangle((sc['z_cible'], sc['m_cible_lo']),
                                    max(0, X1-sc['z_cible']),
                                    sc['m_cible_hi']-sc['m_cible_lo'], facecolor=CB["vert"],
                                    alpha=0.20, edgecolor=CB["vert"], linewidth=1.0,
                                    zorder=1))
            iax.axvline(sc['k_env'], color=CB["vermillon"], linestyle="--", linewidth=1.0)
            iax.axhline(sc['k_social'], color=CB["orange"], linestyle="-.", linewidth=1.0)
            iax.axvline(sc['z_cible'], color=CB["vert"], linestyle=":", linewidth=1.0)
            iax.plot(x, M, color=col, linewidth=1.3, zorder=3)
            iax.plot([x[0]], [M[0]], marker="o", markersize=5, markerfacecolor="white",
                     markeredgecolor=col, markeredgewidth=1.3, zorder=5)
            iax.set_xlim(X0, X1); iax.set_ylim(iy0, iy1)
            iax.tick_params(labelsize=5.5)
            iax.set_title("agrandissement : $M\\in[%.2f\\,;\\,%.2f]$" % (iy0, iy1),
                          fontsize=5.8, pad=2)
            for sp in iax.spines.values():
                sp.set_visible(True); sp.set_linewidth(0.7); sp.set_color("#888888")
    fig.tight_layout()
    save_figure(fig, "figure_06_portraits_IED_M")

# ── Figures complémentaires (annexe) ──
def fig_delai_complementaire():
    """S1 — effet du délai τ sur IED(t), scénario E. Complète T6 : le tableau ne
    donne que max|ΔIED| et l'état d'oscillation, pas la forme de l'écart."""
    base = SCN["E_delai"]
    taus = [0, 3, 6, 10, 15]
    cols = [CB["gris"], CB["bleu"], CB["cyan"], CB["orange"], CB["vermillon"]]
    lss = ["-", "--", "-.", ":", (0, (5, 1, 1, 1))]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.4, 3.7))
    b0 = dict(base); b0['tau'] = 0; r0 = run(b0, 0.05); t0 = [w["t"] for w in r0["rows"]]
    for tau, c, ls in zip(taus, cols, lss):
        bt = dict(base); bt['tau'] = tau; rt = run(bt, 0.05)
        t = [w["t"] for w in rt["rows"]]
        ax1.plot(t, rt["ied_traj"], color=c, linestyle=ls, label="$\\tau=%d$" % tau)
        n = min(len(t), len(t0))
        ax2.plot(t[:n], [rt["ied_traj"][k]-r0["ied_traj"][k] for k in range(n)],
                 color=c, linestyle=ls, label="$\\tau=%d$" % tau)
    # seuil environnemental propre au scénario E : repère interne au domaine
    # [0,70], annoté ; aucun élément graphique hors du domaine de définition
    ax1.plot([0.0, 70.0], [base['k_env']]*2, color=CB["vermillon"], linestyle="--",
             linewidth=0.9, alpha=0.55, zorder=1)
    ax1.annotate("$k_{env}$ (E) = %.2f" % base['k_env'], xy=(69.4, base['k_env']),
                 xytext=(0, 3), textcoords="offset points", fontsize=6.8,
                 color=CB["vermillon"], ha="right", va="bottom",
                 bbox=dict(boxstyle="round,pad=0.14", facecolor="white",
                           edgecolor="none", alpha=0.85))
    ax1.set_xlabel("$t$"); ax1.set_ylabel("$\\mathrm{IED}(t)$")
    ax2.set_xlabel("$t$"); ax2.set_ylabel("$\\mathrm{IED}_\\tau(t)-\\mathrm{IED}_0(t)$")
    for a in (ax1, ax2):
        a.grid(True, linestyle=":", linewidth=0.4); a.legend(fontsize=7)
        a.set_xlim(0, 70)          # domaine de définition des trajectoires
    ax2.axhline(0.0, color="#999999", linewidth=0.8)
    fig.tight_layout()
    save_figure(fig, "figure_S1_delai_tau")

def fig_beta_zcible_complementaire():
    """S2 — β vs premier instant d'entrée et durée d'occupation, scénario F.
    Les β sans entrée (None) sont ABSENTS de la courbe et signalés par une
    convention explicite en bas de panneau."""
    sc = SCN["F_acces_Zcible"]
    betas = list(BETA_GRID)
    first, dwell, none_b = [], [], []
    for b in betas:
        p = dict(sc); p['beta'] = b; r = run(p, 0.05)
        dwell.append(r["z_longest_stay"])
        if r["z_first_entry_t"] is None: none_b.append(b)
        else: first.append((b, r["z_first_entry_t"]))
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.4, 3.7))
    ax1.plot([b for b, _ in first], [v for _, v in first], color=CB["bleu"],
             linestyle="-", marker="o", markersize=4)
    if none_b:
        ax1.axvspan(min(none_b)-0.005, max(none_b)+0.005, color=CB["vermillon"],
                    alpha=0.06)
        ax1.text(0.5*(min(none_b)+max(none_b)), 4.0,
                 "aucune entrée sur $[0,T]$\n(temps d'entrée indéfini)",
                 fontsize=6.8, color=CB["vermillon"], ha="center", va="center")
        # encadrement obtenu du moteur, non écrit en dur
        last_no = max(none_b); first_yes = min(b for b, _ in first)
        ax1.axvline(0.5*(last_no+first_yes), color=CB["gris"], linestyle=":",
                    linewidth=1.1)
        ax1.text(0.5*(last_no+first_yes)+0.004, 62,
                 "transition numérique :\n$%.2f<\\beta_{entry}^{\\Delta t}(70;p_{-\\beta})\\leq%.2f$\n"
                 "(dépendante du protocole)" % (last_no, first_yes),
                 fontsize=6.8, color=CB["gris"], va="top")
    ax1.set_ylim(0, 72)
    ax1.set_xlabel("$\\beta$"); ax1.set_ylabel("premier instant d'entrée dans $Z_{cible}$")
    ax2.plot(betas, dwell, color=CB["vert"], linestyle="-", marker="s", markersize=4)
    pk = max(range(len(betas)), key=lambda i: dwell[i])
    ax2.annotate("maximum sur la grille\n($\\beta=%.2f$)" % betas[pk],
                 xy=(betas[pk], dwell[pk]), xytext=(0, -34),
                 textcoords="offset points", fontsize=6.8, ha="center",
                 color="#00543D",
                 arrowprops=dict(arrowstyle="-|>", color="#00543D", linewidth=0.8))
    ax2.text(0.03, 0.96, "occupation mesurée sur la grille admissible ;\naucune attractivité n'est inférée",
             transform=ax2.transAxes, fontsize=6.8, ha="left", va="top", color=CB["gris"])
    ax2.set_xlabel("$\\beta$"); ax2.set_ylabel("plus long séjour continu dans $Z_{cible}$")
    for a in (ax1, ax2):
        a.grid(True, linestyle=":", linewidth=0.4)
    fig.tight_layout()
    save_figure(fig, "figure_S2_beta_Zcible")

def generate_all_figures():
    """Génère toutes les figures si matplotlib est disponible ; sinon signale
    explicitement l'omission sans interrompre le protocole."""
    print("\n[FIGURES] génération graphique"); bar()
    if not HAS_MPL:
        print("  matplotlib INDISPONIBLE (%s)." % MPL_ERR)
        print("  Les figures n'ont PAS été générées. Tableaux et CSV sont complets et inchangés.")
        print("  Installation : pip install matplotlib")
        return False
    setup_plot_style()
    del FIG_FILES[:]
    for fn in (fig_architecture_programme, fig_architecture_triangulaire,
               fig_geometrie_k_zcible, fig_c31_statique, fig_trajectoires_ied,
               fig_portraits_ied_m, fig_delai_complementaire,
               fig_beta_zcible_complementaire):
        fn()
    ok = True
    for path in FIG_FILES:
        exists = os.path.exists(path) and os.path.getsize(path) > 0
        ok = ok and exists
        print("  [FIG] figures/%s%s" % (os.path.basename(path),
                                        "" if exists else "   *** MANQUANT OU VIDE ***"))
    return ok

if __name__ == "__main__":
    print("Démonstrateur EH — moteur de référence, sorties reproductibles")
    t_robustesse(); t_monetaire(); t_pointsfixes()
    t_c31_statique(); t_c31_dynamique(); t_delai(); t_beta(); t_registres()
    t_admissibilite(); t_certificat_B()
    figs_ok = generate_all_figures()
    print("\nFIN — tous les tableaux régénérés en une exécution.")
    print("  CSV     : ./csv/")
    if HAS_MPL:
        print("  Figures : ./figures/  (PNG 300 dpi + PDF vectoriel) — génération %s"
              % ("réussie" if figs_ok else "INCOMPLÈTE : voir les avertissements ci-dessus"))
    else:
        print("  Figures : non générées (matplotlib absent) — les tableaux et les CSV")
        print("            ci-dessus sont complets et strictement inchangés.")
