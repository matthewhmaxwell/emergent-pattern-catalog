"""N13 plug-in: all-to-all noisy voter model with population-wide shocks (arXiv 2609.39858, v1).

Paper model ("Exogenous shocks in the voter model" and "Exogenous shocks in the noisy voter model"). N agents with
opinions s_i in {-1, +1} on the complete graph; magnetisation m = (1/N) sum_i s_i (eq. 1), n = number holding +1,
m = 2n/N - 1. Discrete time, Delta t = 1/N per elementary step. Each step:
  * with probability r an exogenous shock drives the system instantaneously to the reference magnetisation m_r
    ("delta shocks", eq. 10: P_r = (1 - r) P + r R, R sending every state to n_r);
  * otherwise (probability 1 - r) one noisy-voter update: a random agent, with probability eta, changes its opinion
    regardless of the others; with probability 1 - eta it copies the opinion of a random agent (eq. 27:
    p_up = (1 - theta)[eta + (1 - eta) theta], p_down = theta[eta + (1 - eta)(1 - theta)], theta = n/N).
Without shocks the stationary distribution changes from bimodal to unimodal at eta_c = 1/N. Fig. 6 (N = 500,
m_0 = m_r = 0.3) maps the (r, eta) plane; the trimodal state (peaks at m = -1, m_r, +1) is panel d), whose marker
sits at eta = 10^-3.5 = 3.16e-4 (= 0.158/N) and r = 10^-4.5 = 3.16e-5 (read from the log axes of the phase diagram:
the markers a)-f) lie on half-decade grid values, e.g. b) at eta = 1e-2, r = 1e-4; the text gives no numbers).

Headline (plan row N13: the paper's trimodal setting, keeping the paper's scaling in N):
  eta = 0.158/N (eta N = 10^-3.5 * 500 = 0.1581; census N = 200: eta = 7.91e-4; area4 N = 800: 1.98e-4);
  r = 10^-4.5 (500/N)^2 per elementary step (r N^2 = 7.91 fixed; census N = 200: 1.98e-4; area4 N = 800: 1.24e-5).
  Why r scales as 1/N^2: in the continuous-time description (Appendix, eq. 32) the shock rate per unit time is N r
  (Delta t = 1/N per step) and the stationary solutions depend on it only through N times that rate, i.e. the
  factors N r / ((n+1)(n+2) + N r) in eq. (34) and N r / (2 + N r) in eq. (35) with r per unit time, which is
  N^2 * (per-step r); the noise enters through eta N (eta_c = 1/N). Keeping eta N and N^2 r fixed keeps the
  paper's regime at the census size (consensus mass 4.2% at N = 200 vs 3.6% at N = 500).
  m_r = 0.3 (65% hold +1), start m_0 = m_r (as in Fig. 6), opinions +1 / -1 stored as node types 1 / 0.
  The exact stationary distribution of these rules is trimodal at N = 200, 500 and 800 (checked with the chain and
  the binned-prominence criterion below).

Census mapping and choices:
  * one census step = one sweep = N elementary steps of the paper's chain (one unit of the paper's time, Delta t
    = 1/N per step); each step is a shock with probability r, otherwise one agent update;
  * the agent update copies a randomly chosen OTHER agent (plan rule 2, card text; the paper's eq. 27 uses theta =
    n/N, i.e. self included; the difference is a factor N/(N - 1) in the copy rate, O(1/N));
  * "flips its opinion at random" (card) = the chosen agent changes its opinion with probability eta (eq. 27);
  * shock: the paper fixes only the macroscopic jump (on the complete graph it does not matter which agents change);
    following the card ("opinions are instantly reassigned"), a uniformly random set of exactly n_r = 0.65 N agents
    is given +1 and the rest -1;
  * output: "nt" (S, frames, NN) int8 opinions in {0, 1}, frames every 2 sweeps; "adj0" the complete graph without
    self-loops; no phases;
  * knock-out (plan): rule 2 (copying) removed; noise flips and shocks are kept (an update that would have copied
    leaves the agent unchanged).

Fidelity check (plan: three peaks of the stationary mean-opinion distribution, one near 0.3, at the paper's N = 500,
eta = 3.16e-4, r = 3.16e-5, m_r = 0.3). Tolerance (fixed before running), pass = both of:
  (a) the exact stationary PMF of the census rules' Markov chain and of the paper's eq. 27 chain, binned as in (b),
      has exactly three peaks of relative prominence >= 0.1, at the two consensus states and in the bin of m_r
      (state-level local maxima are reported too, as context);
  (b) a long simulation with this module's own agent-level update code (M independent replicas from m_0 = 0.3,
      burn-in, then m sampled every sweep) gives a binned density (consensus states as their own cells, interior
      states in bins of 25 states centred on n_r) with exactly three peaks of relative prominence >= 0.1, at the two
      consensus states and in the bin of m_r, and a Kolmogorov-Smirnov distance <= 0.02 between the simulated and
      the exact state-level distributions.

presence(): the same distribution over a census run (all frames of all seeds): binned density, peaks, shares of time
at consensus and near m_r, next to the exact stationary distribution at the census N for the same rules.
"""
import time

import numpy as np

from . import base

CID = "planted:N13"
LAYERS = "N"
K = {"N": 2}
SPEED = 0.0
RULES = [("NOISE_FLIP", False),            # 1: chosen agent changes opinion with probability eta (one-way)
         ("COPY_RANDOM_OTHER", True),      # 2: otherwise copies a random other agent (I)
         ("POPULATION_SHOCK", False)]      # 3: with probability r per step, opinions reset to m = m_r (external)
ETA_N = 10 ** -3.5 * 500                  # eta * N at the trimodal marker d) of Fig. 6 (eta = 10^-3.5 at N = 500)
HEADLINE = {
    "eta_times_N": ETA_N,                 # 0.1581; eta = ETA_N / NN (census 7.91e-4, area4 1.98e-4)
    "eta": ETA_N / 200,
    "r_times_N2": 10 ** -4.5 * 500 ** 2,  # 7.91; r = 10^-4.5 per step at N = 500 (Fig. 6 marker d)), r ~ 1/N^2
    "r": 10 ** -4.5 * (500 / 200) ** 2,   # 1.98e-4 per elementary step at NN = 200 (area4: 1.24e-5)
    "m_r": 0.3,                           # Fig. 6 (m_0 = m_r = 0.3)
    "NN": 200,
}
STEP = "one sweep = NN elementary steps of the paper's chain (Delta t = 1/N each): shock with probability r, else one agent update"
RULE_TEXT = ("World: 200 agents, each holding opinion +1 or -1, every agent able to interact with every other. Rules "
             "each step: pick an agent; with probability eta = 0.00079 (0.16 times the inverse population) it flips "
             "its opinion at random; otherwise it copies the opinion of a randomly chosen other agent. Separately, "
             "with probability r = 0.000198 per step, a population-wide shock occurs: opinions are instantly "
             "reassigned so that the average opinion jumps to 0.3 (65% holding +1). Start with 65% holding +1 and "
             "run 1000 sweeps of 200 steps. Measure the long-run distribution of the average opinion.")

PAPER_N = 500
PAPER_ETA = 10 ** -3.5
PAPER_R = 10 ** -4.5


def _params(NN):
    """eta, per-step r and n_r at population NN (eta N and r N^2 fixed at the paper's N = 500 values)."""
    return ETA_N / NN, HEADLINE["r_times_N2"] / NN ** 2, int(round((1 + HEADLINE["m_r"]) / 2 * NN))


# ------------------------------------------------------------------------------------------------------ dynamics
def _draws(gens, S, N):
    if isinstance(gens, np.random.Generator):
        return gens.integers(0, N, (S, N)), gens.random((S, N)), gens.integers(0, N - 1, (S, N))
    I = np.stack([g.integers(0, N, N) for g in gens]); U = np.stack([g.random(N) for g in gens])
    Jo = np.stack([g.integers(0, N - 1, N) for g in gens])
    return I, U, Jo


def _reset(s, row, nr, g):
    N = s.shape[1]
    s[row] = 0; s[row, g.permutation(N)[:nr]] = 1


def _run(s, sweeps, eta, r, nr, copy_on, gens, rec_every=None, rec_n=False):
    """agent-level dynamics on (S, N) int8 opinions (1 = +1, 0 = -1), in place. gens: list (one per row) or one
    Generator for all rows. Returns (frames list, n per sweep array or None)."""
    S, N = s.shape; ar = np.arange(S)
    pf = r + (1.0 - r) * eta                      # u < r: shock; r <= u < pf: noise flip; else copy (or nothing)
    frames = [s.copy()] if rec_every else []
    ns = np.empty((sweeps, S), np.int32) if rec_n else None
    for t in range(1, sweeps + 1):
        I, U, Jo = _draws(gens, S, N)
        Jo = Jo + (Jo >= I)                       # a uniformly random OTHER agent
        F = U < pf
        shocks = np.argwhere(U < r)               # (row, step) pairs, rare
        shock_at = {}
        for row, k in shocks: shock_at.setdefault(int(k), []).append(int(row))
        for k in range(N):
            i = I[:, k]
            cur = s[ar, i]
            if copy_on:
                s[ar, i] = np.where(F[:, k], 1 - cur, s[ar, Jo[:, k]])
            else:
                s[ar, i] = np.where(F[:, k], 1 - cur, cur)
            if k in shock_at:                     # the shock replaces this step's update in those rows
                for row in shock_at[k]:
                    _reset(s, row, nr, gens if isinstance(gens, np.random.Generator) else gens[row])
        if rec_every and t % rec_every == 0: frames.append(s.copy())
        if rec_n: ns[t - 1] = s.sum(1)
    return frames, ns


def _init(S, N, nr, gens):
    s = np.zeros((S, N), np.int8)
    for row in range(S): _reset(s, row, nr, gens if isinstance(gens, np.random.Generator) else gens[row])
    return s


def simulate(seed0, steps=1000, S=3, size="census", keep=None):
    on = base.kept(keep, RULES)
    z = base.SIZES[size]; NN, W = z["NN"], z["W"]
    eta, r, nr = _params(NN)
    gens = base.rngs(seed0, S)
    s = _init(S, NN, nr, gens)                    # m_0 = m_r (Fig. 6)
    frames, _ = _run(s, steps, eta, r, nr, 2 in on, gens, rec_every=base.REC["nt"])
    nt = np.stack(frames, 1).astype(np.int8)
    adj0 = np.broadcast_to(~np.eye(NN, dtype=bool), (S, NN, NN)).copy()
    return {"nt": nt, "adj0": adj0, "meta": base.meta(steps, W)}


# ---------------------------------------------------------------------------------------------- exact chain, peaks
def stationary(N, eta, r, nr, other=True):
    """exact stationary PMF over n = 0..N of the discrete-time chain P_r = (1 - r) P + r R (eqs. 10, 27-28).
    other=True: copy target is a random other agent (census rule); False: theta = n/N as in eq. (27)."""
    n = np.arange(N + 1); th = n / N
    cp = (n / (N - 1), (N - n) / (N - 1)) if other else (th, 1 - th)
    up = (1 - th) * (eta + (1 - eta) * cp[0])
    dn = th * (eta + (1 - eta) * cp[1])
    P = np.zeros((N + 1, N + 1))                  # row-stochastic P[i, j] = Pr(i -> j)
    P[n[:-1], n[:-1] + 1] = up[:-1]; P[n[1:], n[1:] - 1] = dn[1:]; P[n, n] += 1 - up - dn
    A = (np.eye(N + 1) - (1 - r) * P).T
    b = np.zeros(N + 1); b[nr] = r
    pi = np.linalg.solve(A, b)
    return pi / pi.sum()


def _bins(N, nr, w=None):
    """cells: [0], interior bins of w states centred on nr (partial bins at the edges), [N]."""
    if w is None: w = 2 * int(round(N / 40)) + 1
    h = w // 2
    lo = nr - h - w * ((nr - h - 1) // w)         # first full-bin start >= 1
    edges = [1] + list(range(lo, N, w)) if lo > 1 else list(range(lo, N, w))
    cells = [(0, 0)]
    for a, b in zip(edges, edges[1:] + [N]):
        if b - 1 >= a: cells.append((a, min(b - 1, N - 1)))
    cells.append((N, N))
    return cells


def _binned(pmf_states, N, nr, w=None):
    cells = _bins(N, nr, w)
    dens = np.array([pmf_states[a:b + 1].sum() / (b - a + 1) for a, b in cells])
    mid = np.array([(2 * (a + b) / 2.0) / N - 1 for a, b in cells])
    return mid, dens, cells


def _peaks(d, rel_prom=0.0):
    """strict local maxima of a 1D sequence whose topographic prominence (height above the higher of the lowest
    points between it and the nearest higher point on each side; one-sided at the ends) is >= rel_prom * height."""
    d = np.asarray(d, float); n = len(d); out = []
    for i in range(n):
        l = d[i - 1] if i > 0 else -np.inf; rr = d[i + 1] if i < n - 1 else -np.inf
        if not (d[i] > l and d[i] > rr): continue
        bases = []
        for step in (-1, 1):
            j = i + step
            if not 0 <= j < n: continue
            mn = d[j]
            while 0 <= j < n and d[j] <= d[i]:
                mn = min(mn, d[j]); j += step
            bases.append(mn)
        prom = d[i] - max(bases) if bases else d[i]
        if prom >= rel_prom * d[i]: out.append(i)
    return out


def _state_peaks(pi):
    N = len(pi) - 1
    return [k for k in range(N + 1) if pi[k] > (pi[k - 1] if k > 0 else -1) and pi[k] > (pi[k + 1] if k < N else -1)]


def fidelity(budget_s=7200, M=2000, burn=1000, sweeps=4000):
    t0 = time.time()
    N, eta, r = PAPER_N, PAPER_ETA, PAPER_R
    nr = int(round((1 + HEADLINE["m_r"]) / 2 * N))
    res = {"check": "peaks of the stationary distribution of m at the trimodal setting (Fig. 6 d)",
           "setting": {"N": N, "eta": eta, "r": r, "m_r": HEADLINE["m_r"], "m_0": HEADLINE["m_r"]},
           "tolerance": {"a": "exact PMFs (census rules, eq. 27), binned as in b: exactly 3 peaks with relative "
                              "prominence >= 0.1 at -1, +1 and the m_r bin",
                         "b": "simulated binned density (25-state bins + consensus cells): exactly 3 peaks with "
                              "relative prominence >= 0.1 at -1, +1 and the m_r bin; KS distance to exact <= 0.02"},
           "paper": {"peaks": 3, "middle_peak": 0.3}}
    # (a) exact chains
    pi_c = stationary(N, eta, r, nr, other=True); pi_p = stationary(N, eta, r, nr, other=False)
    pk_c = [2 * k / N - 1 for k in _state_peaks(pi_c)]; pk_p = [2 * k / N - 1 for k in _state_peaks(pi_p)]
    ex = {}
    ok_a = True
    for tag, pi_ in (("census_rules", pi_c), ("paper_eq27", pi_p)):
        mid_x, dens_x, cells_x = _binned(pi_, N, nr, w=25)
        pkx = _peaks(dens_x, rel_prom=0.1)
        ib = next(j for j, (a, b) in enumerate(cells_x) if a <= nr <= b)
        okx = len(pkx) == 3 and 0 in pkx and len(cells_x) - 1 in pkx and ib in pkx
        ex[tag] = {"binned_peaks_m": [float(mid_x[i]) for i in pkx], "pass": bool(okx)}
        ok_a &= okx
    res["exact"] = {"binned": ex, "census_rules_state_peaks_m": pk_c, "paper_eq27_state_peaks_m": pk_p,
                    "pmf_at_peaks": [float(pi_c[int(round((m + 1) * N / 2))]) for m in pk_c],
                    "consensus_mass": float(pi_c[0] + pi_c[-1])}
    # (b) agent-level simulation with this module's update code; replicas scaled to the budget
    g = np.random.default_rng([base.seed_of(CID), 39858])
    tp = time.time(); s = _init(200, N, nr, g); _run(s, 5, eta, r, nr, True, g)
    per = (time.time() - tp) / 5 / 200                   # seconds per replica-sweep
    fit = int(0.85 * (budget_s - (time.time() - t0)) / max(per * (burn + sweeps), 1e-12))
    M = int(max(50, min(M, fit)))
    s = _init(M, N, nr, g)
    _run(s, burn, eta, r, nr, True, g)
    _, ns = _run(s, sweeps, eta, r, nr, True, g, rec_n=True)
    cnt = np.bincount(ns.ravel(), minlength=N + 1).astype(float); pmf_sim = cnt / cnt.sum()
    ks = float(np.abs(np.cumsum(pmf_sim) - np.cumsum(pi_c)).max())
    mid, dens, cells = _binned(pmf_sim, N, nr, w=25)
    pk = _peaks(dens, rel_prom=0.1)
    pk_m = [float(mid[i]) for i in pk]
    ibin = next(j for j, (a, b) in enumerate(cells) if a <= nr <= b)
    ok_b = (len(pk) == 3 and 0 in pk and len(cells) - 1 in pk and ibin in pk and ks <= 0.02)
    mid_e, dens_e, _ = _binned(pi_c, N, nr, w=25)
    res["simulation"] = {"replicas": M, "burn_sweeps": burn, "sample_sweeps": sweeps, "ks_distance": ks,
                         "binned_peaks_m": pk_m, "consensus_mass": float(pmf_sim[0] + pmf_sim[-1]),
                         "pmf_sim_at_exact_peaks": [float(pmf_sim[int(round((m + 1) * N / 2))]) for m in pk_c],
                         "binned_density_sim": dens.round(6).tolist(), "binned_density_exact": dens_e.round(6).tolist(),
                         "bin_centres_m": mid.round(4).tolist()}
    res.update({"a_pass": bool(ok_a), "b_pass": bool(ok_b), "pass": bool(ok_a and ok_b),
                "secs": round(time.time() - t0, 1)})
    return res


# -------------------------------------------------------------------------------------------------------- presence
def presence(out):
    nt = np.asarray(out["nt"]); S, F, N = nt.shape
    eta, r, nr = _params(N)
    n = nt.astype(np.int32).sum(-1)                         # (S, F)
    cnt = np.bincount(n.ravel(), minlength=N + 1).astype(float); pmf = cnt / cnt.sum()
    mid, dens, cells = _binned(pmf, N, nr)
    pk = _peaks(dens, rel_prom=0.1)
    pi = stationary(N, eta, r, nr, other=True)
    mid_e, dens_e, _ = _binned(pi, N, nr)
    m = 2 * n / N - 1
    return {"measure": "distribution of the mean opinion m over the census run (all frames, all seeds)",
            "N": int(N), "frames": int(F), "binned_peaks_m": [float(mid[i]) for i in pk],
            "n_peaks": len(pk), "share_consensus": float(np.mean(np.abs(m) == 1.0)),
            "share_near_m_r": float(np.mean(np.abs(m - HEADLINE["m_r"]) <= 0.05)),
            "per_seed_mean_m": m.mean(1).round(4).tolist(),
            "binned_density": dens.round(5).tolist(), "bin_centres_m": mid.round(3).tolist(),
            "exact_stationary_peaks_m": [2 * k / N - 1 for k in _state_peaks(pi)],
            "exact_binned_density": dens_e.round(5).tolist(),
            "exact_consensus_mass": float(pi[0] + pi[-1]),
            "expected_shocks_per_seed": float(r * N * (F - 1) * out["meta"]["rec"]["nt"]),
            "note": "one census run is short compared with the stationary mixing of the distribution"}
