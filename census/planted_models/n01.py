"""N01 plug-in: all-to-all pulse-coupled oscillators at linear charging (arXiv 2609.01668, v2).

Paper model (Sec. 2, eq. 1): N identical oscillators on a complete graph, voltage x_i in [0, 1],
    dx_i/dt = 1 - gamma * x_i,
studied at the critical, linear charging curve gamma = 0 (Sec. 3), so x_i(t) = x_i(0) + t between firings. When
x_i = 1 the unit fires, resets to 0 and delivers a pulse Delta = k/N to all others, where k is the size of its cluster
(the group firing together); any oscillator pushed to x_j + Delta >= 1 merges into the firing cluster (absorption).
Initial voltages are i.i.d. Uniform(0, 1). The paper's exact recurrence (Sec. 3: next firer f = argmax (b_i, i),
"it absorbs every j allowed by the rule above, and each survivor moves to b_j' = b_j - b_f + N + s_f") fixes the
two details the prose leaves open, and this module follows it: (a) the pulse in one event is the firer's mass s_f
only (absorbed units join the cluster and reset with it but do not add a second pulse in the same event: no
avalanche); (b) absorbed units take the firer's phase (reset to 0 together with it).

Headline (plan row N01): gamma = 0, kick k/N, absorption, complete graph of NN = 200 nodes (area4: 800).

Census mapping and choices:
  * One census step = one free charging period (plan exception for N01): the time an uncoupled oscillator needs to
    charge from 0 to 1, i.e. 1 time unit at gamma = 0. Inside it the dynamics is integrated exactly, event by event
    (every firing and absorption), in float64 on cluster voltages. Clusters hold identical voltages; ties between
    distinct clusters have probability zero.
  * Output (graph world, census conventions): "Nph" = 2*pi*voltage as float32, frames every 5 census steps
    (phase convention of census.sim, phases in [0, 2pi)); "nt" int8 zeros (one node type, k = 1), frames every 2
    steps; "adj0" the complete graph without self-loops. The voltage at a frame time is the state after all events
    at or before that time (frames start at the random initial voltages, census step 0).
  * Once the state is persistent (a single cluster, or K equal clusters that complete a full round of K firings
    with no absorption) it is exactly periodic with period 1/K, so voltages at integer times no longer change; the
    remaining frames are copies of the first frame after that point. Unequal clusters always end in an absorption
    (Sec. 3), so the event loop terminates.
  * Knock-out (plan): kick = 0. Rule 1 alone gives x_i(t) = frac(x_i(0) + t), integrated in closed form (the same
    event dynamics with zero pulse; no unit ever reads another).

Fidelity check (plan: probability of full synchrony at N = 4, 5, 6; paper 0.887, 0.99968, 0.906 = eq. 4,
A_{N,1}/N^N = 227/256, 3124/3125, 42275/46656). Two parts, both run through this module's own event dynamics:
  (1) exact enumeration as in the paper (Sec. 3): each of the N^N words a in {0..N-1}^N is represented by the start
      x_i = (a_i + (i + 1/2)/N)/N (oscillator i has the i-th smallest fractional part); the simulator's outcome is
      constant on each word's cell, so the number of words ending in one cluster must equal A_{N,1} exactly. It is
      counted twice: with a vectorised batch version of the event rule, and start by start through _evolve (the
      code path simulate() uses);
  (2) Monte Carlo from i.i.d. uniform starts (n_mc starts per N, default 10,000,000).
  Tolerance (fixed before running): pass if both enumeration counts equal the paper's A_{N,1} for N = 4, 5, 6 and
  each Monte Carlo share lies within 4 binomial standard errors (computed from the paper's A_{N,1}/N^N) of it.

presence(): the plan records N01's new behaviour as a property of many starts across N, not observable in one census
run. presence() reports, per census run, whether the last frame is fully synchronised, the final cluster sizes, and
whether all clusters are equal in size (the only persistent outcomes at gamma = 0).
"""
import itertools
import time

import numpy as np

from . import base

CID = "planted:N01"
LAYERS = "N"
K = {"N": 1}
SPEED = 0.0
RULES = [("CHARGE_AND_FIRE", False),        # 1: linear charging to threshold 1, fire and reset (one-way)
         ("PULSE_KICK_ABSORB", True)]       # 2: firing group of k kicks all others by k/N; pushed to 1 => joins (I)
HEADLINE = {
    "gamma": 0.0,          # linear charging, Sec. 3 (the critical curve)
    "threshold": 1.0,      # Sec. 2
    "reset": 0.0,          # Sec. 2
    "kick_coeff": 1.0,     # pulse Delta = kick_coeff * k / N, Sec. 2 (knock-out: 0)
    "NN": 200,             # census graph size (plan: graphs 200 nodes; area4 800)
}
STEP = ("one free charging period (time 1 at dx/dt = 1, threshold 1); exact event-driven integration of every "
        "firing and absorption event inside it")
RULE_TEXT = ("World: N = 200 identical oscillators, each coupled to every other (no space). Each has a voltage between 0 "
             "and 1, drawn uniformly at random at the start. Rules: (1) all voltages rise at the same constant speed "
             "(straight-line charging, no leak); (2) when a voltage reaches 1 the oscillator fires and resets to 0, "
             "together with any oscillators already firing with it as a group; (3) a firing group of k oscillators "
             "kicks every other voltage up by k/N; (4) any oscillator pushed to 1 or above by a kick joins the firing "
             "group and fires with it from then on. Run for 1000 charging periods from one random start and record "
             "whether all 200 end as one group.")

TWO_PI = 2.0 * np.pi
PAPER_A = {4: (227, 256), 5: (3124, 3125), 6: (42275, 46656)}       # eq. (4)
PAPER_ROUNDED = {4: 0.887, 5: 0.99968, 6: 0.906}                     # card / plan values


# ----------------------------------------------------------------------------------------------- single realisation
def _evolve(x0, kick, times, until_persistent=False, max_events=1_000_000):
    """exact event-driven run of one population; returns voltages (len(times), N) at the requested sorted times.
    until_persistent=True (with times empty) runs only until the state is persistent and returns its final sizes.

    State: cluster voltages v (all < 1), cluster masses m, oscillator -> cluster labels lab."""
    N = len(x0)
    out = np.empty((len(times), N))
    if kick == 0.0:                                   # rule 1 alone: closed form
        for q, T in enumerate(times):
            out[q] = np.mod(x0 + T, 1.0)
        return out, {"events": None, "persistent_at": 0.0, "sizes": [1] * N}
    v = np.asarray(x0, float).copy(); m = np.ones(N); lab = np.arange(N)
    t = 0.0; q = 0; streak = 0; events = 0; persistent_at = None
    nT = len(times)

    def more():
        return q < nT or (until_persistent and persistent_at is None)
    while more() and events < max_events:
        f = int(np.argmax(v)); dt = 1.0 - v[f]
        while q < nT and times[q] < t + dt:          # frames strictly before the next firing
            out[q] = (v + (times[q] - t))[lab]; q += 1
            if persistent_at is not None and q < nT:
                out[q:] = out[q - 1]; q = nT          # exact periodicity: integer-time frames repeat
        if not more():
            break
        t += dt; events += 1
        vn = v + dt + kick * m[f] / N
        absorbed = vn >= 1.0
        absorbed[f] = False
        vn[f] = 0.0
        if absorbed.any():
            vn[absorbed] = 0.0
            m[f] += m[absorbed].sum()
            keep = ~absorbed
            newidx = np.cumsum(keep) - 1
            newidx[absorbed] = newidx[f]
            lab = newidx[lab]; v = vn[keep]; m = m[keep]
            streak = 0
        else:
            v = vn; streak += 1
        if persistent_at is None:
            Kc = len(m)
            if Kc == 1 or (streak >= Kc and np.all(m == m[0])):
                persistent_at = t
                if until_persistent and q >= nT:
                    break
    return out, {"events": events, "persistent_at": persistent_at, "sizes": sorted(m.astype(int).tolist(), reverse=True)}


def simulate(seed0, steps=1000, S=3, size="census", keep=None):
    on = base.kept(keep, RULES)
    kick = HEADLINE["kick_coeff"] if 2 in on else 0.0
    z = base.SIZES[size]; NN, W = z["NN"], z["W"]
    rec_ph, rec_nt = base.REC["Nph"], base.REC["nt"]
    times = np.arange(0, steps + 1, rec_ph, dtype=float)
    Nph = np.empty((S, base.n_frames(steps, "Nph"), NN), np.float32)
    for s, g in enumerate(base.rngs(seed0, S)):
        x0 = g.random(NN)                                  # i.i.d. Uniform(0, 1), Sec. 2
        x, _ = _evolve(x0, kick, times)
        Nph[s] = (TWO_PI * np.mod(x, 1.0)).astype(np.float32)
    nt = np.zeros((S, base.n_frames(steps, "nt"), NN), np.int8)
    adj0 = np.broadcast_to(~np.eye(NN, dtype=bool), (S, NN, NN)).copy()
    return {"Nph": Nph, "nt": nt, "adj0": adj0, "meta": base.meta(steps, W)}


# ------------------------------------------------------------------------------------------ batch outcomes (small N)
def _batch_outcome(x0, kick=1.0, max_events=100000):
    """final state of many independent small populations (rows of x0), run with the same event rule as _evolve but on
    oscillator voltages (members of a cluster carry identical float values). Returns (K_final, size_final, unresolved)
    where size_final = common cluster size for persistent equal-cluster states."""
    x = np.asarray(x0, float).copy(); M, N = x.shape
    Kf = np.zeros(M, int); act = np.arange(M); streak = np.zeros(M, int)
    for _ in range(max_events):
        if act.size == 0: break
        X = x[act]
        vmax = X.max(1, keepdims=True)
        fire = X == vmax
        k = fire.sum(1, keepdims=True)
        Xn = X + (1.0 - vmax) + kick * k / N
        absorbed = (~fire) & (Xn >= 1.0)
        Xn[fire | absorbed] = 0.0
        st = np.where(absorbed.any(1), 0, streak[act] + 1)
        cnt = (Xn[:, :, None] == Xn[:, None, :]).sum(2)
        Kc = np.rint((1.0 / cnt).sum(1)).astype(int)
        equal = cnt.min(1) == cnt.max(1)
        done = (Kc == 1) | (equal & (st >= Kc))
        x[act] = Xn; streak[act] = st
        Kf[act[done]] = Kc[done]
        act = act[~done]
    return Kf, N // np.maximum(Kf, 1), int(act.size)


def _enumerate(N):
    """number of words (out of N^N) whose representative start ends fully synchronised (paper eq. 3-4)."""
    words = np.array(list(itertools.product(range(N), repeat=N)), float)
    x0 = (words + (np.arange(N) + 0.5) / N) / N
    nsync = 0; unres = 0; Kcount = {}
    for c in range(0, len(x0), 200000):
        Kf, _, u = _batch_outcome(x0[c:c + 200000]); unres += u
        nsync += int((Kf == 1).sum())
        for kk, cc in zip(*np.unique(Kf, return_counts=True)): Kcount[int(kk)] = Kcount.get(int(kk), 0) + int(cc)
    return nsync, N ** N, unres, Kcount


def _enumerate_evolve(N):
    """the same enumeration run start by start through _evolve, the code path simulate() uses."""
    nsync = 0; unres = 0; Kcount = {}
    none = np.empty(0)
    for a in itertools.product(range(N), repeat=N):
        x0 = (np.array(a, float) + (np.arange(N) + 0.5) / N) / N
        _, inf = _evolve(x0, 1.0, none, until_persistent=True, max_events=100000)
        if inf["persistent_at"] is None: unres += 1; continue
        kk = len(inf["sizes"]); Kcount[kk] = Kcount.get(kk, 0) + 1
        nsync += int(kk == 1)
    return nsync, N ** N, unres, Kcount


def _monte_carlo(N, n, rng, chunk=250000):
    nsync = 0; unres = 0; Kcount = {}
    for c in range(0, n, chunk):
        mm = min(chunk, n - c)
        Kf, _, u = _batch_outcome(rng.random((mm, N))); unres += u
        nsync += int((Kf == 1).sum())
        for kk, cc in zip(*np.unique(Kf, return_counts=True)): Kcount[int(kk)] = Kcount.get(int(kk), 0) + int(cc)
    return nsync, unres, Kcount


def fidelity(budget_s=7200, n_mc=10_000_000, Ns=(4, 5, 6)):
    """plan check: P(full synchrony) at N = 4, 5, 6 by exact enumeration and Monte Carlo through the event dynamics."""
    t0 = time.time()
    rng = np.random.default_rng([base.seed_of(CID), 1668])
    res = {"check": "probability of full synchrony from uniform random starts, N = 4, 5, 6 (paper eq. 4)",
           "setting": {"gamma": 0.0, "kick": "k/N", "absorption": True},
           "tolerance": "enumeration counts (batch and simulate()'s _evolve) == paper A_{N,1}; "
                        "|P_MC - A_{N,1}/N^N| <= 4 binomial SE computed from the paper's A_{N,1}/N^N",
           "per_N": {}}
    # cost estimate from a small pilot, scale n_mc down if it would not fit the budget
    tp = time.time(); _monte_carlo(6, 20000, rng); per = (time.time() - tp) / 20000
    n_fit = int(0.8 * max(budget_s - 60, 1) / max(per * len(Ns), 1e-9))
    n = int(min(n_mc, n_fit)); res["n_mc"] = n
    ok = True
    for N in Ns:
        A, NN_, unres_e, Ke = _enumerate(N)                  # batch event code
        A2, _, unres_e2, Ke2 = _enumerate_evolve(N)          # simulate()'s event code, start by start
        s, unres_m, Km = _monte_carlo(N, n, rng)
        pa, pb = PAPER_A[N]
        p_paper = pa / pb; p_mc = s / n
        se = np.sqrt(p_paper * (1 - p_paper) / n)
        e_ok = (A == pa and NN_ == pb and unres_e == 0)
        e2_ok = (A2 == pa and unres_e2 == 0)
        m_ok = abs(p_mc - p_paper) <= 4 * se and unres_m == 0
        res["per_N"][N] = {"enumeration_batch": f"{A}/{NN_}", "enumeration_evolve": f"{A2}/{NN_}",
                           "paper": f"{pa}/{pb}", "paper_rounded": PAPER_ROUNDED[N],
                           "P_mc": p_mc, "se_mc_from_paper": se, "z": (p_mc - p_paper) / se,
                           "final_cluster_counts_enum": Ke, "final_cluster_counts_enum_evolve": Ke2,
                           "final_cluster_counts_mc": Km, "unresolved": unres_e + unres_e2 + unres_m,
                           "enum_pass": bool(e_ok), "enum_evolve_pass": bool(e2_ok), "mc_pass": bool(m_ok)}
        ok &= e_ok and e2_ok and m_ok
    res["pass"] = bool(ok); res["secs"] = round(time.time() - t0, 1)
    return res


# -------------------------------------------------------------------------------------------------------- presence
def presence(out):
    """final grouping of each census run (from the Nph frames only): fully synchronised or not, cluster sizes."""
    ph = np.asarray(out["Nph"])
    S, F, NN = ph.shape
    runs = []
    for s in range(S):
        last = ph[s, -1]
        _, cnt = np.unique(np.round(last.astype(np.float64), 5), return_counts=True)
        sizes = sorted(cnt.tolist(), reverse=True)
        settled = None
        for f in range(F - 1, -1, -1):              # earliest frame from which the frames stop changing
            if not np.array_equal(ph[s, f], last): break
            settled = f
        runs.append({"fully_synchronised": len(sizes) == 1, "n_clusters": len(sizes), "cluster_sizes": sizes[:20],
                     "equal_sizes": len(set(sizes)) == 1,
                     "settled_by_step": None if settled is None else int(settled * out["meta"]["rec"]["Nph"])})
    return {"measure": "final grouping per run (synchrony / equal-size clusters)",
            "note": "plan: the new behaviour (P_sync set by the prime factors of N) is a property of many starts "
                    "across N and is not observable in one census run",
            "share_synchronised": float(np.mean([r["fully_synchronised"] for r in runs])),
            "runs": runs}
