"""Interaction knock-out test (round 6; the EPC instrument's collapse-ablation, applied to every flag).

A flagged pattern counts as EMERGENT only if it depends on interaction: the same program re-run with its interaction
rules removed (same seeds) must either lose the pattern or look clearly different. Otherwise the "pattern" came from
the initial condition, one-way conversion, or independent parts — not from interaction.

Interaction rules = rules that read other entities (neighbour copying / counting, alignment, attraction, repulsion,
quorum, toppling, games, rewiring, coupling, agent<->cell sensing, gradient climbing, autocatalysis between fields).
Kept in the knock-out: one-way rules (spontaneous flips, decay, feed, emission, deposition). Phase coupling is set to
K = 0 rather than removed, so the phase view still exists to compare.

Decision per view (final form; "does the detected emergence signal disappear?"), PAIRED over ALL seeds:
  evidence e = max(generic screen score, min(1, consensus gain / 0.4)) per seed — continuous channels only; the
  yes/no model-free complexity flag is NOT evidence (placebo test, 2026-10-02: it is a coin flip on several kinds
  of interaction-free data, and one random yes among 3 seeds fakes an evidence drop).
  d_s = e_real(s) - e_knockout(s) for every seed s (same seeds in both runs).
  interaction-driven  <=>  knock-out view absent
      OR  mean(d) >= DROP_MIN  AND  mean(d) >= Z_MIN * max(sd(d) / sqrt(S), 0.02)                    [route A]
      OR  for some textbook ORDER measure f, with g_s = direction * (f_real(s) - f_knockout(s)):
          mean(g) >= ORDER_MIN  AND  mean(g) >= 4 * max(sd(g) / sqrt(S), 0.01)  AND  relative change >= 0.25  [route B]
  DROP_MIN = 0.25, Z_MIN = 2, ORDER_MIN = 0.1.
  Order measures only (organisation, not activity or side effects such as links moved): lattice Moran's I, largest
  domain, correlation length; agent polarization, local alignment, clustering (lower nearest-neighbour ratio),
  type segregation; network type modularity, consensus gain; phase order r and local r; field Moran's I, peak
  sharpness.
  Lessons that shaped it: the knock-out's own screen FLAG must not be used (round 10); "any feature changed" lets
  trivial side effects through (round 10c); uncoupled oscillators need the order route (round 11); comparing only
  the seeds where the REAL run is emergent is a selection bias (final gate: a coin-flip screen on network noise made
  a null look interaction-driven) -> the test is paired over all seeds. With 3 seeds the paired test cannot protect
  against a screen channel that fires at random (second final gate + placebo: model-free complexity on oscillator
  phases) -> every channel must be quiet on interaction-free data, which census.placebo now measures directly.
"""
import numpy as np
from census import grammar_g1 as G

ONE_WAY = {"C.FLIP", "N.FLIP", "F.DECAY", "F.FEED", "C.EMIT", "A.DEPOSIT"}
D_MIN = 0.25


G2_ONE_WAY_ACTS = {"SET", "SET_NEXT", "EMIT", "DEPOSIT", "SCALE", "RELAX1", "TURN_L", "TURN_R", "SLOW"}


def _one_way(r):
    """a rule that reads no other entity: G1 one-way templates; G2 = unconditional rule with a non-reading action
    (mirrors G1: flips, decay, feed, emission, deposition; plus G2's unconditional cycle / turn / slow)."""
    if r.tmpl != "G2": return r.tmpl in ONE_WAY
    return r.cond == "ALWAYS" and r.act in G2_ONE_WAY_ACTS


def knockout_program(p):
    rules = []
    for r in p.rules:
        if _one_way(r): rules.append(r)
        elif r.tmpl in ("C.KURA", "N.KURA"): rules.append(G.Rule(r.tmpl, (("K", 0.0),)))
        elif r.tmpl == "G2" and r.act == "PHASE_COUPLE":                    # keep the oscillators, cut the coupling
            rules.append(type(r)(r.target, "any", "ALWAYS", (), "PHASE_COUPLE", (("K", 0.0),), 1.0))
    return G.make(p.layers, p.k, p.nf, p.D, p.agent, rules)


def is_interaction_free(p):
    return all(_one_way(r) for r in p.rules)


def fp_distance(a, b):
    keys = sorted(set(a) | set(b))
    if not keys: return 0.0
    return float(np.median([abs(a.get(k, 0.0) - b.get(k, 0.0)) / (abs(a.get(k, 0.0)) + abs(b.get(k, 0.0)) + 0.05) for k in keys]))


def per_seed_view_results(out, p):
    """{view: [ {emergent, fp} per seed ]} using the census screen + fingerprint (same thinning as the census)."""
    from census import filter as FL
    from census.fingerprint import fingerprint
    res = {}
    S = next(v.shape[0] for k, v in out.items() if k != "meta" and hasattr(v, "shape"))
    for s in range(S):
        for name, hist, md in FL.views(out, p, s):
            th = FL._screen_hist(name, hist, "adj" in out)
            if name == "C.aval": sc = FL.screen_avalanche(hist)
            elif name == "N":                                   # node states: network channels need the graph
                adj = out["adj0"][s] if "adj0" in out else hist[-1]["adjacency"]
                sc = FL.screen(th, adj=adj, network=True)
            else: sc = FL.screen(th, complexity=not name.endswith("phase"))   # phases: generic score only
            fp = fingerprint(name, th, adj0=out["adj0"][s].astype(np.int64) if "adj0" in out else None) if name != "C.aval" else {}
            # evidence for the knock-out comparison: CONTINUOUS channels only. The yes/no model-free complexity
            # flag still counts for "emergent" (non-phase views) but not here: it fires at in-between rates on
            # interaction-free data (placebo: 40% of uncoupled oscillators, 256 / 942 field nulls, static agents),
            # and one random yes/no among 3 seeds is enough to fake an evidence drop.
            ev = max(float(sc.get("em_score") or 0.0), min(1.0, float(sc.get("consensus_gain") or 0.0) / 0.4))
            res.setdefault(name, []).append({"emergent": bool(sc["emergent"]), "fp": fp, "em_score": sc.get("em_score"),
                                             "evidence": ev, "em_kind": sc.get("em_kind"),
                                             "is_complex": bool(sc.get("is_complex")),
                                             "consensus_gain": float(sc.get("consensus_gain") or 0.0)})
    return res


DROP_MIN, Z_MIN, ORDER_MIN = 0.25, 2.0, 0.1
ORDER = {"g_moran": 1, "g_largest_domain": 1, "g_corrlen": 1, "a_polar": 1, "a_local_align": 1, "a_nn_ratio": -1,
         "a_type_segregation": 1, "n_type_modularity": 1, "n_consensus_gain": 1, "p_r": 1, "p_local_r": 1,
         "f_moran": 1, "f_peak_sharp": 1}


def interaction_driven(real, ko):
    """real, ko: outputs of per_seed_view_results -> {view: [bool per seed]} (one decision per view, paired seeds)."""
    out = {}
    for v, seeds in real.items():
        kv = ko.get(v)
        if kv is None:
            out[v] = [True] * len(seeds); continue
        S = min(len(seeds), len(kv)); driven = False
        if S:
            d = np.array([seeds[s]["evidence"] - kv[s]["evidence"] for s in range(S)])
            driven = d.mean() >= DROP_MIN and d.mean() >= Z_MIN * max(d.std() / np.sqrt(S), 0.02)
            for f, sgn in ORDER.items():
                if driven: break
                if not all(f in seeds[s]["fp"] and f in kv[s]["fp"] for s in range(S)): continue
                x = np.array([seeds[s]["fp"][f] for s in range(S)]); y = np.array([kv[s]["fp"][f] for s in range(S)])
                g = sgn * (x - y)
                if (g.mean() >= ORDER_MIN and g.mean() >= 4 * max(g.std() / np.sqrt(S), 0.01)
                        and abs(x.mean() - y.mean()) / (abs(x.mean()) + abs(y.mean()) + 0.05) >= 0.25):
                    driven = True
        out[v] = [bool(driven)] * len(seeds)
    return out


def run_with_knockout(p, seed0, **kw):
    """simulate p and its knock-out (same seeds) -> (raw output, real per-seed results, driven flags per view/seed,
    raw knock-out output or None when nothing was emergent)."""
    from census import sim
    out = sim.run(p, seed0=seed0, **kw); real = per_seed_view_results(out, p)
    if not any(sum(1 for x in seeds if x["emergent"]) >= 2 for seeds in real.values()):
        # no view is emergent in >= 2 seeds -> nothing can be flagged whatever the knock-out says: skip it
        return out, real, {v: [False] * len(s) for v, s in real.items()}, None
    kp = knockout_program(p); ko_out = sim.run(kp, seed0=seed0, **kw)
    ko = per_seed_view_results(ko_out, kp)
    return out, real, interaction_driven(real, ko), ko_out
