"""Interaction knock-out test (round 6; the EPC instrument's collapse-ablation, applied to every flag).

A flagged pattern counts as EMERGENT only if it depends on interaction: the same program re-run with its interaction
rules removed (same seeds) must either lose the pattern or look clearly different. Otherwise the "pattern" came from
the initial condition, one-way conversion, or independent parts — not from interaction.

Interaction rules = rules that read other entities (neighbour copying / counting, alignment, attraction, repulsion,
quorum, toppling, games, rewiring, coupling, agent<->cell sensing, gradient climbing, autocatalysis between fields).
Kept in the knock-out: one-way rules (spontaneous flips, decay, feed, emission, deposition). Phase coupling is set to
K = 0 rather than removed, so the phase view still exists to compare.

Decision per view (round 10; noise-aware — round 9 extended negatives showed a noisy program whose knock-out
differed from it by no more than two seeds of the same program differ from each other):
  interaction-driven  <=>  knock-out view absent  OR  d_between > max(D_MIN, 2 * d_within)
  (the knock-out's own screen flag is NOT used: the generic screen fires erratically on pure noise, so "the
  knock-out was not flagged" is not evidence of interaction — round-10 lesson)
  d(a, b)   = median over fingerprint features of |a - b| / (|a| + |b| + 0.05);  D_MIN = 0.25
  d_between = median over seeds of d(real_s, ko_s)          (seeds where the real run is emergent)
  d_within  = median over seed pairs of d(real_s, real_t)   (seed-to-seed variation of the real program)
The same decision applies to every seed of the view.
"""
import numpy as np
from census import grammar_g1 as G

ONE_WAY = {"C.FLIP", "N.FLIP", "F.DECAY", "F.FEED", "C.EMIT", "A.DEPOSIT"}
D_MIN = 0.25


def knockout_program(p):
    rules = []
    for r in p.rules:
        if r.tmpl in ONE_WAY: rules.append(r)
        elif r.tmpl in ("C.KURA", "N.KURA"): rules.append(G.Rule(r.tmpl, (("K", 0.0),)))
    return G.make(p.layers, p.k, p.nf, p.D, p.agent, rules)


def is_interaction_free(p):
    return all(r.tmpl in ONE_WAY for r in p.rules)


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
            sc = FL.screen(th) if name != "C.aval" else FL.screen_avalanche(hist)
            fp = fingerprint(name, th, adj0=out["adj0"][s].astype(np.int64) if "adj0" in out else None) if name != "C.aval" else {}
            res.setdefault(name, []).append({"emergent": bool(sc["emergent"]), "fp": fp, "em_score": sc.get("em_score")})
    return res


def interaction_driven(real, ko):
    """real, ko: outputs of per_seed_view_results. -> {view: [bool per seed]} (one decision per view, noise-aware)."""
    out = {}
    for v, seeds in real.items():
        kv = ko.get(v)
        if kv is None:
            out[v] = [True] * len(seeds); continue
        em = [s for s in range(len(seeds)) if seeds[s]["emergent"]]
        both = [s for s in em if s < len(kv)]
        d_between = float(np.median([fp_distance(seeds[s]["fp"], kv[s]["fp"]) for s in both])) if both else 1.0
        pairs = [(a, b) for i, a in enumerate(em) for b in em[i + 1:]]
        d_within = float(np.median([fp_distance(seeds[a]["fp"], seeds[b]["fp"]) for a, b in pairs])) if pairs else 0.0
        out[v] = [d_between > max(D_MIN, 2 * d_within)] * len(seeds)
    return out


def run_with_knockout(p, seed0, **kw):
    """simulate p and its knock-out (same seeds) -> (raw output, real per-seed results, driven flags per view/seed,
    raw knock-out output or None when nothing was emergent)."""
    from census import sim
    out = sim.run(p, seed0=seed0, **kw); real = per_seed_view_results(out, p)
    if not any(x["emergent"] for seeds in real.values() for x in seeds):
        return out, real, {v: [False] * len(s) for v, s in real.items()}, None
    kp = knockout_program(p); ko_out = sim.run(kp, seed0=seed0, **kw)
    ko = per_seed_view_results(ko_out, kp)
    return out, real, interaction_driven(real, ko), ko_out
