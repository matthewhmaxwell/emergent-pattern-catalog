"""Interaction knock-out test (round 6; the EPC instrument's collapse-ablation, applied to every flag).

A flagged pattern counts as EMERGENT only if it depends on interaction: the same program re-run with its interaction
rules removed (same seeds) must either lose the pattern or look clearly different. Otherwise the "pattern" came from
the initial condition, one-way conversion, or independent parts — not from interaction.

Interaction rules = rules that read other entities (neighbour copying / counting, alignment, attraction, repulsion,
quorum, toppling, games, rewiring, coupling, agent<->cell sensing, gradient climbing, autocatalysis between fields).
Kept in the knock-out: one-way rules (spontaneous flips, decay, feed, emission, deposition). Phase coupling is set to
K = 0 rather than removed, so the phase view still exists to compare.

Decision per view and seed (fixed before round 6):
  interaction-driven  <=>  knock-out view absent  OR  knock-out view not emergent  OR  d > D_MIN
  d = median over fingerprint features of |real - ko| / (|real| + |ko| + 0.05);  D_MIN = 0.25
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
    """real, ko: outputs of per_seed_view_results. -> {view: [bool per seed]} (only meaningful where real is emergent)."""
    out = {}
    for v, seeds in real.items():
        kv = ko.get(v); flags = []
        for s, rs in enumerate(seeds):
            if kv is None or s >= len(kv) or not kv[s]["emergent"]: flags.append(True)
            else: flags.append(fp_distance(rs["fp"], kv[s]["fp"]) > D_MIN)
        out[v] = flags
    return out


def run_with_knockout(p, seed0, **kw):
    """simulate p and its knock-out (same seeds) -> (raw output, real per-seed results, driven flags per view/seed)."""
    from census import sim
    out = sim.run(p, seed0=seed0, **kw); real = per_seed_view_results(out, p)
    if not any(x["emergent"] for seeds in real.values() for x in seeds):
        return out, real, {v: [False] * len(s) for v, s in real.items()}
    kp = knockout_program(p)
    ko = per_seed_view_results(sim.run(kp, seed0=seed0, **kw), kp)
    return out, real, interaction_driven(real, ko)
