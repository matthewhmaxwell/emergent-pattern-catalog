"""Interaction knock-out test (round 6; the EPC instrument's collapse-ablation, applied to every flag).

A flagged pattern counts as EMERGENT only if it depends on interaction: the same program re-run with its interaction
rules removed (same seeds) must either lose the pattern or look clearly different. Otherwise the "pattern" came from
the initial condition, one-way conversion, or independent parts — not from interaction.

Interaction rules = rules that read other entities (neighbour copying / counting, alignment, attraction, repulsion,
quorum, toppling, games, rewiring, coupling, agent<->cell sensing, gradient climbing, autocatalysis between fields).
Kept in the knock-out: one-way rules (spontaneous flips, decay, feed, emission, deposition). Phase coupling is set to
K = 0 rather than removed, so the phase view still exists to compare.

Decision per view (round 10c; noise-aware, feature-wise):
  interaction-driven  <=>  knock-out view absent  OR  some fingerprint feature f changes clearly:
      |mean_real(f) - mean_ko(f)| >= Z_MIN * noise(f)   AND   relative change >= REL_MIN
  over the seeds where the real run is emergent; noise(f) = max(sd_real, sd_ko, 5% of the feature's magnitude, 1e-3).
  Z_MIN = 4, REL_MIN = 0.25. Round-10 lessons: the knock-out's own screen flag is erratic on pure noise (not used);
  a median over features is dominated by features the interaction does not touch (speed, graph shape), so the test
  looks for ANY feature the interaction clearly changes.
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


Z_MIN, REL_MIN = 4.0, 0.25


def interaction_driven(real, ko):
    """real, ko: outputs of per_seed_view_results -> {view: [bool per seed]} (one decision per view)."""
    out = {}
    for v, seeds in real.items():
        kv = ko.get(v)
        if kv is None:
            out[v] = [True] * len(seeds); continue
        em = [s for s in range(len(seeds)) if seeds[s]["emergent"] and s < len(kv)]
        driven = False
        if em:
            keys = sorted(set().union(*[seeds[s]["fp"] for s in em]))
            for f in keys:
                a = np.array([seeds[s]["fp"].get(f, 0.0) for s in em]); b = np.array([kv[s]["fp"].get(f, 0.0) for s in em])
                diff = abs(a.mean() - b.mean()); mag = abs(a.mean()) + abs(b.mean())
                noise = max(a.std(), b.std(), 0.05 * mag, 1e-3)
                if diff >= Z_MIN * noise and diff / (mag + 0.05) >= REL_MIN:
                    driven = True; break
        out[v] = [driven] * len(seeds)
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
