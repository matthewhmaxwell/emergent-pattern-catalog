"""Interaction knock-out test (round 6; the EPC instrument's collapse-ablation, applied to every flag).

A flagged pattern counts as EMERGENT only if it depends on interaction: the same program re-run with its interaction
rules removed (same seeds) must either lose the pattern or look clearly different. Otherwise the "pattern" came from
the initial condition, one-way conversion, or independent parts — not from interaction.

Interaction rules = rules that read other entities (neighbour copying / counting, alignment, attraction, repulsion,
quorum, toppling, games, rewiring, coupling, agent<->cell sensing, gradient climbing, autocatalysis between fields).
Kept in the knock-out: one-way rules (spontaneous flips, decay, feed, emission, deposition). Phase coupling is set to
K = 0 rather than removed, so the phase view still exists to compare.

Decision per view (round 11; "does the detected emergence signal disappear?"):
  evidence e = max(generic emergence score, 1 if model-free complexity fired, min(1, consensus gain / 0.4)) — the
  same signals the census screen flags on, as one continuous number per seed.
  interaction-driven  <=>  knock-out view absent  OR
      mean_real(e) - mean_ko(e) >= max(DROP_MIN, Z_MIN * noise),  noise = max(sd_real(e), sd_ko(e), 0.02)
  over the seeds where the real run is emergent. DROP_MIN = 0.25, Z_MIN = 2.
  OR (route B, round 11b) some textbook ORDER measure is clearly higher with interaction than without:
      direction * (mean_real - mean_ko) >= max(ORDER_MIN, 4 * noise)  AND  relative change >= 0.25; ORDER_MIN = 0.1.
  Order measures only (organisation, not activity or side effects such as links moved): lattice Moran's I, largest
  domain, correlation length; agent polarization, local alignment, clustering (lower nearest-neighbour ratio),
  type segregation; network type modularity, consensus gain; phase order r and local r; field Moran's I, peak
  sharpness. Round-11 lesson: uncoupled oscillators rotate so regularly that the generic screen scores them high,
  so for sync only the order measure separates real from knock-out.
  Lessons: the knock-out's own screen FLAG is erratic on pure noise (round 10) -> compare seed-averaged evidence
  with its spread instead; "any fingerprint feature changed" lets trivial side effects through (round 10c: one early
  rewiring step changed the graph while the flag came from forced conversion, present in the knock-out too).
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
            ev = max(float(sc.get("em_score") or 0.0), 1.0 if sc.get("is_complex") else 0.0,
                     min(1.0, float(sc.get("consensus_gain") or 0.0) / 0.4))
            res.setdefault(name, []).append({"emergent": bool(sc["emergent"]), "fp": fp, "em_score": sc.get("em_score"),
                                             "evidence": ev})
    return res


DROP_MIN, Z_MIN, ORDER_MIN = 0.25, 2.0, 0.1
ORDER = {"g_moran": 1, "g_largest_domain": 1, "g_corrlen": 1, "a_polar": 1, "a_local_align": 1, "a_nn_ratio": -1,
         "a_type_segregation": 1, "n_type_modularity": 1, "n_consensus_gain": 1, "p_r": 1, "p_local_r": 1,
         "f_moran": 1, "f_peak_sharp": 1}


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
            a = np.array([seeds[s]["evidence"] for s in em]); b = np.array([kv[s]["evidence"] for s in em])
            noise = max(a.std(), b.std(), 0.02)
            driven = (a.mean() - b.mean()) >= max(DROP_MIN, Z_MIN * noise)
            for f, sgn in ORDER.items():
                if driven: break
                if not all(f in seeds[s]["fp"] and f in kv[s]["fp"] for s in em): continue
                x = np.array([seeds[s]["fp"][f] for s in em]); y = np.array([kv[s]["fp"][f] for s in em])
                gain = sgn * (x.mean() - y.mean()); nz = max(x.std(), y.std(), 0.01)
                if gain >= max(ORDER_MIN, 4 * nz) and abs(x.mean() - y.mean()) / (abs(x.mean()) + abs(y.mean()) + 0.05) >= 0.25:
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
