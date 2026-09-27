"""Census filter (pilot draft v0): simulator output -> per-layer observation views -> screen -> known filter.

Principle: a view carries only what an outside observer could measure from the raw state (grid, positions, headings,
phases, adjacency, fields) plus geometry (box size, grid size, frame spacing). It never carries mechanism flags or
rule parameters — that would hand the detector the answer.

  screen       emergent  := generic_emergence score >= 0.5  OR  model-free complexity (surrogate + dead-state gated)
  known filter MATCH(P)  := the 36-detector battery, match_min_tier = "confirmation" (the T2c OOD operating point)
  per view     MATCH(P) if P matches in >= 2 of 3 seeds; UNCLASSIFIED if emergent in >= 2 seeds and no such match;
               otherwise NONE.  A program is flagged UNCLASSIFIED if any of its views is.
"""
import numpy as np
from .sim import W, NA, NN

_BATTERY = _FNS = None


def _load():
    global _BATTERY, _FNS
    if _FNS is None:
        from analysis.battery_profile import build_detector_fns
        from epc.phase2a.battery import Battery
        _BATTERY, _FNS = Battery.load(), build_detector_fns()
    return _BATTERY, _FNS


def _phase_frames(ph, rec):
    out = []
    for f in range(ph.shape[0]):
        th = np.mod(ph[f], 2 * np.pi)
        out.append({"theta": th.ravel(), "r": float(np.abs(np.exp(1j * th).mean())), "phases": th, "step": f * rec})
    return out


def views(out, prog, s):
    """[(name, history, metadata)] for seed s."""
    rec = out["meta"]["rec"]; V = []
    if "C" in out:
        C = out["C"][s].astype(np.int64); k = max(prog.k["C"], 2); hist = []
        for f in range(C.shape[0]):
            act = float((C[f] != C[f - 1]).mean()) if f else 0.0
            hist.append({"grid": C[f], "grid_dims": (W, W), "n_states": k, "step": f * rec, "activity_density": act})
        V.append(("C", hist, {"rows": W, "cols": W, "boundary": "periodic", "neighborhood": "moore",
                              "substrate_type": "lattice_2d", "n_states": k}))
    if "aval" in out:
        a = out["aval"][s].astype(float)
        V.append(("C.aval", [{"avalanche_sizes": a[a > 0], "activity": a, "step": 0}], {}))
    if "Cph" in out:
        V.append(("C.phase", _phase_frames(out["Cph"][s], rec), {"N": W * W, "substrate_type": "oscillator"}))
    if "pos" in out:
        v = prog.agent[0]; hist = []
        for f in range(out["pos"].shape[1]):
            h = np.mod(out["head"][s, f], 2 * np.pi)
            fr = {"positions": out["pos"][s, f], "velocities": v * np.stack([np.cos(h), np.sin(h)], -1),
                  "headings": h, "box_size": float(W), "step": f * rec}
            if prog.k["A"] >= 2: fr["labels"] = out["at"][s, f].astype(np.int64)
            hist.append(fr)
        V.append(("A", hist, {"n_particles": NA, "box_size": float(W), "dt": 1.0,
                              "space_type": "continuous_2d_periodic"}))
    if "nt" in out:
        k = prog.k["N"]; hist = []
        adj = out["adj"][s] if "adj" in out else None
        for f in range(out["nt"].shape[1]):
            A = (adj[f] if adj is not None else out["adj0"][s]).astype(np.int64)
            op = out["nt"][s, f].astype(float) / max(k - 1, 1)
            hist.append({"adjacency": A, "opinions": op, "step": f * rec})
        V.append(("N", hist, {"N": NN, "substrate_type": "network"}))
    if "Nph" in out:
        V.append(("N.phase", _phase_frames(out["Nph"][s], rec), {"N": NN, "substrate_type": "oscillator"}))
    if "F" in out:
        for fi in range(out["F"].shape[2]):
            hist = [{"field": out["F"][s, f, fi], "step": f * rec} for f in range(out["F"].shape[1])]
            V.append((f"F{fi}", hist, {"rows": W, "cols": W, "boundary": "periodic", "substrate_type": "field_2d"}))
    return V


def screen(history):
    from epc.phase2a.emergence import generic_emergence
    from epc.phase2a.novelty_tripwire import model_free_complexity
    try:
        em = generic_emergence(history, seed=0)
    except Exception as e:
        em = {"score": 0.0, "kind": f"error:{type(e).__name__}"}
    try:
        mf = model_free_complexity(history)
    except Exception as e:
        mf = {"is_complex": False, "C": None, "psi": None, "struct": None, "collapsed": None}
    sc = float(em.get("score", 0.0) or 0.0)
    return {"em_score": round(sc, 4), "em_kind": em.get("kind"), "is_complex": bool(mf.get("is_complex")),
            "C": mf.get("C"), "psi": mf.get("psi"), "collapsed": mf.get("collapsed"),
            "emergent": bool(sc >= 0.5 or mf.get("is_complex"))}


def known(history, metadata):
    from analysis.battery_profile import profile_observation
    B, F = _load()
    prof = profile_observation(history, metadata, battery=B, detector_fns=F, match_min_tier="confirmation")
    fired = [(r["pattern_id"], r["detector_tier"]) for r in prof["profile"] if r.get("fired_detected")]
    v = prof["verdict"]
    return {"verdict": v["verdict"], "pattern": v.get("pattern_id"), "tier": v.get("detector_tier"),
            "demoted": v.get("demoted_match"), "fired": fired}


def evaluate(out, prog, battery_all=False):
    """per view: per-seed screen (+ known filter on emergent seeds), then the view verdict."""
    S = out["C"].shape[0] if "C" in out else next(v.shape[0] for k, v in out.items() if k != "meta" and hasattr(v, "shape"))
    per = {}
    for s in range(S):
        for name, hist, md in views(out, prog, s):
            r = screen(hist) if name != "C.aval" else {"emergent": True, "em_score": None, "em_kind": "avalanche-bundle"}
            if r["emergent"] or battery_all: r.update(known(hist, md))
            per.setdefault(name, []).append(r)
    res = {}
    for name, rs in per.items():
        em = sum(r["emergent"] for r in rs)
        pats = [r.get("pattern") for r in rs if r.get("verdict") == "MATCH"]
        top = max(set(pats), key=pats.count) if pats else None
        if top and pats.count(top) >= 2: v = "MATCH"
        elif em >= 2: v = "UNCLASSIFIED"
        else: v = "NONE"
        res[name] = {"verdict": v, "pattern": top if v == "MATCH" else None, "emergent_seeds": em, "seeds": rs}
    return res
