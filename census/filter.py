"""Census filter (pilot draft v0): simulator output -> per-layer observation views -> screen -> known filter.

Principle: a view carries only what an outside observer could measure from the raw state (grid, positions, headings,
phases, adjacency, fields) plus geometry (box size, grid size, frame spacing). It never carries mechanism flags or
rule parameters — that would hand the detector the answer.

  screen       emergent  := generic_emergence score >= 0.5  OR  model-free complexity (surrogate + dead-state gated)
  known filter MATCH(P)  := the 36-detector battery, match_min_tier = "confirmation" (the T2c OOD operating point),
               run in two stages for cost: (1) FAST — every detector with reduced null counts on <= FAST_FRAMES frames,
               only to see which detectors fire at all; (2) CONFIRM — each fired detector re-run with its validated
               settings on <= CONFIRM_FRAMES frames; MATCH needs a validated firing at >= confirmation. The fast stage can
               only cause misses (known -> "unclassified", which triage re-checks with the full battery), never a
               false MATCH.
  cost control the screen runs on a thinned copy of each view (<= SCREEN_FRAMES frames); the battery sees full resolution
  absorption   a lattice run that absorbs into a single type also gets a "transient" view (frames up to absorption),
               so coarsening-to-consensus is judged on its dynamics rather than on the dead final state
  per view     MATCH(P) if P matches in >= 2 of 3 seeds; UNCLASSIFIED if emergent in >= 2 seeds and no such match;
               otherwise NONE.  A program is flagged UNCLASSIFIED if any of its views is.
"""
import numpy as np
from .sim import W, NA, NN

_BATTERY = _FNS = None
SCREEN_FRAMES, FAST_FRAMES, CONFIRM_FRAMES = 101, 125, 251
_FAST = None
_TIER = {"definitive": 3, "confirmation": 2, "screening": 1, "none": 0}


def _load():
    global _BATTERY, _FNS
    if _FNS is None:
        from analysis.battery_profile import build_detector_fns
        from epc.phase2a.battery import Battery
        _BATTERY, _FNS = Battery.load(), build_detector_fns()
    return _BATTERY, _FNS


def _fast_fns():
    """validated factories with reduced null counts (19 permutations / null runs, 2 P15 variations)."""
    global _FAST
    if _FAST is None:
        import inspect
        import analysis.run_phase2a_panel as R
        from analysis.battery_profile import PATTERNS
        red = {"n_permutations": 19, "n_null_runs": 19, "n_boot": 100, "n_variations": 2}
        _FAST = {}
        for pid in PATTERNS:
            mk = getattr(R, f"make_{pid}_detector_fn", None)
            if mk is None: continue
            sig = inspect.signature(mk).parameters
            _FAST[pid.upper()] = mk(**{k: v for k, v in red.items() if k in sig})
    return _FAST


def _thin(hist, n):
    return hist[::max(1, -(-len(hist) // n))] if len(hist) > 1 else hist


def _phase_frames(ph, rec):
    out = []
    for f in range(ph.shape[0]):
        th = np.mod(ph[f], 2 * np.pi)
        out.append({"theta": th.ravel(), "r": float(np.abs(np.exp(1j * th).mean())), "phases": th, "step": f * rec})
    return out


def _grid_hist(C, k, rec):
    hist = []
    for f in range(C.shape[0]):
        act = float((C[f] != C[f - 1]).mean()) if f else 0.0
        hist.append({"grid": C[f], "grid_dims": (W, W), "n_states": k, "step": f * rec, "activity_density": act})
    return hist


def views(out, prog, s):
    """[(name, history, metadata)] for seed s."""
    rec = out["meta"]["rec"]; V = []
    if "C" in out:
        C = out["C"][s].astype(np.int64); k = max(prog.k["C"], 2)
        md = {"rows": W, "cols": W, "boundary": "periodic", "neighborhood": "moore", "substrate_type": "lattice_2d",
              "n_states": k}
        hist = _grid_hist(C, k, rec["C"])
        if any(r.tmpl == "C.GAME" or getattr(r, "act", None) == "IMITATE_BEST" and r.target == "C" for r in prog.rules):
            for h in hist: h["coop_fraction"] = float((h["grid"] == 0).mean())    # grammar fixes type 0 = cooperate
        V.append(("C", hist, md))
        ab = out["meta"].get("absorbed_at")
        if ab is not None and len(np.unique(C[-1])) == 1 and ab // rec["C"] >= 20:
            V.append(("C.transient", _grid_hist(C[:ab // rec["C"] + 1], k, rec["C"]), md))
    if "aval" in out:
        a = out["aval"][s].astype(float); d = out["aval_dur"][s].astype(float); m = a > 0
        V.append(("C.aval", [{"avalanche_sizes": a[m], "avalanche_durations": d[m], "activity": a, "step": 0}], {}))
    if "Cph" in out:
        V.append(("C.phase", _phase_frames(out["Cph"][s], rec["Cph"]), {"N": W * W, "substrate_type": "oscillator"}))
    if "pos" in out:
        v = prog.agent[0]; hist = []; B = float(out["meta"].get("WA", W))
        for f in range(out["pos"].shape[1]):
            h = np.mod(out["head"][s, f], 2 * np.pi)
            fr = {"positions": out["pos"][s, f], "velocities": v * np.stack([np.cos(h), np.sin(h)], -1),
                  "headings": h, "box_size": B, "step": f * rec["pos"]}
            if prog.k["A"] >= 2: fr["labels"] = out["at"][s, f].astype(np.int64)
            hist.append(fr)
        V.append(("A", hist, {"n_particles": NA, "box_size": B, "dt": 1.0,
                              "space_type": "continuous_2d_periodic"}))
    if "nt" in out:
        k = prog.k["N"]; hist = []
        adj = out["adj"][s] if "adj" in out else None; step = rec["adj"] if adj is not None else rec["nt"]
        ratio = step // rec["nt"]
        for f in range(adj.shape[0] if adj is not None else out["nt"].shape[1]):
            A = (adj[f] if adj is not None else out["adj0"][s]).astype(np.int64)
            op = out["nt"][s, f * ratio].astype(float) / max(k - 1, 1)
            hist.append({"adjacency": A, "opinions": op, "step": f * step})
        V.append(("N", hist, {"N": NN, "substrate_type": "network"}))
    if "Nph" in out:
        V.append(("N.phase", _phase_frames(out["Nph"][s], rec["Nph"]), {"N": NN, "substrate_type": "oscillator"}))
    if "F" in out:
        for fi in range(out["F"].shape[2]):
            hist = [{"field": out["F"][s, f, fi].astype(np.float64), "step": f * rec["F"]} for f in range(out["F"].shape[1])]
            V.append((f"F{fi}", hist, {"rows": W, "cols": W, "boundary": "periodic", "substrate_type": "field_2d"}))
    return V


def _screen_hist(name, hist, rewiring):
    """cost-bounded copy of a view for the screen + fingerprint (the battery stage always gets full resolution)."""
    if name == "N":
        if not rewiring:                                     # static graph: the dynamics live in the node states
            return [{"opinions": h["opinions"], "step": h["step"]} for h in _thin(hist, SCREEN_FRAMES)]
        return _thin(hist, 25)
    if name == "A": return _thin(hist, 51)
    return _thin(hist, SCREEN_FRAMES)


def screen_avalanche(history):
    """avalanche bundles have no generic lens: emergent iff the SOC detector (P14) fires at >= screening."""
    from epc.phase2a.panel import _detected
    try:
        ok = bool(_detected(_fast_fns()["P14"](history, {})))
    except Exception:
        ok = False
    return {"emergent": ok, "em_score": 1.0 if ok else 0.0, "em_kind": "avalanche(P14-screen)"}


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
    from epc.phase2a.panel import _detected, _verdict
    B, FULL = _load(); FAST = _fast_fns()
    h1, h2 = _thin(history, FAST_FRAMES), _thin(history, CONFIRM_FRAMES)
    fired = []
    for pid, fn in FAST.items():
        try:
            res = fn(h1, metadata)
            if _detected(res): fired.append((pid, _verdict(res)))
        except Exception:
            pass
    confirmed = []
    for pid, _ in fired:
        try:
            res = FULL[pid](h2, metadata)
            if _detected(res): confirmed.append((pid, _verdict(res)))
        except Exception:
            pass
    confirmed.sort(key=lambda x: _TIER.get(x[1], 0), reverse=True)
    top = confirmed[0] if confirmed else None
    if top and _TIER.get(top[1], 0) >= 2:
        return {"verdict": "MATCH", "pattern": top[0], "tier": top[1], "fired": fired, "confirmed": confirmed}
    return {"verdict": "NO-MATCH", "pattern": None, "tier": top[1] if top else None, "fired": fired,
            "confirmed": confirmed}


def evaluate(out, prog, battery_all=False, battery=True):
    """per view: per-seed screen (+ fingerprint and, if battery, the known filter on emergent seeds), then the view
    verdict. battery=False is the cluster-first census mode: no per-program battery; views are EMERGENT or NONE and
    carry fingerprints for clustering."""
    from .fingerprint import fingerprint
    S = out["C"].shape[0] if "C" in out else next(v.shape[0] for k, v in out.items() if k != "meta" and hasattr(v, "shape"))
    per = {}
    for s in range(S):
        for name, hist, md in views(out, prog, s):
            thin = _screen_hist(name, hist, "adj" in out)
            r = screen(thin) if name != "C.aval" else screen_avalanche(hist)
            if r["emergent"] and name != "C.aval":
                r["fp"] = fingerprint(name, thin, adj0=out["adj0"][s].astype(np.int64) if "adj0" in out else None)
            if battery and (r["emergent"] or battery_all): r.update(known(hist, md))
            per.setdefault(name, []).append(r)
    res = {}
    for name, rs in per.items():
        em = sum(r["emergent"] for r in rs)
        pats = [r.get("pattern") for r in rs if r.get("verdict") == "MATCH"]
        top = max(set(pats), key=pats.count) if pats else None
        if top and pats.count(top) >= 2: v = "MATCH"
        elif em >= 2: v = "UNCLASSIFIED" if battery else "EMERGENT"
        else: v = "NONE"
        res[name] = {"verdict": v, "pattern": top if v == "MATCH" else None, "emergent_seeds": em, "seeds": rs}
    return res
