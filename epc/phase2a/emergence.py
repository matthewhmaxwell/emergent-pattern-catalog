"""Generic, substrate-agnostic emergence indicator (Milestone C, T2b).

The "none-of-the-above" primitive: a signal that *something* self-organizing is
happening, independent of any specific catalogued pattern. Combined with the
battery's per-pattern calibrated confidences, it yields the three-way verdict:
  - a specific detector fires            -> MATCH (that pattern)
  - emergence high, no detector fires    -> EMERGENT-UNCLASSIFIED  (the novel case)
  - emergence low,  no detector fires    -> NO-EMERGENCE

Principle (self-calibrating, like the detectors' own nulls, lifted to a generic
structural statistic): the late-window observation is (a) more structured than an
order-destroying shuffle of itself (null_z) AND (b) more structured than its own
early window (order_gain). The structure statistic is chosen by data kind.

COVERAGE (validated 2026-06-15, positive vs random-null separation):
  WORKS — field/grid spatial autocorrelation (P1/P3/P13), phase order (P9/P10
    Kuramoto r), point clustering (P5 flocking). E_pos ~0.75-0.98 vs E_null
    ~0.02-0.07.
  KNOWN GAPS — morphologies whose structure is NOT point-clustering / field-
    autocorrelation / phase-order in the FINAL frame: rotational (P6 milling),
    banded (P7 lanes), directional/gradient (P17), transient waves whose final
    state is uniform (P22 SIR), and distributional/abstract vectors needing a
    type-aware null (P28 wealth, P16 +/-1 state). A universal indicator over all
    morphologies needs a battery of generic structure measures (clustering +
    alignment + rotation + banding + gradient + temporal-wave + distributional
    concentration), each with an appropriate null -- a T2b continuation, not yet
    built. Use generic_emergence today for the WORKS morphologies; treat a LOW
    score on a GAP morphology as INCONCLUSIVE, not "no emergence".
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple

import numpy as np


def _frame_array(frame: Dict[str, Any]) -> Tuple[Optional[str], Optional[np.ndarray]]:
    if not isinstance(frame, dict):
        return None, None
    for fk in ("field", "grid"):
        if fk in frame:
            a = np.asarray(frame[fk], dtype=float)
            if a.ndim == 2:
                return "field", a
    if "positions" in frame:
        a = np.asarray(frame["positions"], dtype=float)
        if a.ndim == 2 and a.shape[1] >= 2:
            return "positions", a[:, :2]
    for pk in ("phases", "theta"):
        if pk in frame:
            return "phases", np.asarray(frame[pk], dtype=float).ravel()
    for k in ("opinions", "state", "fraction_on", "attendance", "x", "task_assignments", "wealth"):
        if k in frame:
            a = np.asarray(frame[k], dtype=float).ravel()
            if a.size >= 1:
                return "vector", a
    return None, None


def _moran(field: np.ndarray) -> float:
    f = field - field.mean()
    denom = float((f * f).sum())
    if denom <= 1e-12:
        return 0.0
    num = float((f[:-1, :] * f[1:, :]).sum() + (f[:, :-1] * f[:, 1:]).sum())
    n_pairs = f[:-1, :].size + f[:, :-1].size
    return (num / n_pairs) / (denom / f.size)


def _clustering(pos: np.ndarray) -> float:
    # inverse mean nearest-neighbour distance, scale-normalised by domain extent
    n = pos.shape[0]
    if n < 3:
        return 0.0
    d = np.sqrt(((pos[:, None, :] - pos[None, :, :]) ** 2).sum(-1))
    np.fill_diagonal(d, np.inf)
    nn = d.min(1)
    extent = max((float(np.ptp(pos[:, 0])) + float(np.ptp(pos[:, 1]))) / 2.0, 1e-9)
    return float(extent / (nn.mean() + 1e-9))


def _order_vector(v: np.ndarray, bins: int = 16) -> float:
    h, _ = np.histogram(v, bins=bins)
    p = h / max(h.sum(), 1)
    nz = p[p > 0]
    ent = -(nz * np.log(nz)).sum() / np.log(bins)
    return 1.0 - float(ent)            # concentration = order


def _order_phases(ph: np.ndarray) -> float:
    return float(abs(np.exp(1j * ph).mean()))   # Kuramoto r


def _order_orient(theta: np.ndarray) -> float:
    """Orientational order from a set of headings: max of polar and nematic.

    polar  = |<e^{iθ}>|   (1 for a common direction — flocking)
    nematic= |<e^{2iθ}>|  (1 for a common AXIS, head-tail symmetric — nematic
             alignment, where +θ and -θ cancel under the polar measure)
    Taking the max catches both polar swarms and apolar (nematic) order, which a
    clustering/position measure cannot see. Random headings → ~1/sqrt(N).
    """
    th = np.asarray(theta, dtype=float).ravel()
    if th.size == 0:
        return 0.0
    polar = abs(np.exp(1j * th).mean())
    nematic = abs(np.exp(2j * th).mean())
    return float(max(polar, nematic))


def _local_phase_order(theta: np.ndarray, w: int = 5) -> float:
    """Mean LOCAL phase coherence on a 2-D phase grid: <|<e^{iθ}>_neighbourhood|>.
    High for locally-aligned phase fields (twisted states, traveling/spiral waves) even
    when GLOBAL Kuramoto r ~ 0; ~1/sqrt(window) for incoherent (uncorrelated) phases.
    This is the local-vs-global order distinction (the P33 active-nematic lesson) for
    oscillator lattices, which the global-r phase channel (_order_phases on the raveled
    field) cannot see."""
    t = np.asarray(theta, dtype=float)
    if t.ndim != 2 or t.size < 9:
        return 0.0
    try:
        from scipy.ndimage import uniform_filter
    except Exception:
        return 0.0
    c = uniform_filter(np.cos(t), w, mode="wrap")
    s = uniform_filter(np.sin(t), w, mode="wrap")
    return float(np.sqrt(c ** 2 + s ** 2).mean())


def _phase_grid_from_frame(frame: Dict[str, Any]) -> Optional[np.ndarray]:
    """The 2-D phase grid from a frame ('phases'/'theta'), or None. Distinct from
    _frame_array, which RAVELS phases (losing the spatial structure local order needs).
    Returns None for 1-D phase arrays (unstructured oscillators) -> channel stays inert."""
    if not isinstance(frame, dict):
        return None
    for pk in ("phases", "theta"):
        if pk in frame:
            a = np.asarray(frame[pk], dtype=float)
            return a if a.ndim == 2 else None
    return None


def _structure(kind: str, a: np.ndarray) -> float:
    if kind == "field":
        return _moran(a)
    if kind == "positions":
        return _clustering(a)
    if kind == "phases":
        return _order_phases(a)
    if kind == "orientation":
        return _order_orient(a)
    if kind == "local_phase":
        return _local_phase_order(a)
    return _order_vector(a)


def _shuffle(kind: str, a: np.ndarray, rng) -> np.ndarray:
    if kind == "orientation":
        return rng.uniform(-np.pi, np.pi, size=a.shape)
    if kind == "field":
        flat = a.ravel().copy(); rng.shuffle(flat); return flat.reshape(a.shape)
    if kind == "positions":
        lo = a.min(0); hi = a.max(0)
        return rng.uniform(lo, hi, size=a.shape)
    if kind == "phases":
        return rng.uniform(-np.pi, np.pi, size=a.shape)
    if kind == "local_phase":                       # destroy local order, keep value set
        flat = a.ravel().copy(); rng.shuffle(flat); return flat.reshape(a.shape)
    flat = a.copy(); rng.shuffle(flat); return flat


def _orientation_from_frame(frame: Dict[str, Any]) -> Optional[np.ndarray]:
    """Per-agent headings θ from a velocity field, if present. Returns None when
    the frame carries no usable velocity/heading information."""
    if not isinstance(frame, dict):
        return None
    for vk in ("velocities", "velocity", "vel"):
        if vk in frame:
            v = np.asarray(frame[vk], dtype=float)
            if v.ndim == 2 and v.shape[1] >= 2:
                speed = np.hypot(v[:, 0], v[:, 1])
                if np.any(speed > 1e-9):
                    return np.arctan2(v[:, 1], v[:, 0])
    for hk in ("headings", "orientation", "orientations", "angles"):
        if hk in frame:
            h = np.asarray(frame[hk], dtype=float).ravel()
            if h.size:
                return h
    return None


def _score_series(kind: str, series: List[np.ndarray], rng,
                  n_null: int) -> Dict[str, float]:
    """Emergence score for one structure channel: late-window structure vs an
    order-destroying shuffle (null_z) AND vs the early window (order_gain)."""
    w = max(1, len(series) // 5)
    early = float(np.mean([_structure(kind, a) for a in series[:w]]))
    late = float(np.mean([_structure(kind, a) for a in series[-w:]]))
    late_arr = series[-1]
    nulls = np.array([_structure(kind, _shuffle(kind, late_arr, rng)) for _ in range(n_null)])
    nmu, nsd = float(nulls.mean()), float(nulls.std())
    null_z = (late - nmu) / nsd if nsd > 1e-9 else (10.0 if late - nmu > 1e-6 else 0.0)
    order_gain = late - early
    z_term = 1.0 / (1.0 + np.exp(-(null_z - 2.0)))          # ~0 below 2 SD, ~1 above
    g_term = 1.0 / (1.0 + np.exp(-(order_gain) * 8.0))       # >0.5 when structure grew
    score = float(z_term * (0.5 + 0.5 * g_term))             # null-excess gated, growth-weighted
    return {"score": round(score, 4), "order_gain": round(order_gain, 4),
            "null_z": round(float(null_z), 3)}


def generic_emergence(history: List[Dict[str, Any]], n_null: int = 50,
                      seed: int = 0) -> Dict[str, Any]:
    """Return {score in [0,1], kind, order_gain, null_z, n_frames}.

    Scores the dominant spatial/phase/vector channel, AND — when the frames carry
    velocity/heading data — an ORIENTATION channel (polar + nematic order), taking
    the stronger of the two. The orientation channel catches apolar/nematic order
    that a clustering measure misses (rods that align without clustering), and can
    only RAISE the score (no null carries velocities, so no false-novel risk)."""
    rng = np.random.default_rng(seed)
    if not history:
        return {"score": 0.0, "kind": None, "order_gain": 0.0, "null_z": 0.0, "n_frames": 0}
    frames = [(_frame_array(f)) for f in history]
    frames = [(k, a) for k, a in frames if k is not None]
    if len(frames) >= 2:
        kind = frames[-1][0]
        same = [a for k, a in frames if k == kind]
        prim = _score_series(kind, same, rng, n_null)
        best = {**prim, "kind": kind, "n_frames": len(same)}
        # Orientation channel (velocity/heading-bearing systems)
        oris = [_orientation_from_frame(f) for f in history]
        oris = [o for o in oris if o is not None]
        if len(oris) >= 2:
            orie = _score_series("orientation", oris, np.random.default_rng(seed + 1), n_null)
            if orie["score"] > best["score"]:
                best = {**orie, "kind": "orientation", "n_frames": len(oris)}
    else:
        # No recognized morphology observable (e.g. a sandpile that emits only
        # 'activity'/'avalanche_sizes'). Fall through — the info-theoretic channels
        # below, especially heavy-tail which reads raw frames, may still detect it.
        kind = None
        best = {"score": 0.0, "kind": None, "order_gain": 0.0, "null_z": 0.0, "n_frames": 0}

    # Local phase-order channel — locally-coherent / globally-incoherent phase fields
    # (twisted states, traveling/spiral waves) that the GLOBAL Kuramoto-r phase channel
    # MISSES (it reads the RAVELED field, so it sees only global coherence). The P33
    # local-vs-global order lesson, for oscillator lattices. Reads the 2-D phase grid;
    # scores local order vs a phase-shuffle null (shuffling destroys local order -> no
    # phase null carries local structure, so no false-novel risk). The unstructured 1-D
    # uncoupled-phase null is 1-D, so _phase_grid_from_frame returns None -> inert there.
    pgrids = [g for g in (_phase_grid_from_frame(f) for f in history) if g is not None]
    if len(pgrids) >= 2:
        lpo = _score_series("local_phase", pgrids, np.random.default_rng(seed + 2), n_null)
        if lpo["score"] > best["score"]:
            best = {**lpo, "kind": "local-phase-order", "n_frames": len(pgrids)}

    # Synergy / causal-emergence channel (Psi_CE): catches "greater-than-the-parts"
    # collective effects (XOR-type, connectivity) that the morphology channels miss.
    # Psi_CE is intrinsically specific (it subtracts per-part redundancy, so noise
    # scores negative; the blind-spot audit measured zero null false-positives), so
    # a small positive threshold is used and it can only RAISE the score.
    try:
        from epc.phase2a.info_channels import (graph_structure_score, heavy_tail_score,
                                                micro_macro, oscillation_score, psi_ce_best)
        M, cands = micro_macro(history)
        if M is not None:
            psi, feat = psi_ce_best(M, cands)
            # Synergy / causal-emergence channel. Only contribute when Psi_CE is
            # CLEARLY positive (gate at 0.08; nulls score ~-1.7) to avoid a floor
            # artifact at Psi_CE≈0 for static systems.
            if psi == psi and psi > 0.08:
                syn = float(min(0.999, 0.55 + 0.44 * (1.0 - np.exp(-(psi - 0.08) * 6.0))))
                if syn > best["score"]:
                    best = {"score": round(syn, 4),
                            "order_gain": round(float(psi), 4), "null_z": 0.0,
                            "kind": f"synergy(psi_ce:{feat})", "n_frames": M.shape[0]}
            # Temporal channel — sustained oscillation / limit cycle / propagating
            # or chaotic dynamics, via spectral peak-to-mean of the collective
            # signal. GATED OFF for phase-kind data: phase synchronization is
            # already handled by the order-parameter channel, and a raw wrapped-
            # phase mean produces a spurious peak. Threshold 15 sits well above the
            # non-phase null maximum (~9, Brownian drift) and below the emergent
            # range (limit cycle ~49, chaos ~64, front ~35).
            if kind != "phases":
                osc = oscillation_score(M.mean(1))
                if osc > 15.0:
                    oscore = float(min(0.95, 0.55 + 0.4 * (1.0 - np.exp(-(osc - 15.0) * 0.04))))
                    if oscore > best["score"]:
                        best = {"score": round(oscore, 4),
                                "order_gain": round(float(osc), 4), "null_z": 0.0,
                                "kind": "temporal(spectral-peak)", "n_frames": M.shape[0]}
        # Heavy-tail / self-organized-criticality channel. Reads RAW frames (size
        # observables like 'activity'/'avalanche_sizes' are invisible to the
        # spatial/phase channels — a sandpile scores 0 on all of them). Fires when
        # the event-size distribution is power-law (LLR>0) over ≥1.3 decades.
        ht = heavy_tail_score(history)
        if ht is not None:
            llr, decades = ht
            if llr > 0.0 and decades >= 1.3:
                hts = float(min(0.95, 0.55 + 0.2 * (decades - 1.3)))
                if hts > best["score"]:
                    best = {"score": round(hts, 4), "order_gain": round(float(decades), 3),
                            "null_z": 0.0, "kind": "heavy-tail(power-law)",
                            "n_frames": len(history)}
        # Network-topology channel (reads an 'adjacency' observable — a substrate
        # the spatial/phase channels cannot ingest): emergent community structure
        # (spectral modularity vs an ER null) or a scale-free degree distribution.
        gs = graph_structure_score(history)
        if gs is not None:
            zmod, qmod, scale_free = gs
            if zmod > 3.0 or scale_free:
                gscore = float(min(0.95, 0.6 + 0.05 * max(zmod - 3.0, 0.0)))
                if scale_free:
                    gscore = max(gscore, 0.7)
                if gscore > best["score"]:
                    best = {"score": round(gscore, 4), "order_gain": round(float(qmod), 3),
                            "null_z": round(float(zmod), 3),
                            "kind": "network(modularity/scale-free)", "n_frames": len(history)}
        # Bond-orientational (hexatic) channel: emergent local psi6 order from a
        # disordered point set (entropy-driven crystallization). Reads 'positions';
        # fires only on a GAIN in local hexatic order reaching a hexatic floor. A
        # crystal has near-uniform density (the clustering channel misses it), and
        # absolute psi6~0.5 is shared with random close packing — only the GAIN is
        # specific. Verified to NOT fire on flocking, disordered gas, dense clumps
        # (Keller-Segel), or random points (gain <= +0.003 there).
        from epc.metrics.positional_order import local_psi6
        pos_frames = [np.asarray(f["positions"]) for f in history
                      if isinstance(f, dict) and "positions" in f]
        if len(pos_frames) >= 8:
            allp = np.vstack(pos_frames)
            Lh = float(allp.max() - allp.min()) or None
            e6 = float(np.mean([local_psi6(p, Lh) for p in pos_frames[1:4]]))
            l6 = float(np.mean([local_psi6(p, Lh) for p in pos_frames[-5:]]))
            if l6 >= 0.42 and (l6 - e6) > 0.04:
                hx = float(min(0.95, 0.55 + 2.0 * ((l6 - e6) - 0.04)))
                if hx > best["score"]:
                    best = {"score": round(hx, 4), "order_gain": round(float(l6 - e6), 4),
                            "null_z": 0.0, "kind": "hexatic(psi6-gain)",
                            "n_frames": len(pos_frames)}
        # Conserved-scalar power-law channel: a Pareto (power-law) tail in the
        # CROSS-SECTIONAL distribution of a conserved scalar resource ('wealth'),
        # emergent from agent heterogeneity (P36). The heavy-tail/SOC channel above
        # reads event sizes OVER TIME and misses this; the spatial/phase channels
        # cannot ingest a bare value vector. Fires only on a clean power-law tail
        # (Pareto index ~1, power-law preferred over exponential, not condensed), so
        # the Gamma / Boltzmann / Yard-Sale-condensation / null look-alikes do NOT fire.
        from epc.metrics.wealth_concentration import (
            gini as _gini_w, hill_tail_alpha as _hill_a,
            pareto_ks_distance as _pks, exp_ks_distance as _eks)
        wfr = [np.asarray(f["wealth"], dtype=float) for f in history
               if isinstance(f, dict) and "wealth" in f]
        wfr = [w for w in wfr if w.size >= 100 and w.sum() > 0]
        if len(wfr) >= 4:
            late_w = wfr[(3 * len(wfr)) // 4:]
            a_, ksp_, kse_, ms_ = [], [], [], []
            for w in late_w:
                al, wmin, _ = _hill_a(w)
                ms_.append(float(w.max() / w.sum()))
                if np.isfinite(al) and np.isfinite(wmin) and wmin > 0:
                    a_.append(al); ksp_.append(_pks(w, al, wmin)); kse_.append(_eks(w, wmin))
            if a_:
                al = float(np.median(a_)); ksp = float(np.nanmedian(ksp_))
                adv = float(np.nanmedian(kse_)) - ksp; ms = float(np.mean(ms_))
                gw = float(np.mean([_gini_w(w) for w in late_w]))
                if (0.6 <= al <= 1.8) and adv > 0.12 and ms < 0.20 and gw > 0.45 and ksp < 0.08:
                    pw = float(min(0.95, 0.55 + 1.2 * (adv - 0.12)))
                    if pw > best["score"]:
                        best = {"score": round(pw, 4), "order_gain": round(adv, 4),
                                "null_z": 0.0, "kind": "heavy-tail(pareto-wealth)",
                                "n_frames": len(wfr)}
        # Coexistence-oscillation channel: sustained multi-species oscillation under
        # field-mediated resource competition (P37, Huisman-Weissing). Reads
        # 'abundances'+'resources'; fires only on multi-species coexistence (>=3 persist)
        # whose abundance vector keeps oscillating (does not settle) -- invisible to the
        # spatial/phase/temporal channels, which do not ingest an abundance+resource
        # series. Stable coexistence (settles -> osc~0) and exclusion (few survivors) do NOT fire.
        from epc.metrics.competition_dynamics import coexistence_oscillation as _coex
        _cm = _coex(history)
        if _cm is not None and int(_cm["n_persist"]) >= 3 and _cm["osc_cv"] > 0.05:
            cs = float(min(0.95, 0.55 + 2.0 * (_cm["osc_cv"] - 0.05)))
            if cs > best["score"]:
                best = {"score": round(cs, 4), "order_gain": round(float(_cm["osc_cv"]), 4),
                        "null_z": 0.0, "kind": "coexistence-oscillation",
                        "n_frames": len(history)}
    except Exception:
        pass

    return {"score": best["score"], "kind": best["kind"],
            "order_gain": best["order_gain"], "null_z": best["null_z"],
            "n_frames": best["n_frames"]}


def three_way_verdict(top_calibrated_confidence: Optional[float],
                      emergence_score: float,
                      present_threshold: float = 0.90,
                      emergence_threshold: float = 0.50) -> Dict[str, Any]:
    """Combine the battery's best per-pattern calibrated confidence with the
    generic emergence score into the T2b three-way verdict.

      MATCH                 - a specific detector is calibrated 'present' (>= thr)
      EMERGENT-UNCLASSIFIED - no detector reaches 'present' but emergence is high
                              (the novel-pattern / discovery signal)
      NO-EMERGENCE          - no detector fires and emergence is low

    A low emergence_score is only trustworthy for the morphologies the indicator
    covers (see module COVERAGE); on a GAP morphology a non-MATCH with low
    emergence should be surfaced as INCONCLUSIVE by the caller.
    """
    tc = top_calibrated_confidence if top_calibrated_confidence is not None else 0.0
    if tc >= present_threshold:
        verdict = "MATCH"
    elif emergence_score >= emergence_threshold:
        verdict = "EMERGENT-UNCLASSIFIED"
    else:
        verdict = "NO-EMERGENCE"
    return {"verdict": verdict, "top_confidence": round(tc, 4),
            "emergence_score": round(float(emergence_score), 4)}
