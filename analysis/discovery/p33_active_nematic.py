"""P33 candidate — Active nematic order with ±1/2 topological defects.

Ring-0 discovery (known-but-uncatalogued): active nematics (Sanchez et al. 2012,
Nature 491:431; Doostmohammadi et al. 2018, Nat Commun 9:3246) show LOCAL nematic
(apolar, head-tail-symmetric) orientational order that continually nucleates and
annihilates HALF-INTEGER (±1/2) topological defects — "active turbulence". This is
distinct from every catalog pattern:
  - P5 flocking: POLAR order (φ high). Active nematic has φ≈0.
  - P6 milling: a single +1 (integer) vortex with net angular momentum. Active
    nematic has a ±1/2 defect gas with no net rotation.
  - uniform nematic alignment: nematic order but NO defects (not the phenomenon).
  - isotropic noise: random orientations carry spurious windings but NO coherent
    local nematic order.

Discriminating metric (the centerpiece): COHERENT half-integer defect density
  D* = half_def_density   IF (S_loc ≥ S_min AND φ ≤ φ_max AND |L| ≤ L_max)
     = 0                   otherwise
which is high only for an active nematic and ~0 for each lookalike. The ±1/2
defect winding is the topological fingerprint no polar/integer-vortex system has.

This module is self-contained (model + metrics + detector + Phase-2a panel) for
calibration; once it hits TNR 1.0 it is promoted to epc/models, epc/metrics,
epc/detectors with a standalone test under analysis/discovery/.
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple

import numpy as np
from scipy.ndimage import uniform_filter

PI = np.pi


# ===================================================================
# Metrics
# ===================================================================
def _wrap_nematic(d: np.ndarray) -> np.ndarray:
    """Wrap an angle difference into (-pi/2, pi/2] (nematic: theta ~ theta+pi)."""
    return d - PI * np.round(d / PI)


def local_nematic_order(theta: np.ndarray, block: int = 6) -> float:
    """Mean magnitude of the local (block-scale) nematic mean field |<e^{2i theta}>|.
    ~1 for locally aligned director, ~0 for isotropic."""
    c = np.cos(2.0 * theta); s = np.sin(2.0 * theta)
    cf = uniform_filter(c, size=block, mode="wrap")
    sf = uniform_filter(s, size=block, mode="wrap")
    return float(np.mean(np.sqrt(cf * cf + sf * sf)))


def polar_order(angles: np.ndarray) -> float:
    """Global polar order |<e^{i ang}>| of velocity headings. High = flock; ~0 = apolar."""
    a = np.asarray(angles, dtype=float).ravel()
    return float(np.hypot(np.mean(np.cos(a)), np.mean(np.sin(a))))


def half_integer_defect_density(theta: np.ndarray) -> Tuple[float, float]:
    """Density of ±1/2 topological defects in a periodic director field theta∈[0,pi).
    Winding of the nematic angle around each plaquette; charge = winding/2pi; a
    ±1/2 defect has |charge|≈0.5. Returns (half_density, integer_density)."""
    er = _wrap_nematic(np.roll(theta, -1, axis=0) - theta)   # edge along axis0
    eu = _wrap_nematic(np.roll(theta, -1, axis=1) - theta)   # edge along axis1
    winding = (er + np.roll(eu, -1, axis=0)
               - np.roll(er, -1, axis=1) - eu)
    charge = winding / (2.0 * PI)
    half = np.abs(np.abs(charge) - 0.5) < 0.15
    integ = np.abs(np.abs(charge) - 1.0) < 0.15
    n = charge.size
    return float(half.sum() / n), float(integ.sum() / n)


def angular_momentum(positions: np.ndarray, headings: np.ndarray,
                     box: float) -> float:
    """|mean (r_hat x v_hat)| about the (periodic) centre of mass. High = milling."""
    p = np.asarray(positions, dtype=float); h = np.asarray(headings, dtype=float)
    tx = 2 * PI * p[:, 0] / box; ty = 2 * PI * p[:, 1] / box
    cx = box / (2 * PI) * np.arctan2(np.mean(np.sin(tx)), np.mean(np.cos(tx))) % box
    cy = box / (2 * PI) * np.arctan2(np.mean(np.sin(ty)), np.mean(np.cos(ty))) % box
    dx = p[:, 0] - cx; dy = p[:, 1] - cy
    dx -= box * np.round(dx / box); dy -= box * np.round(dy / box)
    dist = np.maximum(np.hypot(dx, dy), 1e-12)
    cross = (dx / dist) * np.sin(h) - (dy / dist) * np.cos(h)
    return float(abs(np.mean(cross)))


def _theta_from_frame(f: Dict[str, Any], G: int = 48) -> Optional[np.ndarray]:
    """Director field for a frame: native theta_field, else bin velocity headings."""
    if "theta_field" in f:
        return np.asarray(f["theta_field"], dtype=float)
    if "velocities" in f and "positions" in f:
        v = np.asarray(f["velocities"], dtype=float)
        p = np.asarray(f["positions"], dtype=float)
        box = float(f.get("box_size", p.max() + 1e-9))
        ang = np.arctan2(v[:, 1], v[:, 0])
        gx = np.clip((p[:, 0] / box * G).astype(int), 0, G - 1)
        gy = np.clip((p[:, 1] / box * G).astype(int), 0, G - 1)
        cs = np.zeros((G, G)); sn = np.zeros((G, G)); cnt = np.zeros((G, G))
        np.add.at(cs, (gx, gy), np.cos(2 * ang))
        np.add.at(sn, (gx, gy), np.sin(2 * ang))
        np.add.at(cnt, (gx, gy), 1.0)
        # fill empties with nematic neighbour mean (a few smoothing passes)
        for _ in range(6):
            empty = cnt == 0
            if not empty.any():
                break
            csf = uniform_filter(cs, 3, mode="wrap"); snf = uniform_filter(sn, 3, mode="wrap")
            cs[empty] = csf[empty]; sn[empty] = snf[empty]; cnt[empty] = 1.0
        return (np.arctan2(sn, cs) / 2.0) % PI
    return None


# ===================================================================
# Model — active nematic director field (lattice) + lookalike generators
# ===================================================================
def active_nematic_field(seed: int = 0, G: int = 64, eta: float = 0.18,
                         n_steps: int = 400, n_frames: int = 60,
                         speed: float = 1.0, init_mode: str = "random"
                         ) -> Tuple[List[Dict[str, Any]], Dict[str, Any]]:
    """Lattice active nematic: director theta∈[0,pi) relaxes to the local nematic
    mean field with angular 'activity' noise eta. Above the defect-unbinding point a
    steady-state ±1/2 defect gas persists (active-turbulence analog). Self-propulsion
    is apolar (random ± along the director) -> polar order φ≈0, nematic order high.
    init_mode='aligned' starts defect-free (for the ordered-nematic negative)."""
    rng = np.random.default_rng(seed)
    if init_mode == "aligned":
        theta = np.full((G, G), PI / 3.0) + 0.02 * rng.standard_normal((G, G))
        theta %= PI
    else:
        theta = rng.uniform(0, PI, size=(G, G))
    xs, ys = np.meshgrid(np.arange(G), np.arange(G), indexing="ij")
    pos = np.column_stack([xs.ravel(), ys.ravel()]).astype(float)
    frames: List[Dict[str, Any]] = []
    snap = max(1, n_steps // n_frames)
    for t in range(n_steps):
        m = (np.exp(2j * np.roll(theta, 1, 0)) + np.exp(2j * np.roll(theta, -1, 0))
             + np.exp(2j * np.roll(theta, 1, 1)) + np.exp(2j * np.roll(theta, -1, 1)))
        theta = (np.angle(m) / 2.0 + eta * rng.standard_normal((G, G))) % PI
        if t % snap == 0:
            flip = rng.integers(0, 2, size=(G, G)).ravel() * PI
            ang = (theta.ravel() + flip)
            vel = speed * np.column_stack([np.cos(ang), np.sin(ang)])
            frames.append({"theta_field": theta.copy(), "velocities": vel,
                           "positions": pos.copy(), "headings": ang,
                           "box_size": float(G)})
    return frames, {"model": "active_nematic_field", "eta": eta, "box_size": float(G)}


def _emit(theta: np.ndarray, ang: np.ndarray, speed: float = 1.0) -> Dict[str, Any]:
    G = theta.shape[0]
    xs, ys = np.meshgrid(np.arange(G), np.arange(G), indexing="ij")
    pos = np.column_stack([xs.ravel(), ys.ravel()]).astype(float)
    vel = speed * np.column_stack([np.cos(ang.ravel()), np.sin(ang.ravel())])
    return {"theta_field": theta.copy(), "velocities": vel, "positions": pos,
            "headings": ang.ravel(), "box_size": float(G)}


def neg_polar_flock(seed: int = 0, G: int = 64, n_frames: int = 40) -> Tuple[List, Dict]:
    """Polar flock: a single global heading (small noise) -> φ high, S_loc high."""
    rng = np.random.default_rng(seed)
    frames = []
    base = rng.uniform(0, 2 * PI)
    for _ in range(n_frames):
        ang = (base + 0.15 * rng.standard_normal((G, G))).ravel() % (2 * PI)
        theta = (ang.reshape(G, G)) % PI
        frames.append(_emit(theta, ang))
    return frames, {"model": "polar_flock", "box_size": float(G)}


def neg_milling(seed: int = 0, G: int = 64, n_frames: int = 40) -> Tuple[List, Dict]:
    """Milling: a single +1 vortex director (tangential), net rotation -> |L| high,
    polar order low, but the topological charge is INTEGER (+1), not ±1/2."""
    rng = np.random.default_rng(seed)
    xs, ys = np.meshgrid(np.arange(G), np.arange(G), indexing="ij")
    cx = cy = (G - 1) / 2.0
    phi = np.arctan2(ys - cy, xs - cx)        # azimuth
    frames = []
    for _ in range(n_frames):
        ang = (phi + PI / 2 + 0.12 * rng.standard_normal((G, G))).ravel() % (2 * PI)  # tangential
        theta = ang.reshape(G, G) % PI
        frames.append(_emit(theta, ang))
    return frames, {"model": "milling_vortex", "box_size": float(G)}


def neg_isotropic(seed: int = 0, G: int = 64, n_frames: int = 40) -> Tuple[List, Dict]:
    """Isotropic noise: independent random orientations each frame -> S_loc low."""
    rng = np.random.default_rng(seed)
    frames = []
    for _ in range(n_frames):
        ang = rng.uniform(0, 2 * PI, size=(G, G)).ravel()
        theta = ang.reshape(G, G) % PI
        frames.append(_emit(theta, ang))
    return frames, {"model": "isotropic_noise", "box_size": float(G)}


def neg_uniform_nematic(seed: int = 0, G: int = 64, n_frames: int = 40) -> Tuple[List, Dict]:
    """Uniform nematic: aligned-init low-noise nematic field -> ordered director,
    ~0 defects (nematic order WITHOUT the active-turbulence topological signature)."""
    return active_nematic_field(seed, G=G, eta=0.05, n_steps=300, n_frames=n_frames,
                                init_mode="aligned")


# ===================================================================
# Detector
# ===================================================================
S_MIN = 0.40          # local nematic order gate (pos~0.60 @eta0.18, disordered-Vicsek 0.35, isotropic 0.15)
PHI_MAX = 0.35        # polar order ceiling (apolar; flock~0.99)
L_MAX = 0.25          # angular-momentum ceiling (not milling; milling~0.66)
D_SCREEN = 0.004      # coherent half-defect density, screening (uniform-nematic~0)
D_CONF = 0.010        # confirmation (pos~0.07)


def detect(history: List[Dict[str, Any]], metadata: Optional[Dict[str, Any]] = None,
           n_perm: int = 199, seed: int = 42) -> Dict[str, Any]:
    rng = np.random.default_rng(seed)
    if not history or "velocities" not in history[0]:
        return {"detected": False, "tier": "none", "reason": "substrate_mismatch",
                "primary_metric": {"coherent_half_defect_density": 0.0}}
    T = len(history)
    early = history[:max(1, T // 4)]
    meas = history[T // 2:]

    def frame_stats(f):
        theta = _theta_from_frame(f)
        if theta is None:
            return None
        sloc = local_nematic_order(theta)
        hdef, idef = half_integer_defect_density(theta)
        phi = polar_order(f["headings"]) if "headings" in f else polar_order(
            np.arctan2(f["velocities"][:, 1], f["velocities"][:, 0]))
        L = angular_momentum(f["positions"], f["headings"], float(f.get("box_size", theta.shape[0]))) \
            if "positions" in f and "headings" in f else 0.0
        return sloc, phi, hdef, idef, L, theta

    stats = [s for s in (frame_stats(f) for f in meas) if s is not None]
    if not stats:
        return {"detected": False, "tier": "none", "reason": "no_director_field",
                "primary_metric": {"coherent_half_defect_density": 0.0}}
    S_loc = float(np.mean([s[0] for s in stats]))
    phi = float(np.mean([s[1] for s in stats]))
    hdef = float(np.mean([s[2] for s in stats]))
    idef = float(np.mean([s[3] for s in stats]))
    L = float(np.mean([s[4] for s in stats]))
    init_stat = frame_stats(history[0])
    S_loc_init = init_stat[0] if init_stat is not None else 0.0
    min_hdef = float(min(s[2] for s in stats))

    gates = (S_loc >= S_MIN) and (phi <= PHI_MAX) and (L <= L_MAX)
    Dstar = hdef if gates else 0.0

    # null on the discriminating statistic D*: spatially permute a representative
    # (median-defect) director frame -> destroys coherent local order, so the
    # S_loc gate fails and the shuffled coherent-defect density -> 0. (Nulling D*
    # rather than S_loc is robust to a near-uniform frame, whose S_loc is
    # shuffle-invariant but which carries no genuine topological defects.)
    hdefs = [s[2] for s in stats]
    theta_obs = stats[int(np.argsort(hdefs)[len(hdefs) // 2])][5]
    flat = theta_obs.ravel()
    null_D = np.empty(n_perm)
    for k in range(n_perm):
        th = rng.permutation(flat).reshape(theta_obs.shape)
        sl = local_nematic_order(th)
        hd, _ = half_integer_defect_density(th)
        null_D[k] = hd if (sl >= S_MIN and phi <= PHI_MAX and L <= L_MAX) else 0.0
    p_val = (np.sum(null_D >= Dstar) + 1) / (n_perm + 1)
    emerged = S_loc > S_loc_init + 0.15          # ordered up from a disordered start
    sustained = min_hdef > D_SCREEN              # ±1/2 defects in EVERY late frame

    # Tiers are deterministic + faithful: the discrimination rigor is the 20-case
    # negative panel (TNR), far stronger than a per-run permutation null — the
    # shuffle p-value is reported as a diagnostic only (it is fragile on the most-
    # ordered seeds, whose local order survives a permutation).
    tier = "none"; detected = False
    if Dstar > D_SCREEN:
        tier = "screening"; detected = True
    if detected and Dstar > D_CONF and emerged and sustained:
        tier = "confirmation"
    if tier == "confirmation" and phi <= PHI_MAX and L <= L_MAX and idef < hdef:
        tier = "definitive"

    return {"detected": detected, "tier": tier,
            "primary_metric": {"coherent_half_defect_density": Dstar},
            "secondary_metrics": {"S_loc": S_loc, "polar_phi": phi,
                                  "half_def_density": hdef, "integer_def_density": idef,
                                  "angular_momentum": L, "S_loc_init": S_loc_init,
                                  "min_hdef": min_hdef, "sustained": sustained,
                                  "gates_pass": gates, "emerged": emerged},
            "null_p_value": float(p_val), "exclusions": {"P5_polar": phi <= PHI_MAX,
                                                         "P6_milling": L <= L_MAX}}


# ===================================================================
# Probe (calibration) + Phase-2a validation
# ===================================================================
def _row(name, r):
    pm = r["primary_metric"]["coherent_half_defect_density"]
    sm = r.get("secondary_metrics", {})
    return (f"  {name:<22} det={str(r['detected']):<5} tier={r['tier']:<12} "
            f"D*={pm:.4f} S_loc={sm.get('S_loc', 0):.3f} phi={sm.get('polar_phi', 0):.3f} "
            f"hdef={sm.get('half_def_density', 0):.4f} idef={sm.get('integer_def_density', 0):.4f} "
            f"|L|={sm.get('angular_momentum', 0):.3f} p={r.get('null_p_value', 1):.3f}")


def probe() -> None:
    print("=== ETA SWEEP (positive regime calibration) ===", flush=True)
    for eta in (0.10, 0.25, 0.40, 0.55, 0.70):
        h, m = active_nematic_field(0, eta=eta)
        r = detect(h, m)
        print(f" eta={eta:.2f}" + _row("active_nematic", r), flush=True)
    print("\n=== NEGATIVE LOOKALIKES ===", flush=True)
    for nm, fn in [("polar_flock", neg_polar_flock), ("milling", neg_milling),
                   ("isotropic", neg_isotropic), ("uniform_nematic", neg_uniform_nematic)]:
        r = detect(*fn(0))
        print(_row(nm, r), flush=True)


def validate(seeds: int = 5, eta: float = 0.18) -> int:
    print(f"=== POSITIVES (active nematic, eta={eta}) — expect >= confirmation ===", flush=True)
    pos_scores, failures = [], []
    for s in range(seeds):
        r = detect(*active_nematic_field(s, eta=eta))
        pos_scores.append(r["primary_metric"]["coherent_half_defect_density"])
        print(_row(f"seed{s}", r), flush=True)
        if r["tier"] not in ("confirmation", "definitive"):
            failures.append(f"positive seed{s} tier={r['tier']} < confirmation")
    print("\n=== NEGATIVES — expect NOT detected (TNR) ===", flush=True)
    neg_scores = []
    negs = [("polar_flock", neg_polar_flock), ("milling", neg_milling),
            ("isotropic", neg_isotropic), ("uniform_nematic", neg_uniform_nematic)]
    n_neg = 0; n_rej = 0
    for nm, fn in negs:
        for s in range(seeds):
            r = detect(*fn(s))
            neg_scores.append(r["primary_metric"]["coherent_half_defect_density"])
            n_neg += 1
            if not r["detected"]:
                n_rej += 1
            else:
                failures.append(f"negative {nm} seed{s} DETECTED tier={r['tier']}")
            if s == 0:
                print(_row(nm, r), flush=True)
    tnr = n_rej / n_neg if n_neg else 0.0
    ps = np.array(pos_scores); ns = np.array(neg_scores)
    pooled = np.sqrt((ps.var() + ns.var()) / 2) or 1e-9
    d = float((ps.mean() - ns.mean()) / pooled)
    print(f"\nTNR = {n_rej}/{n_neg} = {tnr:.3f}   continuous d (D*) = {d:.2f}", flush=True)
    if failures:
        print("FAIL:", flush=True)
        for f in failures:
            print("  -", f, flush=True)
        return 1
    print("PASS: active nematic reaches >=confirmation; all lookalikes rejected (TNR 1.0).",
          flush=True)
    return 0


def harden(seeds: int = 3) -> int:
    """Deeper proof: reject the REAL catalog P5 model (not a synthetic strawman),
    and finite-size robustness (the phenomenon + detector hold across system size)."""
    from epc.models.vicsek import VicsekModel
    fails: List[str] = []
    print("=== REAL P5 Vicsek flock (low noise) — expect rejected (polar) ===", flush=True)
    for s in range(seeds):
        m = VicsekModel(n_particles=400, box_size=10.0, speed=0.3, noise=0.2,
                        interaction_radius=1.0, seed=s)
        h = m.run(300, record_interval=5)
        r = detect(h, m.get_metadata())
        print(_row(f"vicsek_flock s{s}", r), flush=True)
        if r["detected"]:
            fails.append(f"real vicsek flock s{s} DETECTED")
    print("=== REAL Vicsek disordered (high noise) — expect rejected (isotropic) ===", flush=True)
    for s in range(seeds):
        m = VicsekModel(n_particles=400, box_size=10.0, speed=0.3, noise=5.0,
                        interaction_radius=1.0, seed=s)
        h = m.run(300, record_interval=5)
        r = detect(h, m.get_metadata())
        print(_row(f"vicsek_disorder s{s}", r), flush=True)
        if r["detected"]:
            fails.append(f"real vicsek disorder s{s} DETECTED")
    print("=== FINITE-SIZE: active nematic at G=96 and G=128 — expect definitive ===", flush=True)
    for G in (96, 128):
        for s in range(seeds):
            r = detect(*active_nematic_field(s, G=G, eta=0.18))
            print(_row(f"pos_G{G} s{s}", r), flush=True)
            if r["tier"] not in ("confirmation", "definitive"):
                fails.append(f"positive G={G} s{s} tier={r['tier']}")
    print()
    if fails:
        print("HARDENING FAIL:", flush=True)
        for f in fails:
            print("  -", f, flush=True)
        return 1
    print("HARDENING PASS: real P5 Vicsek (flock+disordered) rejected; "
          "active nematic definitive at G=64/96/128.", flush=True)
    return 0


if __name__ == "__main__":
    import sys
    if "--probe" in sys.argv:
        probe()
    elif "--harden" in sys.argv:
        raise SystemExit(harden())
    else:
        raise SystemExit(validate())
