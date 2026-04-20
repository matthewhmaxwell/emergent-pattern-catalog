"""Density-phase-separation metrics for P2 (MIPS) detection.

Primitives:
  - particle_local_density: count neighbors in r_cg disk, divide by area
  - phase_coexistence_fractions: (f_gas, f_liquid) at a rho_star threshold
  - two_phase_coexistence_score: min(f_gas, f_liquid) — P2 primary metric
  - density_speed_anticorrelation: Pearson r(rho_local, |v|)
  - constant_speed_null_score: would the same positions + constant v
    (no density feedback) give a two_phase signal? Surrogate null.

All primitives are pure functions over (state_history, params). They do
not own any threshold logic — the detector (p2_mips.py) applies the
tiered thresholds.

Design notes:
  - f_gas uses rho < rho_star/2 as the "dilute" gate (Fily-Marchetti).
  - f_liquid uses rho > rho_star as the "dense/jammed" gate.
  - The two-phase score handles the Vicsek / D'Orsogna / stuck-ABP false
    positives empirically: Vicsek flocking gives f_liquid~1, f_gas~0 (no
    gas phase); stuck ABP gives f_liquid~1, f_gas~0 too; MIPS shows BOTH.
  - Local density is computed with PERIODIC BC when box_size is provided.
    If periodic is False, the cKDTree is built without boxsize and
    particles near boundaries underestimate density (boundary artifact).
    Downstream detector rejects non-periodic substrates by default.

Reference:
  Fily, Y. & Marchetti, M. C. (2012). Athermal Phase Separation of
  Self-Propelled Particles with No Alignment. PRL 108, 235702.
"""
from __future__ import annotations

from typing import Any, Optional

import numpy as np
from numpy.typing import NDArray
from scipy.spatial import cKDTree


def particle_local_density(
    positions: NDArray[np.float64],
    r_cg: float = 1.0,
    box_size: Optional[float] = None,
) -> NDArray[np.float64]:
    """Count neighbors in r_cg disk (including self), normalize by area.

    Parameters
    ----------
    positions : (N, 2) float array
    r_cg : coarse-graining radius (default 1.0 = particle diameter sigma).
    box_size : if not None, apply periodic BC with square box side L.

    Returns (N,) array of local number densities per unit area.
    """
    if positions.ndim != 2 or positions.shape[1] != 2:
        raise ValueError(
            f"positions must be (N, 2), got shape {positions.shape}"
        )
    if r_cg <= 0:
        raise ValueError(f"r_cg must be positive, got {r_cg}")
    if positions.shape[0] == 0:
        return np.zeros(0)
    if box_size is not None:
        tree = cKDTree(positions, boxsize=box_size)
    else:
        tree = cKDTree(positions)
    counts = tree.query_ball_point(positions, r=r_cg, return_length=True)
    counts = np.asarray(counts, dtype=np.float64)
    area = np.pi * r_cg**2
    return counts / area


def particle_speeds(velocities: NDArray[np.float64]) -> NDArray[np.float64]:
    """Magnitude of per-particle velocity vectors."""
    if velocities.ndim != 2 or velocities.shape[1] != 2:
        raise ValueError(
            f"velocities must be (N, 2), got shape {velocities.shape}"
        )
    return np.sqrt(np.sum(velocities**2, axis=1))


def phase_coexistence_fractions(
    rho_local: NDArray[np.float64],
    rho_star: float,
) -> tuple[float, float]:
    """Return (f_gas, f_liquid) under a given slowdown threshold rho_star.

    f_gas: fraction with rho < rho_star/2 (dilute phase, particles moving freely)
    f_liquid: fraction with rho > rho_star (dense phase, particles near-stalled)

    These are the two steady-state phases in Fily-Marchetti MIPS. Regimes
    that show either one without the other (pure flocking: all liquid; or
    dilute gas: all gas) are NOT phase-separated.
    """
    if rho_star <= 0:
        raise ValueError(f"rho_star must be positive, got {rho_star}")
    if rho_local.size == 0:
        return 0.0, 0.0
    f_gas = float((rho_local < rho_star / 2.0).mean())
    f_liquid = float((rho_local > rho_star).mean())
    return f_gas, f_liquid


def two_phase_coexistence_score(
    rho_local: NDArray[np.float64],
    rho_star: float,
) -> float:
    """Primary P2 metric: min(f_gas, f_liquid).

    Range [0, 0.5]. Value > 0.10 is strong two-phase coexistence; > 0.03
    is screening-level; == 0 means one-phase (flocking, uniform, or stuck).
    """
    f_gas, f_liq = phase_coexistence_fractions(rho_local, rho_star)
    return min(f_gas, f_liq)


def density_speed_anticorrelation(
    rho_local: NDArray[np.float64],
    speeds: NDArray[np.float64],
) -> float:
    """Pearson r(rho_local, |v|). Returns 0.0 if either is constant.

    Genuine MIPS: -0.95 < r < -0.3. Values r <= -0.99 usually indicate
    the dilute-Poisson artifact (few discrete rho values, few v values)
    rather than a real density-velocity mechanism — the detector treats
    |r| >= 0.99 as a rejection of DEFINITIVE tier but not of CONFIRMATION
    when the primary is strong.
    """
    if rho_local.size != speeds.size:
        raise ValueError(
            f"rho_local ({rho_local.size}) and speeds ({speeds.size}) "
            f"must have same length"
        )
    if rho_local.size < 3:
        return 0.0
    if rho_local.std() < 1e-12 or speeds.std() < 1e-12:
        return 0.0
    # Pearson correlation
    rho_c = rho_local - rho_local.mean()
    v_c = speeds - speeds.mean()
    num = float((rho_c * v_c).sum())
    den = float(np.sqrt((rho_c**2).sum() * (v_c**2).sum()))
    if den < 1e-12:
        return 0.0
    return num / den


def collect_density_and_speed_history(
    state_history: list[dict[str, Any]],
    box_size: Optional[float],
    r_cg: float = 1.0,
    periodic: bool = True,
    burn_in: int = 0,
    max_frames: Optional[int] = None,
) -> tuple[NDArray[np.float64], NDArray[np.float64], int]:
    """Flatten per-particle density and speed across measurement frames.

    Returns (all_rho, all_speeds, n_frames_used).

    If a snapshot already carries a ``local_density`` field (ABP does),
    it is used directly; otherwise we compute via particle_local_density.
    If a snapshot carries a ``speeds`` field, that is used; otherwise
    we compute magnitudes from ``velocities``.

    burn_in : int
        Skip the first `burn_in` snapshots.
    max_frames : int or None
        If set, subsample to at most this many evenly-spaced frames
        from the post-burn window (used to keep memory bounded on
        large histories).
    """
    if burn_in >= len(state_history):
        return np.zeros(0), np.zeros(0), 0
    frames = state_history[burn_in:]
    if max_frames is not None and len(frames) > max_frames:
        idx = np.linspace(0, len(frames) - 1, max_frames, dtype=int)
        frames = [frames[i] for i in idx]

    rho_parts = []
    v_parts = []
    use_box = box_size if periodic else None
    for snap in frames:
        if "local_density" in snap:
            rho = np.asarray(snap["local_density"], dtype=np.float64)
        else:
            pos = snap.get("positions")
            if pos is None:
                continue
            rho = particle_local_density(
                np.asarray(pos, dtype=np.float64),
                r_cg=r_cg,
                box_size=use_box,
            )
        if "speeds" in snap:
            spd = np.asarray(snap["speeds"], dtype=np.float64)
        else:
            vel = snap.get("velocities")
            if vel is None:
                continue
            spd = particle_speeds(np.asarray(vel, dtype=np.float64))
        if rho.shape != spd.shape:
            continue
        rho_parts.append(rho)
        v_parts.append(spd)
    if not rho_parts:
        return np.zeros(0), np.zeros(0), 0
    return (
        np.concatenate(rho_parts),
        np.concatenate(v_parts),
        len(rho_parts),
    )


def constant_speed_surrogate_null(
    positions_history: list[NDArray[np.float64]],
    rho_star: float,
    r_cg: float = 1.0,
    box_size: Optional[float] = None,
) -> float:
    """Surrogate null: what would primary be if all speeds were v0?

    This answers: is the observed two-phase structure a consequence of
    the POSITIONS alone (geometric clustering), independent of the
    density-speed feedback?

    We just recompute the two_phase score on the observed positions
    WITHOUT weighting by speed. Since the primary is purely
    position-based already, this returns the SAME score — which means
    the primary itself is mechanism-agnostic.

    This is why the detector needs the separate density_speed
    anticorrelation confirmation gate: without it, positional clustering
    from ANY mechanism (flocking, attraction) would trigger the primary.

    Returns the positional two-phase score averaged across frames.
    """
    if not positions_history:
        return 0.0
    scores = []
    use_box = box_size  # pass None to disable periodic
    for pos in positions_history:
        pos = np.asarray(pos, dtype=np.float64)
        if pos.ndim != 2 or pos.shape[1] != 2 or pos.shape[0] == 0:
            continue
        rho = particle_local_density(pos, r_cg=r_cg, box_size=use_box)
        scores.append(two_phase_coexistence_score(rho, rho_star))
    if not scores:
        return 0.0
    return float(np.mean(scores))


def mechanistic_null_test(
    primary_observed: float,
    model_metadata: dict[str, Any] | None,
) -> dict[str, Any]:
    """Check mechanistic-null conditions from model metadata.

    Returns a dict with keys:
      - null_rejects_mips: bool — True if metadata indicates the model
        HAS the ingredients for MIPS (density-dependent speed, no
        attraction, no alignment). If False, we cannot rule MIPS in
        (DEFINITIVE requires True).
      - has_density_dependent_speed: bool (from metadata, or None if absent)
      - has_attraction_rule: bool
      - has_alignment_rule: bool
      - metadata_used: bool
      - reason: str — explanatory text for null outcome

    This does NOT run a constant-speed simulation (that would require
    re-running the model, which the detector framework does not do).
    Instead it uses the model's self-reported rule flags. A model that
    advertises `has_attraction_rule=True` cannot produce a DEFINITIVE
    P2 even if the empirical primary is strong.
    """
    if model_metadata is None:
        return {
            "null_rejects_mips": False,
            "has_density_dependent_speed": None,
            "has_attraction_rule": None,
            "has_alignment_rule": None,
            "metadata_used": False,
            "reason": "no metadata available — cannot verify mechanism",
        }
    has_ddv = model_metadata.get("has_density_dependent_speed")
    has_attr = model_metadata.get("has_attraction_rule")
    has_align = model_metadata.get("has_alignment_rule")
    if has_ddv is None and has_attr is None and has_align is None:
        return {
            "null_rejects_mips": False,
            "has_density_dependent_speed": None,
            "has_attraction_rule": None,
            "has_alignment_rule": None,
            "metadata_used": False,
            "reason": "metadata present but no MIPS rule flags — "
                     "cannot verify mechanism",
        }
    # Explicit mechanism: density-dependent speed + no attraction + no alignment
    ok = (has_ddv is True) and (has_attr is False) and (has_align is False)
    if ok:
        reason = "mechanism confirmed: density-dependent speed + no attraction + no alignment"
    elif has_attr is True:
        reason = f"attraction_rule=True rejects MIPS mechanism"
    elif has_align is True:
        reason = f"alignment_rule=True rejects MIPS mechanism"
    elif has_ddv is False:
        reason = f"density-dependent speed absent — no MIPS mechanism"
    else:
        reason = "mechanism flags ambiguous"
    return {
        "null_rejects_mips": bool(ok),
        "has_density_dependent_speed": has_ddv,
        "has_attraction_rule": has_attr,
        "has_alignment_rule": has_align,
        "metadata_used": True,
        "reason": reason,
    }
