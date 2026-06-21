"""T2c OOD systems — held-out systems the detectors were NOT built on.

Four arms (the instrument's discrimination must hold on all):

  NULL        non-emergent systems -> must verdict NO-EMERGENCE.
  NOVELTY     genuinely emergent ABMs whose phenomenon is NOT among the 32
              -> must verdict EMERGENT-UNCLASSIFIED, and crucially NOT a false MATCH.
  RECOGNITION independent re-implementations of catalog phenomena (alternative
              models the detector was never tuned on) -> must MATCH the right pattern.
  PROBE       coverage-limit probes (clearly labeled, reported separately, NOT in
              headline rates): cases that expose known gaps of the generic
              emergence indicator (structure-without-growth; connectivity-only
              transitions the autocorrelation statistic cannot see).

Every builder takes a seed and returns (history, metadata). History is a list of
observation-bundle frames using the canonical keys the detectors/emergence read
(`grid`/`field` 2-D arrays, `positions` Nx2, `phases`, scalar vectors).

The NULL/NOVELTY/PROBE arms are self-contained minimal ABMs implemented here so
their format and ground truth are unambiguous. The RECOGNITION arm imports the
repo's alternative-implementation models (BooleanGRN, FractionThreshold, El Farol,
Integral controller, MultiBasinGRN, ThresholdUnit) — the same 7/7 cross-model set.
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple

import numpy as np

Frame = Dict[str, Any]
History = List[Frame]
Built = Tuple[History, Optional[Dict[str, Any]]]


# ----------------------------------------------------------------------------
# NULL arm — must read NO-EMERGENCE
# ----------------------------------------------------------------------------
def null_spatial_noise(seed: int = 0, N: int = 200, L: float = 30.0,
                       n_frames: int = 60) -> Built:
    """Independent uniform points re-drawn every frame: no structure, no growth."""
    rng = np.random.default_rng(seed)
    return [{"positions": rng.uniform(0, L, size=(N, 2))} for _ in range(n_frames)], None


def null_random_walk(seed: int = 0, N: int = 200, L: float = 60.0,
                     n_frames: int = 80, step: float = 0.8) -> Built:
    """Independent Brownian walkers: diffuse apart, never cluster."""
    rng = np.random.default_rng(seed)
    pos = rng.uniform(0, L, size=(N, 2))
    frames = [{"positions": pos.copy()}]
    for _ in range(n_frames - 1):
        pos = pos + rng.normal(0, step, size=(N, 2))
        frames.append({"positions": pos.copy()})
    return frames, None


def null_uncoupled_phases(seed: int = 0, N: int = 300, n_frames: int = 120,
                          dt: float = 0.05) -> Built:
    """Oscillators with independent natural frequencies and ZERO coupling:
    the Kuramoto order parameter stays at its finite-N floor — no synchrony."""
    rng = np.random.default_rng(seed)
    omega = rng.normal(0.0, 1.0, size=N)
    theta = rng.uniform(-np.pi, np.pi, size=N)
    frames = []
    for _ in range(n_frames):
        theta = (theta + omega * dt) % (2 * np.pi)
        frames.append({"phases": theta.copy()})
    return frames, None


def null_frozen_noise(seed: int = 0, N: int = 200, L: float = 30.0,
                      n_frames: int = 60) -> Built:
    """A single unstructured (dispersed) configuration, held fixed: low structure,
    zero order gain."""
    rng = np.random.default_rng(seed)
    pos = rng.uniform(0, L, size=(N, 2))
    return [{"positions": pos.copy()} for _ in range(n_frames)], None


# ----------------------------------------------------------------------------
# NOVELTY arm — emergent, out-of-catalog -> must read EMERGENT-UNCLASSIFIED
# ----------------------------------------------------------------------------
def nov_dla(seed: int = 0, L: int = 61, n_particles: int = 360,
            n_frames: int = 80) -> Built:
    """Diffusion-limited aggregation: a single fractal dendrite grows by sticking
    random walkers. Emergent branching structure; NOT any catalog pattern
    (P1 is type-segregation, not fractal growth)."""
    rng = np.random.default_rng(seed)
    grid = np.zeros((L, L), dtype=np.float32)
    c = L // 2
    grid[c, c] = 1.0
    maxr = 2
    frames: History = []
    snap = max(1, n_particles // n_frames)

    def stuck(x: int, y: int) -> bool:
        for dx, dy in ((1, 0), (-1, 0), (0, 1), (0, -1)):
            nx, ny = x + dx, y + dy
            if 0 <= nx < L and 0 <= ny < L and grid[nx, ny] > 0:
                return True
        return False

    def launch() -> Tuple[int, int]:
        rad = min(maxr + 4, L // 2 - 1)
        ang = rng.uniform(0, 2 * np.pi)
        x = int(c + rad * np.cos(ang)); y = int(c + rad * np.sin(ang))
        return min(max(x, 1), L - 2), min(max(y, 1), L - 2)

    for p in range(n_particles):
        x, y = launch()
        for _ in range(8 * L * L):
            if stuck(x, y):
                grid[x, y] = 1.0
                maxr = max(maxr, int(np.hypot(x - c, y - c)))
                break
            d = rng.integers(0, 4)
            x += (1 if d == 0 else -1 if d == 1 else 0)
            y += (1 if d == 2 else -1 if d == 3 else 0)
            if (np.hypot(x - c, y - c) > maxr + 12 or
                    x <= 0 or y <= 0 or x >= L - 1 or y >= L - 1):
                x, y = launch()
        if p % snap == 0:
            frames.append({"grid": grid.copy()})
    frames.append({"grid": grid.copy()})
    return frames, None


def nov_keller_segel(seed: int = 0, N: int = 200, L: float = 50.0, G: int = 64,
                     n_steps: int = 400, n_frames: int = 80) -> Built:
    """Keller-Segel chemotactic aggregation: agents secrete a diffusing
    chemoattractant and climb its SELF-GENERATED gradient, collapsing into
    dense clumps. Mechanistically distinct from P17 (external fixed gradient)
    and P1 (type segregation)."""
    rng = np.random.default_rng(seed)
    pos = rng.uniform(0, L, size=(N, 2))
    chem = np.zeros((G, G), dtype=np.float64)
    frames: History = []
    snap = max(1, n_steps // n_frames)

    def cell(p: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        gx = np.clip((p[:, 0] / L * G).astype(int), 0, G - 1)
        gy = np.clip((p[:, 1] / L * G).astype(int), 0, G - 1)
        return gx, gy

    for t in range(n_steps):
        gx, gy = cell(pos)
        np.add.at(chem, (gx, gy), 1.0)                     # deposit
        lap = (np.roll(chem, 1, 0) + np.roll(chem, -1, 0) +
               np.roll(chem, 1, 1) + np.roll(chem, -1, 1) - 4 * chem)
        chem = (chem + 0.18 * lap) * 0.95                  # diffuse + decay
        gxp = np.clip(gx + 1, 0, G - 1); gxm = np.clip(gx - 1, 0, G - 1)
        gyp = np.clip(gy + 1, 0, G - 1); gym = np.clip(gy - 1, 0, G - 1)
        grad = np.stack([chem[gxp, gy] - chem[gxm, gy],
                         chem[gx, gyp] - chem[gx, gym]], axis=1)
        nrm = np.linalg.norm(grad, axis=1, keepdims=True) + 1e-9
        pos = (pos + 0.6 * grad / nrm + rng.normal(0, 0.25, size=(N, 2))) % L
        if t % snap == 0:
            frames.append({"positions": pos.copy()})
    frames.append({"positions": pos.copy()})
    return frames, None


def nov_active_nematic(seed: int = 0, N: int = 160, L: float = 18.0,
                       n_steps: int = 300, n_frames: int = 80,
                       v: float = 0.06, R: float = 1.2, eta: float = 0.12) -> Built:
    """Self-propelled rods with NEMATIC (head-tail-symmetric) alignment: local
    averaging of 2*orientation. Builds nematic order with no global POLAR
    direction — so the polar flocking order parameter (P5) stays low. A direct
    test of whether the flocking detector over-claims on apolar alignment."""
    rng = np.random.default_rng(seed)
    pos = rng.uniform(0, L, size=(N, 2))
    ang = rng.uniform(0, np.pi, size=N)                    # orientation in [0, pi)
    frames: History = []
    snap = max(1, n_steps // n_frames)
    for t in range(n_steps):
        d2 = ((pos[:, None, :] - pos[None, :, :]) ** 2).sum(-1)
        A = (d2 < R * R).astype(float)
        local = A @ np.exp(2j * ang)                       # nematic mean field
        ang = (np.angle(local) / 2.0 + rng.normal(0, eta, size=N)) % np.pi
        vel = np.stack([np.cos(ang), np.sin(ang)], axis=1) * v * L
        pos = (pos + vel) % L
        if t % snap == 0:
            frames.append({"positions": pos.copy(), "velocities": vel.copy()})
    frames.append({"positions": pos.copy(), "velocities": vel.copy()})
    return frames, None


def nov_phase_ordering(seed: int = 0, G: int = 64, n_steps: int = 2000,
                       n_frames: int = 80, dt: float = 0.1, kappa: float = 1.0) -> Built:
    """Allen-Cahn phase ordering: an initially noisy field coarsens into growing
    +/-1 domains with a smoothly autocorrelated, isotropic morphology. Tests
    whether the Turing detector (P3) over-claims on coarsening domains that are
    NOT a stationary reaction-diffusion wavelength."""
    rng = np.random.default_rng(seed)
    u = rng.uniform(-0.05, 0.05, size=(G, G))
    frames: History = []
    snap = max(1, n_steps // n_frames)

    def lap(f: np.ndarray) -> np.ndarray:
        return (np.roll(f, 1, 0) + np.roll(f, -1, 0) +
                np.roll(f, 1, 1) + np.roll(f, -1, 1) - 4 * f)

    for t in range(n_steps):
        u = u + dt * (-(u ** 3 - u) + kappa * lap(u))
        if t % snap == 0:
            frames.append({"field": u.copy()})
    frames.append({"field": u.copy()})
    return frames, None


# ----------------------------------------------------------------------------
# PROBE arm — coverage-limit demonstrations (reported separately, not scored)
# ----------------------------------------------------------------------------
def probe_static_blob(seed: int = 0, N: int = 200, L: float = 30.0,
                      sigma: float = 2.5, n_frames: int = 60) -> Built:
    """A pre-formed tight cluster held fixed. High spatial structure but ZERO
    order gain: it never self-organized. Exposes that the emergence indicator's
    structure-vs-shuffle term can fire without growth (structure-without-history)."""
    rng = np.random.default_rng(seed)
    pos = np.clip(rng.normal(L / 2, sigma, size=(N, 2)), 0, L)
    return [{"positions": pos.copy()} for _ in range(n_frames)], None


def probe_percolation(seed: int = 0, G: int = 64, n_frames: int = 80,
                      frac: float = 0.62) -> Built:
    """Random site occupation crossing the percolation threshold. The emergent
    object is a SPANNING CLUSTER (a connectivity transition), but the field has
    ~zero spatial autocorrelation, so a Moran-based emergence indicator cannot
    see it. Exposes a morphology coverage gap (connectivity != autocorrelation)."""
    rng = np.random.default_rng(seed)
    order = rng.permutation(G * G)
    grid = np.zeros(G * G, dtype=np.float32)
    n_final = int(frac * G * G)
    snap = max(1, n_final // n_frames)
    frames: History = []
    for k in range(n_final):
        grid[order[k]] = 1.0
        if k % snap == 0:
            frames.append({"grid": grid.reshape(G, G).copy()})
    frames.append({"grid": grid.reshape(G, G).copy()})
    return frames, None


# ----------------------------------------------------------------------------
# RECOGNITION arm — independent re-implementations of catalog phenomena
# ----------------------------------------------------------------------------
def _meta(m: Any) -> Optional[Dict[str, Any]]:
    try:
        return m.get_metadata()
    except Exception:
        return None


def recog_p16_boolean_grn(seed: int = 0) -> Built:
    """P16 associative memory via a Boolean gene-regulatory network (not Hopfield)."""
    from epc.models.hopfield import BooleanGRN, BooleanGRNParams
    m = BooleanGRN(BooleanGRNParams(N=20, P=4, corruption=0.2, seed=seed))
    return m.simulate(n_cues=4), _meta(m)


def recog_p20_fraction_threshold(seed: int = 0) -> Built:
    """P20 quorum sensing via a fraction-threshold switch (not autoinducer ODE)."""
    from epc.models.quorum_sensing import FractionThresholdModel, FractionThresholdParams
    m = FractionThresholdModel(FractionThresholdParams(
        n_agents=200, n_density_levels=24, n_steps_per_density=120, seed=seed))
    return m.simulate(), _meta(m)


def recog_p23_el_farol(seed: int = 0) -> Built:
    """P23 anti-coordination via the El Farol Bar (not the Minority Game)."""
    from epc.models.minority_game import ElFarolBar, ElFarolParams
    m = ElFarolBar(ElFarolParams(n_agents=100, capacity=60, n_rounds=500, seed=seed))
    return m.simulate(), _meta(m)


def recog_p24_integral(seed: int = 0) -> Built:
    """P24 homeostasis via an INTEGRAL controller (detector was tuned on the
    proportional one). Known generalization limit: ranks P24 #1 but the firing
    threshold is calibrated to the proportional transient."""
    from epc.models.homeostasis import IntegralHomeostat, HomeostatParams, PerturbationSchedule
    m = IntegralHomeostat(HomeostatParams(setpoint=10.0, gain=1.0, dt=0.1,
                                          noise_std=0.5, seed=seed))
    hist = m.simulate(n_steps=1000,
                      schedule=PerturbationSchedule(onset=30.0, offset=70.0, amplitude=5.0))
    return hist, _meta(m)


def recog_p25_multibasin(seed: int = 0) -> Built:
    """P25 equifinality via a multi-basin GRN (not the canalized landscape)."""
    from epc.models.canalization import MultiBasinGRN
    m = MultiBasinGRN(n_genes=10, n_ics=20, n_steps=200, seed=seed)
    return m.simulate(), _meta(m)


def recog_p26_threshold_unit(seed: int = 0) -> Built:
    """P26 stochastic resonance via a threshold unit (not the double well)."""
    from epc.models.stochastic_resonance import ThresholdUnit, ThresholdUnitParams
    m = ThresholdUnit(ThresholdUnitParams(
        threshold=1.0, signal_amplitude=0.7, signal_frequency=0.05,
        dt=0.01, n_steps=4000,
        noise_levels=[0.0, 0.2, 0.5, 0.7, 1.0, 1.5, 3.0], seed=seed))
    return m.simulate(), _meta(m)


# ----------------------------------------------------------------------------
# Registry
# ----------------------------------------------------------------------------
# Each entry: name, cls (null|novelty|recognition|probe), expect (pattern id or None),
# build (seed -> Built), stochastic (vary seed across replicates).
SYSTEMS: List[Dict[str, Any]] = [
    # NULL
    {"name": "null_spatial_noise",   "cls": "null", "expect": None, "build": null_spatial_noise,   "stochastic": True},
    {"name": "null_random_walk",     "cls": "null", "expect": None, "build": null_random_walk,     "stochastic": True},
    {"name": "null_uncoupled_phases","cls": "null", "expect": None, "build": null_uncoupled_phases,"stochastic": True},
    {"name": "null_frozen_noise",    "cls": "null", "expect": None, "build": null_frozen_noise,    "stochastic": True},
    # NOVELTY
    {"name": "nov_dla",              "cls": "novelty", "expect": None, "build": nov_dla,            "stochastic": True},
    {"name": "nov_keller_segel",     "cls": "novelty", "expect": None, "build": nov_keller_segel,   "stochastic": True},
    {"name": "nov_active_nematic",   "cls": "novelty", "expect": None, "build": nov_active_nematic, "stochastic": True},
    {"name": "nov_phase_ordering",   "cls": "novelty", "expect": None, "build": nov_phase_ordering, "stochastic": True},
    # RECOGNITION (alt implementations)
    {"name": "recog_p16_boolean_grn",      "cls": "recognition", "expect": "P16", "build": recog_p16_boolean_grn,      "stochastic": True},
    {"name": "recog_p20_fraction_threshold","cls": "recognition","expect": "P20", "build": recog_p20_fraction_threshold,"stochastic": True},
    {"name": "recog_p23_el_farol",         "cls": "recognition", "expect": "P23", "build": recog_p23_el_farol,         "stochastic": True},
    {"name": "recog_p24_integral",         "cls": "recognition", "expect": "P24", "build": recog_p24_integral,         "stochastic": True},
    {"name": "recog_p25_multibasin",       "cls": "recognition", "expect": "P25", "build": recog_p25_multibasin,       "stochastic": True},
    {"name": "recog_p26_threshold_unit",   "cls": "recognition", "expect": "P26", "build": recog_p26_threshold_unit,   "stochastic": True},
    # PROBE (coverage-limit, reported separately)
    {"name": "probe_static_blob",   "cls": "probe", "expect": None, "build": probe_static_blob,   "stochastic": True},
    {"name": "probe_percolation",   "cls": "probe", "expect": None, "build": probe_percolation,   "stochastic": True},
]
