"""Equivalence test for the vectorised P26 driver (repro v2, group C spec P26 section 2).

Feed BistableDoubleWell.simulate() and double_well_fourier_driver() the same
pre-drawn normal increments for 4 trials x 2000 steps at D = 0.15 in the
reproduction frame (a = b = 1, A = 0.1, f = 0.1/2pi, dt = 0.01). The
trajectories must agree to <= 1e-12.
"""

import numpy as np

from epc.models.stochastic_resonance import (
    BistableDoubleWell,
    DoubleWellParams,
    double_well_fourier_driver,
)


class _ReplayRNG:
    """Stands in for np.random.Generator in simulate(): .normal() replays a
    fixed sequence (simulate() draws trial by trial, step by step)."""

    def __init__(self, seq):
        self._seq = list(seq)
        self._i = 0

    def normal(self):
        v = self._seq[self._i]
        self._i += 1
        return v


def test_driver_matches_simulate_trajectories():
    n_trials, n_steps, D = 4, 2000, 0.15
    params = DoubleWellParams(
        a=1.0, b=1.0, signal_amplitude=0.1,
        signal_frequency=0.1 / (2.0 * np.pi), dt=0.01,
        n_steps=n_steps, n_trials=n_trials, noise_levels=[D], seed=0,
    )
    z = np.random.default_rng(20261115).standard_normal((n_trials, n_steps))

    model = BistableDoubleWell(params)
    model.rng = _ReplayRNG(z.ravel())      # trial-major, as simulate() consumes
    hist = model.simulate()
    x_ref = np.array([h["x"] for h in hist]).reshape(n_trials, n_steps).T

    out = double_well_fourier_driver(
        params, D, n_trials=n_trials, n_burn=0, n_meas=n_steps,
        normals=z.T.copy(), chunk=777, return_trajectory=True,
    )
    diff = np.max(np.abs(out["trajectory"] - x_ref))
    assert diff <= 1e-12, diff
    # the trajectory must actually hop between wells for the test to bite
    assert np.ptp(x_ref) > 1.0


def test_driver_fourier_projection_matches_direct_sum():
    params = DoubleWellParams(
        a=1.0, b=1.0, signal_amplitude=0.1,
        signal_frequency=0.1 / (2.0 * np.pi), dt=0.01, seed=0,
    )
    rng = np.random.default_rng(5)
    out = double_well_fourier_driver(
        params, 0.15, n_trials=3, n_burn=500, n_meas=1500, omegas=[0.1, 0.09],
        rng=rng, chunk=256, return_trajectory=True,
    )
    t = np.arange(500, 2000) * params.dt
    seg = out["trajectory"][500:2000]
    for j, w in enumerate([0.1, 0.09]):
        direct = (seg * np.exp(-1j * w * t)[:, None]).sum(axis=0) * params.dt
        assert np.allclose(out["X"][j], direct, rtol=0, atol=1e-9)
