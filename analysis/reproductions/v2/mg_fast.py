"""Faster inner loop for epc.models.minority_game.MinorityGame (S = 2).

`advance(model, n)` plays n rounds and leaves the model in exactly the state that
`model.params.n_rounds = n; model.simulate()` would leave it in: same scores, same
history_bits, same rng state. It returns only the per-round attendance n1, winning
side and history index mu, instead of a list of dicts with a copy of every choice.

Implementation differences from simulate(), none of which changes a result:
  - the round loop runs in C (fastloops.c: mg_rounds);
  - mu is updated incrementally, mu <- ((mu << 1) | w) mod P, which is what
    np.roll(history_bits, -1) followed by _history_index() computes;
  - with tie_break='random', simulate() calls rng.integers(0, 2, N) once per round.
    Generator.integers(0, 2) takes one 32-bit value per coin from the bit generator,
    which keeps the unused half of a 64-bit draw in its own buffer between calls
    (has_uint32). The coins of n successive calls of size N are therefore one
    contiguous stream, and a single call of size n*N, reshaped to (n, N), yields the
    same coins and leaves the same rng state (probed on numpy 2.5.0; verified below).
Bit-identity against simulate() (attendance, winning side, scores, history bits and
rng state, all four option combinations) is checked by check_bit_identical.py.
"""
from __future__ import annotations

import numpy as np

import fastlib


def _tables(model):
    cache = getattr(model, "_repro_v2_strat_t", None)
    if cache is None:
        if model.S != 2:
            raise ValueError("mg_fast supports n_strategies = 2 only")
        cache = np.ascontiguousarray(model.strategies.transpose(2, 0, 1)).astype(np.int8)  # (P, N, 2)
        model._repro_v2_strat_t = cache
    return cache


def history_index(bits) -> int:
    idx = 0
    for b in bits:
        idx = 2 * idx + int(b)
    return idx


def _coins(rng, n_rounds: int, N: int) -> np.ndarray:
    c = rng.integers(0, 2, size=n_rounds * N).reshape(n_rounds, N)
    return np.ascontiguousarray(c, dtype=np.int64)


def advance(model, n_rounds: int):
    """Play n_rounds rounds; return (attendance int64[n], winning int8[n], mu int64[n])."""
    N, P, m = model.N, model.P, model.m
    strat_t = _tables(model)
    tie_random = model.params.tie_break == "random"
    linear = model.params.payoff == "linear"
    if model.scores.dtype != np.float64 or not model.scores.flags["C_CONTIGUOUS"]:
        raise ValueError("scores must be C-contiguous float64")
    coins = _coins(model.rng, n_rounds, N) if tie_random else np.zeros((1, 1), dtype=np.int64)
    mu_io = np.array([history_index(model.history_bits)], dtype=np.int64)
    att = np.empty(n_rounds, dtype=np.int64)
    win = np.empty(n_rounds, dtype=np.int8)
    mus = np.empty(n_rounds, dtype=np.int64)
    fastlib.lib().mg_rounds(
        fastlib.ptr(model.scores, fastlib.ctypes.c_double), fastlib.ptr(strat_t, fastlib.ctypes.c_int8),
        N, P, fastlib.ptr(coins, fastlib.ctypes.c_int64), int(tie_random), int(linear), n_rounds,
        fastlib.ptr(mu_io, fastlib.ctypes.c_int64), fastlib.ptr(att, fastlib.ctypes.c_int64),
        fastlib.ptr(win, fastlib.ctypes.c_int8), fastlib.ptr(mus, fastlib.ctypes.c_int64))
    mu = int(mu_io[0])
    # write back history bits (most significant bit first, as _history_index reads them)
    model.history_bits = np.array([(mu >> (m - 1 - k)) & 1 for k in range(m)], dtype=np.int8)
    return att, win, mus
