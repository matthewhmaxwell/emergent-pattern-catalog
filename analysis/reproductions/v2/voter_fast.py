"""Faster sweep for epc.models.voter.VoterModel: same RNG calls, C update loop.

`sweep(model)` advances the model by one Monte Carlo sweep and leaves grid, rng state
and step count exactly as `model.step()` would, without building the state dict (no
grid copy, no Moran's I). The three rng.integers calls are made in Python exactly as
step() makes them; only the sequential copy loop runs in C (fastloops.c: voter_sweep),
in the same order. Bit-identity is checked by check_bit_identical.py.
"""
from __future__ import annotations

import ctypes

import numpy as np

import fastlib


def sweep(model) -> None:
    R, C = model.rows, model.cols
    N = R * C
    g = model.grid
    if g.dtype != np.int8 or not g.flags["C_CONTIGUOUS"]:
        raise ValueError("grid must be C-contiguous int8")
    offs = getattr(model, "_repro_v2_offs", None)
    if offs is None:
        dr = np.array([o[0] for o in model._n_offsets], dtype=np.int64)
        dc = np.array([o[1] for o in model._n_offsets], dtype=np.int64)
        offs = (dr, dc)
        model._repro_v2_offs = offs
    # identical RNG calls, in the same order, as VoterModel.step()
    site_rs = model.rng.integers(0, R, size=N)
    site_cs = model.rng.integers(0, C, size=N)
    nb_idx = model.rng.integers(0, model._n_neighbors, size=N)
    for arr in (site_rs, site_cs, nb_idx):
        if arr.dtype != np.int64:
            raise ValueError("unexpected rng output dtype")
    fastlib.lib().voter_sweep(fastlib.ptr(g, ctypes.c_int8), R, C,
                              fastlib.ptr(site_rs, ctypes.c_int64), fastlib.ptr(site_cs, ctypes.c_int64),
                              fastlib.ptr(nb_idx, ctypes.c_int64), N,
                              fastlib.ptr(offs[0], ctypes.c_int64), fastlib.ptr(offs[1], ctypes.c_int64))
    model._step_count += 1
