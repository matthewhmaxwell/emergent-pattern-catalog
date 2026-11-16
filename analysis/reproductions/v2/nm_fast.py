"""Generation update of epc.models.nowak_may.NowakMayModel, in C, without the state dict.

`generation(model)` computes the same next grid as NowakMayModel.step() (periodic
boundary): cooperator count in the 3x3 block incl. self, payoff count*1.0 for C and
count*b for D (float64), then the strategy of the highest payoff in the block, self
first and neighbours in the same scan order, strict '>'. It skips _state_dict() (no
grid copy, no Moran's I) and returns the number of cooperators n_C; the model's
coop_fraction equals n_C / (rows * cols) exactly. Bit-identity of the grids (and of
coop_fraction) against model.step() is checked by check_bit_identical.py.
"""
from __future__ import annotations

import ctypes

import numpy as np

import fastlib


def generation(model) -> int:
    if model.boundary != "periodic":
        raise ValueError("nm_fast supports boundary='periodic' only")
    R, C = model.rows, model.cols
    grid = np.ascontiguousarray(model.grid, dtype=np.int8)
    out = np.empty((R, C), dtype=np.int8)
    pay = getattr(model, "_repro_v2_pay", None)
    if pay is None or pay.shape != (R, C):
        pay = np.empty((R, C), dtype=np.float64)
        model._repro_v2_pay = pay
    n_coop = fastlib.lib().nm_generation(fastlib.ptr(grid, ctypes.c_int8), fastlib.ptr(out, ctypes.c_int8),
                                         fastlib.ptr(pay, ctypes.c_double), R, C, float(model.b))
    model.grid = out
    model._step_count += 1
    return int(n_coop)
