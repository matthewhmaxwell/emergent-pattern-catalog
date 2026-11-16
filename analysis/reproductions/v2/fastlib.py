"""Compile fastloops.c with gcc into a temporary directory and load it with ctypes.

Nothing is written into the repository. Flags: -O2 -fPIC -shared (no -ffast-math, so
the float64 arithmetic is IEEE and matches numpy).
"""
from __future__ import annotations

import ctypes
import subprocess
import tempfile
from pathlib import Path

_HERE = Path(__file__).resolve().parent
_LIB = None

i8p = ctypes.POINTER(ctypes.c_int8)
i64p = ctypes.POINTER(ctypes.c_int64)
f64p = ctypes.POINTER(ctypes.c_double)
CFLAGS = ["-O2", "-fPIC", "-shared", "-fno-fast-math"]


def lib():
    global _LIB
    if _LIB is None:
        build = Path(tempfile.mkdtemp(prefix="repro_v2_fastloops_"))
        so = build / "fastloops.so"
        subprocess.run(["gcc", *CFLAGS, "-o", str(so), str(_HERE / "fastloops.c")], check=True)
        L = ctypes.CDLL(str(so))
        L.mg_rounds.argtypes = [f64p, i8p, ctypes.c_int64, ctypes.c_int64, i64p, ctypes.c_int,
                                ctypes.c_int, ctypes.c_int64, i64p, i64p, i8p, i64p]
        L.mg_rounds.restype = None
        L.nm_generation.argtypes = [i8p, i8p, f64p, ctypes.c_int64, ctypes.c_int64, ctypes.c_double]
        L.nm_generation.restype = ctypes.c_int64
        L.voter_sweep.argtypes = [i8p, ctypes.c_int64, ctypes.c_int64, i64p, i64p, i64p,
                                  ctypes.c_int64, i64p, i64p]
        L.voter_sweep.restype = None
        _LIB = L
    return _LIB


def ptr(a, ctype):
    if not a.flags["C_CONTIGUOUS"]:
        raise ValueError("array must be C-contiguous")
    return a.ctypes.data_as(ctypes.POINTER(ctype))
