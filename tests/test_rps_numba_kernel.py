"""Bit-identity test for the Numba inner loop of RPSSpatialModel (repro v2, group C spec P12 section 2).

3 seeds x 50 generations at L = 30, M in {1e-5, 1e-3}: the grid (and the
per-generation reaction counts) of kernel='numba' must equal those of the
reference kernel='python' at every generation.
"""

import numpy as np
import pytest

from epc.models.rps_spatial import RPSSpatialModel

numba = pytest.importorskip("numba")


@pytest.mark.parametrize("mobility", [1e-5, 1e-3])
@pytest.mark.parametrize("seed", [1, 2, 3])
def test_numba_kernel_bit_identical(seed, mobility):
    ref = RPSSpatialModel(rows=30, cols=30, mobility=mobility, seed=seed)
    fast = RPSSpatialModel(rows=30, cols=30, mobility=mobility, seed=seed, kernel="numba")
    s_ref = ref.setup()
    s_fast = fast.setup()
    assert np.array_equal(s_ref["grid"], s_fast["grid"])
    totals = {"n_selections": 0, "n_reproductions": 0, "n_exchanges": 0}
    for g in range(50):
        a = ref.step()
        b = fast.step()
        assert np.array_equal(a["grid"], b["grid"]), f"grid differs at generation {g + 1}"
        for key in ("n_selections", "n_reproductions", "n_exchanges"):
            assert a[key] == b[key], (key, g + 1)
            totals[key] += a[key]
    # the dynamics must have done something for the test to bite
    assert all(v > 0 for v in totals.values()), totals


def test_default_kernel_is_python():
    assert RPSSpatialModel(rows=10, cols=10, mobility=1e-4).kernel == "python"
