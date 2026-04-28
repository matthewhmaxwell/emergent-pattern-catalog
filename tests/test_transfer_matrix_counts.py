"""Transfer-matrix count verification (Sprint 22 D5).

CI-style test that pins the canonical §5.1 counts at both endpoints:
the `scripts/count_transfer_matrix.py` script's `compute_counts()`
function (the source of truth) AND the existing assertions in
`test_orchestration.py` (which are what the registry exposes through
the orchestration API).

Future drift in either endpoint produces a test failure here,
preventing a silent divergence between paper figures, registry, and
script output. The intended response when this test fails is:

  1. Decide whether the registry change was intentional.
  2. If yes, update the EXPECTED constants below to the new measured
     values (run `python scripts/count_transfer_matrix.py` to get
     them) AND update §5.1's prose to match.
  3. If no, revert the registry change.

Sprint 22 D5: addresses Sprint 21 carry-forward #19. The pinned values
were measured against Sprint 21 D1 HEAD (commit 0e267c5).
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

from epc.orchestration import (
    DETECTOR_REGISTRY,
    MODEL_REGISTRY,
    get_compatible_pairs,
)


# Load the count script as a module so the test can call its functions
# directly. The script lives in scripts/, which isn't a package, so we
# load it by file path.
_SCRIPT_PATH = (
    Path(__file__).parent.parent / 'scripts' / 'count_transfer_matrix.py'
)
_spec = importlib.util.spec_from_file_location(
    'count_transfer_matrix', _SCRIPT_PATH,
)
_module = importlib.util.module_from_spec(_spec)
assert _spec.loader is not None
_spec.loader.exec_module(_module)
compute_counts = _module.compute_counts
DISPLAY_FOLDS = _module.DISPLAY_FOLDS


# Expected values, measured against Sprint 21 D1 HEAD. If the registry
# changes, run `python scripts/count_transfer_matrix.py` to get the new
# values, update this dict, and update §5.1 in the paper draft.
EXPECTED: dict[str, int] = {
    'n_models': 20,
    'n_detectors': 19,
    'n_total_cells': 380,
    'n_compatible': 79,
    'n_substrate_mismatch': 274,
    'n_missing_observable': 27,
    'n_total_rejections': 301,
    'n_displayed_rows': 19,
    'n_displayed_cells': 361,
    'n_displayed_compatible': 77,
    'n_displayed_rejections': 284,
}


class TestComputeCountsAgainstExpected:
    """Pin every count emitted by `compute_counts()` to its expected
    value. These are the canonical figures that §5.1 reports."""

    @pytest.fixture(scope='class')
    def counts(self) -> dict:
        return compute_counts()

    @pytest.mark.parametrize('key', sorted(EXPECTED.keys()))
    def test_count_matches_expected(self, counts, key):
        actual = counts[key]
        expected = EXPECTED[key]
        assert actual == expected, (
            f"compute_counts()['{key}'] = {actual}, expected {expected}. "
            f"If this is an intentional registry change, run "
            f"`python scripts/count_transfer_matrix.py`, update EXPECTED "
            f"in this test, and update §5.1's prose to match."
        )


class TestComputeCountsCrossChecksOrchestrationApi:
    """Cross-check that the script's counts are consistent with what the
    orchestration API exposes directly. If either endpoint drifts, both
    fail — and the failure mode tells us which one drifted.

    These tests deliberately do NOT call `compute_counts` and then
    compare against the registry; they each fetch the value from a
    different endpoint and assert equality. That gives independent
    confirmation that the script and the public API agree on the same
    underlying registry."""

    def test_n_models_matches_registry(self):
        counts = compute_counts()
        assert counts['n_models'] == len(MODEL_REGISTRY)

    def test_n_detectors_matches_registry(self):
        counts = compute_counts()
        assert counts['n_detectors'] == len(DETECTOR_REGISTRY)

    def test_n_compatible_matches_get_compatible_pairs(self):
        counts = compute_counts()
        assert counts['n_compatible'] == len(get_compatible_pairs())

    def test_n_total_cells_matches_product(self):
        counts = compute_counts()
        assert (
            counts['n_total_cells']
            == len(MODEL_REGISTRY) * len(DETECTOR_REGISTRY)
        )


class TestComputeCountsInternalConsistency:
    """Internal sanity checks: the count breakdown must sum cleanly at
    both registry and display level, and the fold accounting must be
    coherent."""

    def test_registry_breakdown_sums_to_total(self):
        counts = compute_counts()
        assert (
            counts['n_compatible']
            + counts['n_substrate_mismatch']
            + counts['n_missing_observable']
            == counts['n_total_cells']
        )

    def test_total_rejections_matches_breakdown(self):
        counts = compute_counts()
        assert (
            counts['n_total_rejections']
            == counts['n_substrate_mismatch'] + counts['n_missing_observable']
        )

    def test_displayed_breakdown_sums_to_displayed_total(self):
        counts = compute_counts()
        assert (
            counts['n_displayed_compatible']
            + counts['n_displayed_rejections']
            == counts['n_displayed_cells']
        )

    def test_display_fold_drops_correct_row_count(self):
        counts = compute_counts()
        assert (
            counts['n_displayed_rows']
            == counts['n_models'] - len(DISPLAY_FOLDS)
        )

    def test_display_fold_drops_correct_cell_count(self):
        counts = compute_counts()
        n_dropped = (
            counts['n_total_cells'] - counts['n_displayed_cells']
        )
        assert n_dropped == len(DISPLAY_FOLDS) * counts['n_detectors']

    def test_display_folds_target_real_models(self):
        for kept, dropped in DISPLAY_FOLDS:
            assert kept in MODEL_REGISTRY, (
                f"DISPLAY_FOLDS references unknown 'kept' model: {kept}"
            )
            assert dropped in MODEL_REGISTRY, (
                f"DISPLAY_FOLDS references unknown 'dropped' model: {dropped}"
            )

    def test_display_folds_pair_compatibility_identical(self):
        """Folded variants must have identical compatibility profiles —
        that's what makes folding them in the displayed table sound."""
        from epc.orchestration import check_compatibility
        for kept, dropped in DISPLAY_FOLDS:
            for d in DETECTOR_REGISTRY:
                r_kept = check_compatibility(kept, d)
                r_dropped = check_compatibility(dropped, d)
                assert r_kept.compatible == r_dropped.compatible, (
                    f"Fold violation: {kept} × {d} compatible="
                    f"{r_kept.compatible} but {dropped} × {d} compatible="
                    f"{r_dropped.compatible}. The two variants no longer "
                    f"share a compatibility profile and cannot be folded "
                    f"into a single display row."
                )
