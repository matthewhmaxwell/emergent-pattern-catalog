"""Programmatic transfer-matrix count verification (Sprint 22 D5).

Walks `MODEL_REGISTRY` and `DETECTOR_REGISTRY` and emits the canonical
counts that §5.1 of the paper reports as the transfer-matrix shape:

  - n_models                      registry-level model count
  - n_detectors                   registry-level detector count
  - n_total_cells                 n_models × n_detectors
  - n_compatible                  pairs that pass both substrate and
                                  observable prerequisites
  - n_substrate_mismatch          pairs rejected at the substrate gate
  - n_missing_observable          pairs rejected at the observable gate
                                  (after passing the substrate gate)
  - n_total_rejections            n_substrate_mismatch + n_missing_observable
  - n_displayed_rows              models after collapsing equivalent
                                  Zhang variants (zhang_sequential and
                                  zhang_threaded share substrate,
                                  observables, and primary patterns, so
                                  §5.1 displays them as a single row)
  - n_displayed_cells             n_displayed_rows × n_detectors
  - n_displayed_compatible        compatible cells in the displayed table
                                  (drops one Zhang variant's compatible
                                  cells from the registry-level count)
  - n_displayed_rejections        rejections in the displayed table

This script is the single source of truth for §5.1's figures.
`tests/test_transfer_matrix_counts.py` pins these emit values; if the
registry changes, the test fails first, and §5.1 must be updated to
match the script's output (not the other way around).

Usage:
  python scripts/count_transfer_matrix.py
  python scripts/count_transfer_matrix.py --json
  python -c "from scripts.count_transfer_matrix import compute_counts; print(compute_counts())"

Sprint 22 D5: addresses Sprint 21 carry-forward #19. Look-before-touching
audit at this sprint's clone time found that §5.1's prose claims of
"282 correctly eliminated" and "79 displayed audited cells" are mutually
inconsistent with the registry's actual structure (the correct displayed-
table figures are 284 rejections and 77 compatible). §5.1 prose
correction is a separate, paper-side carry-forward; this script pins
the empirical truth.
"""
from __future__ import annotations

import argparse
import json
import sys
from typing import Any

from epc.orchestration import (
    DETECTOR_REGISTRY,
    MODEL_REGISTRY,
    check_compatibility,
)


# Models that are folded into a single display row in §5.1 because they
# share substrate, observables, and primary patterns. Listed as
# (variant_kept, variant_dropped) — the dropped variant is removed from
# the displayed table while the registry still differentiates them.
DISPLAY_FOLDS: list[tuple[str, str]] = [
    ('zhang_sequential', 'zhang_threaded'),
]


def compute_counts() -> dict[str, Any]:
    """Compute the full set of transfer-matrix counts.

    Returns a dict suitable for JSON serialization. All values are ints
    derived directly from MODEL_REGISTRY / DETECTOR_REGISTRY, with no
    hard-coded numbers.
    """
    n_models = len(MODEL_REGISTRY)
    n_detectors = len(DETECTOR_REGISTRY)
    n_total_cells = n_models * n_detectors

    n_compatible = 0
    n_substrate_mismatch = 0
    n_missing_observable = 0
    for m in MODEL_REGISTRY:
        for d in DETECTOR_REGISTRY:
            r = check_compatibility(m, d)
            if r.compatible:
                n_compatible += 1
            elif r.reason.startswith('substrate_mismatch'):
                n_substrate_mismatch += 1
            elif r.reason.startswith('missing_observable'):
                n_missing_observable += 1
            else:
                # Defensive: any future rejection reason should be
                # surfaced rather than silently miscounted.
                raise ValueError(
                    f'Unrecognized rejection reason for ({m}, {d}): '
                    f'{r.reason!r}. Update DISPLAY_FOLDS or rejection '
                    f'classification in count_transfer_matrix.py.'
                )

    n_total_rejections = n_substrate_mismatch + n_missing_observable
    assert n_compatible + n_total_rejections == n_total_cells, (
        f'Cell-count sanity check failed: compatible ({n_compatible}) + '
        f'rejections ({n_total_rejections}) ≠ total ({n_total_cells})'
    )

    # Displayed-table counts: drop one variant per fold from the rows.
    n_dropped_rows = len(DISPLAY_FOLDS)
    dropped_models = {dropped for _, dropped in DISPLAY_FOLDS}
    n_displayed_rows = n_models - n_dropped_rows
    n_displayed_cells = n_displayed_rows * n_detectors

    n_displayed_compatible = sum(
        1 for m in MODEL_REGISTRY if m not in dropped_models
        for d in DETECTOR_REGISTRY
        if check_compatibility(m, d).compatible
    )
    n_displayed_rejections = n_displayed_cells - n_displayed_compatible

    return {
        'n_models': n_models,
        'n_detectors': n_detectors,
        'n_total_cells': n_total_cells,
        'n_compatible': n_compatible,
        'n_substrate_mismatch': n_substrate_mismatch,
        'n_missing_observable': n_missing_observable,
        'n_total_rejections': n_total_rejections,
        'n_displayed_rows': n_displayed_rows,
        'n_displayed_cells': n_displayed_cells,
        'n_displayed_compatible': n_displayed_compatible,
        'n_displayed_rejections': n_displayed_rejections,
        'display_folds': [
            {'kept': kept, 'dropped': dropped}
            for kept, dropped in DISPLAY_FOLDS
        ],
    }


def format_human_readable(counts: dict[str, Any]) -> str:
    """Format counts as a console-friendly table."""
    lines = [
        '',
        '=== Transfer matrix counts (registry level) ===',
        f"  Models:                  {counts['n_models']}",
        f"  Detectors:               {counts['n_detectors']}",
        f"  Total cells:             {counts['n_total_cells']}",
        f"  Compatible:              {counts['n_compatible']}",
        f"  Substrate-mismatch:      {counts['n_substrate_mismatch']}",
        f"  Missing-observable:      {counts['n_missing_observable']}",
        f"  Total rejections:        {counts['n_total_rejections']}",
        '',
        '=== Transfer matrix counts (displayed in §5.1) ===',
        f"  Display rows:            {counts['n_displayed_rows']}",
        f"  Displayed cells:         {counts['n_displayed_cells']}",
        f"  Displayed compatible:    {counts['n_displayed_compatible']}",
        f"  Displayed rejections:    {counts['n_displayed_rejections']}",
        '',
        '=== Display folds ===',
    ]
    for fold in counts['display_folds']:
        lines.append(f"  {fold['kept']} ← {fold['dropped']}")
    lines.append('')
    return '\n'.join(lines)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description='Emit canonical transfer-matrix counts.'
    )
    parser.add_argument(
        '--json',
        action='store_true',
        help='Emit machine-readable JSON instead of human-readable text.',
    )
    args = parser.parse_args(argv)

    counts = compute_counts()
    if args.json:
        print(json.dumps(counts, indent=2))
    else:
        print(format_human_readable(counts), end='')
    return 0


if __name__ == '__main__':
    sys.exit(main())
