"""P37 — Field-mediated resource competition: sustained competitive oscillation (detector).

v2 — POST-PANEL FIX (2026-11-15, validation rebuild panels v12). The v1 detector
scored D* = median per-species late-window CV gated by >= 3 persisting species.
That is "several species persist and fluctuate", which the full Phase-2a panel
showed is shared by white noise, a monotone drift, a period-2 alternation and
oscillating look-alikes (TNR 12/22, continuous d -0.09;
analysis/outputs/p37_resource_competition_phase2a_panel_v12.json, frozen record).
v2 tests the documented signature instead: species competing through shared
depletable resources that COEXIST THROUGH SUSTAINED OSCILLATION (Huisman &
Weissing 1999, Nature 402:407; 2001, Ecology 82:2682). The fix was designed from
that signature and checked on the five canonical positives only before any panel
negative was scored; the fixed panel result is in-sample (the frozen v1 panel had
been seen). Write-up: docs/validation_rebuild/p33_p37_full_panels.md, "Post-panel
fix".

Measurements (epc.metrics.competition_dynamics.sustained_oscillation_test), over
the last half of the frames:
  coexistence     n_persist >= 3 species (late mean > 1% of the global peak)
  material        >= 3 persisting species with detrended amplitude > 0.05
  oscillation     the community's detrended abundance series concentrate their
                  variance in a dominant periodic band (>= 3 cycles per window,
                  below the Nyquist band, so period-2 alternation is excluded)
                  beyond what n_surrogates permutation surrogates give (same
                  amplitude distribution, temporal order destroyed): p_osc, z_osc.
                  White noise has no band; a linear drift has no detrended
                  amplitude; red noise puts its variance below 3 cycles/window.
  coupling        abundances co-vary with the resource field (max lagged |corr|)
                  beyond multivariate phase-randomised resource surrogates:
                  p_coupling. A 'resources' series unrelated to the abundances,
                  or constant, fails.

Tiers
  screening     coexistence AND material AND p_osc < 0.05 AND p_coupling < 0.05
  confirmation  screening AND p_osc < 0.01 AND p_coupling < 0.01 AND n_persist >= 4
  definitive    confirmation AND median persisting amplitude > 0.10 AND
                >= 3 species sustained (2nd-half amplitude >= 0.5 x 1st-half)

Continuous metric (primary_metric["sustained_oscillation"], the canonical one in
epc/phase2a/continuous_metrics.py from v2 on): z_osc when coexistence, material
and coupling (p < 0.05) hold, else 0. The v1 score is still reported as
primary_metric["coexistence_oscillation"] for comparison; it is NOT the
discriminator any more.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

import numpy as np

from epc.metrics.competition_dynamics import coexistence_oscillation, sustained_oscillation_test

OSC_SCREEN = 0.05     # detrended amplitude floor for a species to count as fluctuating materially
OSC_DEF = 0.10        # definitive: median persisting-species amplitude
MIN_PERSIST = 3       # multi-species coexistence (the exclusion-violation essence)
ALPHA_SCREEN = 0.05   # oscillation and coupling significance, screening
ALPHA_CONF = 0.01     # confirmation


@dataclass
class DetectorResult:
    pattern_id: str
    detected: bool
    tier: str
    confidence: float
    primary_metric: dict
    secondary_metrics: dict
    null_p_value: float = 1.0
    null_type: str = "permutation (oscillation) + phase-randomised resources (coupling)"
    exclusions_checked: list = field(default_factory=list)
    exclusion_results: dict = field(default_factory=dict)
    notes: str = ""


def _none(note: str, legacy: float = 0.0) -> DetectorResult:
    return DetectorResult("P37", False, "none", 0.0,
                          {"sustained_oscillation": 0.0, "coexistence_oscillation": legacy}, {},
                          notes=note)


class P37ResourceCompetitionDetector:
    """Detector for P37 — sustained oscillatory coexistence under resource competition (v2)."""

    def __init__(self, n_surrogates: int = 199, seed: int = 42) -> None:
        self.n_surrogates = n_surrogates
        self.seed = seed

    def detect(self, history: List[Dict[str, Any]],
               metadata: Optional[Dict[str, Any]] = None) -> DetectorResult:
        legacy = coexistence_oscillation(history)
        legacy_score = 0.0
        if legacy is not None and int(legacy["n_persist"]) >= MIN_PERSIST:
            legacy_score = round(float(legacy["osc_cv"]), 4)
        m = sustained_oscillation_test(history, n_surrogates=self.n_surrogates, seed=self.seed,
                                       amp_floor=OSC_SCREEN)
        if m is None:
            return _none("substrate mismatch — P37 needs >= 24 frames of 'abundances'+'resources'.",
                         legacy_score)
        sp = m["species"]
        npr = int(m["n_persist"])
        n_mat = sum(1 for s in sp if s["material"])
        amps = [s["amp"] for s in sp]
        med_amp = float(np.median(amps)) if amps else 0.0
        n_sust = sum(1 for s in sp if s["material"] and s["sustained"])
        p_osc = m["p_osc"] if m["p_osc"] is not None else 1.0
        z_osc = m["z_osc"] if m["z_osc"] is not None else 0.0
        p_cpl = m["p_coupling"] if m["p_coupling"] is not None else 1.0

        coexist = npr >= MIN_PERSIST
        material = n_mat >= MIN_PERSIST
        coupled = p_cpl < ALPHA_SCREEN
        gates = coexist and material and coupled
        Dstar = float(z_osc) if gates else 0.0

        tier, detected, conf = "none", False, 0.0
        if gates and p_osc < ALPHA_SCREEN:
            tier, detected, conf = "screening", True, 0.4
        if detected and p_osc < ALPHA_CONF and p_cpl < ALPHA_CONF and npr >= 4:
            tier, conf = "confirmation", 0.6
        if tier == "confirmation" and med_amp > OSC_DEF and n_sust >= MIN_PERSIST:
            tier, conf = "definitive", 0.9

        return DetectorResult(
            pattern_id="P37", detected=detected, tier=tier, confidence=round(conf, 3),
            primary_metric={"sustained_oscillation": round(Dstar, 4),
                            "coexistence_oscillation": legacy_score},
            secondary_metrics={"n_persist": npr, "n_material": n_mat, "n_sustained": n_sust,
                               "median_amp": round(med_amp, 4),
                               "S_pool": m["S_pool"], "p_osc": p_osc, "z_osc": z_osc,
                               "coupling_corr": m["coupling_corr"], "p_coupling": p_cpl,
                               "n_window": m["n_window"], "n_resources": m["n_resources"],
                               "n_species": m["n_species"],
                               "species_p": [round(s["p"], 4) for s in sp],
                               "species_amp": [round(s["amp"], 4) for s in sp]},
            null_p_value=float(max(p_osc, p_cpl)),
            exclusions_checked=["P11", "P9"],
            exclusion_results={"multi_species_coexistence": coexist,
                               "sustained_oscillation": p_osc < ALPHA_SCREEN,
                               "resource_coupled": coupled},
            notes=f"z_osc={z_osc:.2f} p_osc={p_osc:.3f} p_coupling={p_cpl:.3f} "
                  f"persist={npr}/{m['n_species']} material={n_mat} median_amp={med_amp:.3f}.")
