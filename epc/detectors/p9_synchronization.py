"""P9 — Temporal synchronization detector.

Detects global phase synchronization via the Kuramoto order parameter.
Three-tier detection with uncoupled random-phase null model and P10
(chimera state) exclusion via uniform local order parameter check.

Validated on Kuramoto model:
- Definitive at K=8K_c (r=0.963, 119 T_osc, p=0.005)
- None below K_c (r=0.087, p=0.185)
- P10 excluded with local_r_CV=0.037

Observable scope: State-history only.
"""

from __future__ import annotations

import numpy as np
from dataclasses import dataclass
from typing import Dict, Any, List, Optional, Tuple


@dataclass
class P9DetectorResult:
    """Result of P9 synchronization detection."""
    detected: bool
    tier: str                    # 'none', 'screening', 'confirmation', 'definitive'
    confidence: float
    
    # Primary
    r_mean: float                # Mean order parameter in measurement window
    r_std: float                 # Std of r in measurement window
    n_oscillation_periods: float # Duration in units of T_osc
    
    # Null model
    null_r_mean: float           # Mean r under random-phase null
    null_r_std: float
    p_value: float               # Fraction of null r >= observed r
    n_permutations: int
    
    # P10 exclusion
    local_r_cv: Optional[float]  # CV of local r across spatial windows
    p10_excluded: bool           # True = uniform sync (P9), False = chimera (P10)
    
    # Metadata
    warnings: List[str]


def compute_order_parameter(theta: np.ndarray) -> float:
    """Kuramoto order parameter r = |<exp(iθ)>|."""
    return float(np.abs(np.mean(np.exp(1j * theta))))


def compute_local_r(
    theta: np.ndarray,
    n_windows: int = 10,
) -> np.ndarray:
    """Compute local order parameter in spatial windows.
    
    Oscillators are assumed sorted by natural frequency.
    Divides into n_windows contiguous groups.
    """
    N = len(theta)
    window_size = max(1, N // n_windows)
    local_rs = []
    
    for i in range(0, N, window_size):
        chunk = theta[i:i + window_size]
        if len(chunk) > 1:
            r = float(np.abs(np.mean(np.exp(1j * chunk))))
            local_rs.append(r)
    
    return np.array(local_rs)


def detect_p9(
    state_history: List[Dict[str, Any]],
    n_null_runs: int = 199,
    seed: int = 0,
    model_metadata: Optional[Dict[str, Any]] = None,
) -> P9DetectorResult:
    """Run P9 synchronization detection pipeline.
    
    Parameters
    ----------
    state_history : list of dict
        Each dict must contain 'theta' (N,) array of phases and 'r' float.
    n_null_runs : int
        Number of random-phase null comparisons.
    seed : int
        RNG seed for null model.
    model_metadata : dict, optional
        Model parameters for timescale estimation.
    
    Returns
    -------
    P9DetectorResult
    """
    warnings = []
    
    if len(state_history) < 10:
        return P9DetectorResult(
            detected=False, tier='none', confidence=0.0,
            r_mean=0.0, r_std=0.0, n_oscillation_periods=0.0,
            null_r_mean=0.0, null_r_std=0.0, p_value=1.0,
            n_permutations=n_null_runs, local_r_cv=None,
            p10_excluded=False, warnings=["Too few snapshots"],
        )
    
    # Extract r trajectory
    r_values = np.array([h['r'] for h in state_history])
    N = len(state_history[0]['theta'])
    
    # Use second half as measurement window (after transient)
    mid = len(r_values) // 2
    r_measure = r_values[mid:]
    r_mean = float(np.mean(r_measure))
    r_std = float(np.std(r_measure))
    
    # Estimate oscillation period
    if model_metadata and 'gamma' in model_metadata:
        # For Kuramoto: T_osc ≈ 2π / mean|ω|
        # Approximate: mean|ω| ≈ γ for Lorentzian
        gamma = model_metadata['gamma']
        dt = model_metadata.get('dt', 0.05)
        T_osc_steps = 2 * np.pi / (gamma * dt) if gamma > 0 else 100
    else:
        T_osc_steps = 100  # conservative default
    
    # Recording interval
    if len(state_history) > 1 and 'step' in state_history[0]:
        record_interval = state_history[1]['step'] - state_history[0]['step']
    else:
        record_interval = 1
    
    n_T_osc = (len(r_measure) * record_interval) / T_osc_steps
    
    # === NULL MODEL: random-phase comparison ===
    rng = np.random.default_rng(seed)
    null_rs = np.zeros(n_null_runs)
    
    for i in range(n_null_runs):
        # Random uniform phases → r ≈ 1/√N
        random_theta = rng.uniform(0, 2 * np.pi, N)
        null_rs[i] = compute_order_parameter(random_theta)
    
    null_r_mean = float(np.mean(null_rs))
    null_r_std = float(np.std(null_rs))
    
    # p-value: fraction of null r >= observed r
    p_value = float(np.mean(null_rs >= r_mean))
    if p_value == 0:
        p_value = 1.0 / (n_null_runs + 1)
    
    # === P10 EXCLUSION: local r uniformity check ===
    # Use the last snapshot for spatial analysis
    last_theta = state_history[-1]['theta']
    local_rs = compute_local_r(last_theta, n_windows=10)
    
    if len(local_rs) > 1 and np.mean(local_rs) > 0:
        local_r_cv = float(np.std(local_rs) / np.mean(local_rs))
    else:
        local_r_cv = 0.0
    
    # P10 (chimera): local r is heterogeneous (some high, some low)
    # P9 (full sync): local r is uniformly high → CV is low
    p10_excluded = local_r_cv < 0.2  # uniform → P9, not chimera
    
    # === TIER DETERMINATION ===
    tier = 'none'
    confidence = 0.0
    detected = False
    
    # Screening: r > 0.7 over ≥ 10 T_osc
    if r_mean > 0.7 and n_T_osc >= 10:
        tier = 'screening'
        confidence = 0.35
        detected = True
        
        if p_value < 0.01:
            confidence += 0.10
        confidence = min(confidence, 0.60)
        
        # Confirmation: r > 0.9 over ≥ 50 T_osc AND p < 0.01
        if r_mean > 0.9 and n_T_osc >= 50 and p_value < 0.01:
            tier = 'confirmation'
            confidence = 0.55
            
            if p_value < 0.001:
                confidence += 0.10
            if r_std < 0.05:  # stable sync
                confidence += 0.05
            confidence = min(confidence, 0.85)
            
            # Definitive: ≥ 100 T_osc + P10 exclusion
            if n_T_osc >= 100 and p10_excluded:
                tier = 'definitive'
                confidence = 0.75
                
                if p_value <= 0.005:
                    confidence += 0.10
                if local_r_cv < 0.05:
                    confidence += 0.05
                confidence = min(confidence, 1.00)
    
    return P9DetectorResult(
        detected=detected,
        tier=tier,
        confidence=confidence,
        r_mean=r_mean,
        r_std=r_std,
        n_oscillation_periods=n_T_osc,
        null_r_mean=null_r_mean,
        null_r_std=null_r_std,
        p_value=p_value,
        n_permutations=n_null_runs,
        local_r_cv=local_r_cv,
        p10_excluded=p10_excluded,
        warnings=warnings,
    )
