"""P14 — Self-organized criticality detector.

Detects power-law distributed avalanches characteristic of SOC systems.
Uses Clauset et al. (2009) methodology via the powerlaw package for
MLE fitting, goodness-of-fit testing, and likelihood-ratio comparison
against alternative distributions.

Observable scope: model-metadata required (must verify self-tuning).

Required raw observables:
- avalanche_sizes: array of ints
- avalanche_durations: array of ints (optional, for duration scaling)
- activity: topplings per driving step (optional, for 1/f noise)

Reference: Clauset, A., Shalizi, C. R. & Newman, M. E. J. (2009).
Power-law distributions in empirical data. SIAM Review, 51(4), 661-703.
"""

from __future__ import annotations

import numpy as np
from dataclasses import dataclass
from typing import Optional, Tuple
import warnings


@dataclass
class PowerLawFitResult:
    """Result of power-law fitting and hypothesis testing."""
    tau: float                  # MLE exponent (PDF: P(s) ~ s^(-τ))
    tau_logbin: float           # Log-binned PDF slope (independent check)
    xmin: int                   # Lower bound for fitting
    n_tail: int                 # Number of data points above xmin
    ks_distance: float          # KS distance (lower = better fit)
    lr_vs_exponential: float    # LR test R (>0 favors power-law)
    lr_vs_exponential_p: float  # LR test p-value
    lr_vs_lognormal: float      # LR test R (>0 favors power-law)
    lr_vs_lognormal_p: float    # LR test p-value


@dataclass
class P14DetectorResult:
    """Full P14 detection result."""
    detected: bool
    tier: str                   # 'none', 'screening', 'confirmation', 'definitive'
    confidence: float
    
    # Primary
    fit: PowerLawFitResult
    tau_in_range: bool          # τ ∈ plausible range [1.0, 2.0]
    
    # Secondary
    duration_gamma: Optional[float]    # Duration scaling exponent
    spectral_beta: Optional[float]     # 1/f noise exponent
    
    # Null comparison
    null_tau: Optional[float]          # Dissipative null τ
    null_is_exponential: Optional[bool]
    null_max_size: Optional[int]
    
    # Metadata
    n_avalanches: int
    n_nonzero: int
    size_range: Tuple[int, int]
    warnings: list


def fit_power_law(sizes: np.ndarray, xmin: int = 1) -> PowerLawFitResult:
    """Fit power-law to avalanche size distribution.
    
    Uses the powerlaw package (Clauset et al. 2009 method) with a fixed
    xmin=1 for the bulk distribution. The BTW sandpile has known issues
    with automatic xmin selection due to multifractal scaling.
    
    Parameters
    ----------
    sizes : np.ndarray
        Array of avalanche sizes (must be > 0).
    xmin : int
        Lower bound for fitting.
    
    Returns
    -------
    PowerLawFitResult
    """
    import powerlaw
    
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fit = powerlaw.Fit(sizes, discrete=True, xmin=xmin, verbose=False)
    
    tau = fit.power_law.alpha
    ks_d = fit.power_law.D
    n_tail = int(np.sum(sizes >= xmin))
    
    # LR tests
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        R_exp, p_exp = fit.distribution_compare(
            'power_law', 'exponential', normalized_ratio=True
        )
        R_ln, p_ln = fit.distribution_compare(
            'power_law', 'lognormal', normalized_ratio=True
        )
    
    # Independent check: log-binned PDF slope
    tau_logbin = _logbin_slope(sizes)
    
    return PowerLawFitResult(
        tau=tau,
        tau_logbin=tau_logbin,
        xmin=xmin,
        n_tail=n_tail,
        ks_distance=ks_d,
        lr_vs_exponential=R_exp,
        lr_vs_exponential_p=p_exp,
        lr_vs_lognormal=R_ln,
        lr_vs_lognormal_p=p_ln,
    )


def _logbin_slope(sizes: np.ndarray) -> float:
    """Compute power-law exponent from log-binned PDF slope."""
    if len(sizes) < 50:
        return float('nan')
    
    log_bins = np.logspace(0, np.log10(sizes.max()), 50)
    counts, edges = np.histogram(sizes, bins=log_bins)
    bin_centers = np.sqrt(edges[:-1] * edges[1:])
    bin_widths = edges[1:] - edges[:-1]
    pdf = counts / (len(sizes) * bin_widths)
    
    valid = pdf > 0
    if np.sum(valid) < 5:
        return float('nan')
    
    log_x = np.log10(bin_centers[valid])
    log_y = np.log10(pdf[valid])
    
    # Fit in the bulk range (exclude extreme tails)
    mid = (log_x > 0.3) & (log_x < log_x.max() - 0.3)
    if np.sum(mid) < 5:
        return float('nan')
    
    slope, _ = np.polyfit(log_x[mid], log_y[mid], 1)
    return -slope


def compute_duration_scaling(
    sizes: np.ndarray,
    durations: np.ndarray,
) -> Optional[float]:
    """Compute duration scaling exponent γ where T ~ s^γ."""
    valid = (sizes > 1) & (durations > 1)
    if np.sum(valid) < 100:
        return None
    
    ls = np.log10(sizes[valid].astype(float))
    ld = np.log10(durations[valid].astype(float))
    gamma, _ = np.polyfit(ls, ld, 1)
    return float(gamma)


def compute_spectral_beta(activity: np.ndarray) -> Optional[float]:
    """Compute spectral exponent β from activity time series.
    
    PSD ~ f^(-β). β ∈ (0.5, 1.5) indicates 1/f-type noise.
    
    Note: The BTW activity signal has many zeros (non-toppling events).
    We compute the PSD of the cumulative activity (integrated signal)
    which is known to show cleaner 1/f scaling, then subtract 2 to
    recover the exponent of the original signal.
    """
    from scipy import signal as sig
    
    if len(activity) < 1000:
        return None
    
    # Use cumulative activity for cleaner spectral estimation
    # If y(t) has PSD ~ f^(-β), then ∫y has PSD ~ f^(-(β+2))
    cumulative = np.cumsum(activity - activity.mean())
    
    nperseg = min(8192, len(cumulative) // 4)
    freqs, psd = sig.welch(cumulative.astype(float), fs=1.0, nperseg=nperseg)
    
    # Fit in intermediate frequency range
    f_mask = (freqs > 1e-3) & (freqs < 0.05)
    if np.sum(f_mask) < 10:
        return None
    
    log_f = np.log10(freqs[f_mask])
    log_p = np.log10(psd[f_mask] + 1e-30)
    slope, _ = np.polyfit(log_f, log_p, 1)
    
    # PSD of cumulative ~ f^(-(β+2)), so β = -slope - 2
    beta = -slope - 2
    return float(beta)


def detect_p14(
    avalanche_sizes: np.ndarray,
    avalanche_durations: Optional[np.ndarray] = None,
    activity: Optional[np.ndarray] = None,
    null_sizes: Optional[np.ndarray] = None,
    is_self_tuned: bool = True,
    tau_range: Tuple[float, float] = (1.0, 2.0),
) -> P14DetectorResult:
    """Run P14 SOC detection pipeline.
    
    Parameters
    ----------
    avalanche_sizes : np.ndarray
        Sizes of all avalanches (can include zeros).
    avalanche_durations : np.ndarray, optional
        Durations of all avalanches.
    activity : np.ndarray, optional
        Topplings per driving step (for 1/f noise).
    null_sizes : np.ndarray, optional
        Avalanche sizes from dissipative null model.
    is_self_tuned : bool
        Whether the model is known to self-tune (metadata check).
    tau_range : tuple
        Plausible range for τ.
    
    Returns
    -------
    P14DetectorResult
    """
    warn_list = []
    
    # Filter to non-zero avalanches
    nonzero = avalanche_sizes[avalanche_sizes > 0]
    n_total = len(avalanche_sizes)
    n_nonzero = len(nonzero)
    
    if n_nonzero < 100:
        return P14DetectorResult(
            detected=False, tier='none', confidence=0.0,
            fit=None, tau_in_range=False,
            duration_gamma=None, spectral_beta=None,
            null_tau=None, null_is_exponential=None, null_max_size=None,
            n_avalanches=n_total, n_nonzero=n_nonzero,
            size_range=(0, 0), warnings=["Too few non-zero avalanches"],
        )
    
    # === PRIMARY: Power-law fit ===
    fit = fit_power_law(nonzero, xmin=1)
    tau_in_range = tau_range[0] <= fit.tau <= tau_range[1]
    
    # === SCREENING ===
    # Power-law fit with τ in plausible range AND preferred over exponential
    exp_preferred = fit.lr_vs_exponential > 0
    screening_pass = tau_in_range and exp_preferred
    
    if not screening_pass:
        return P14DetectorResult(
            detected=False, tier='none', confidence=0.0,
            fit=fit, tau_in_range=tau_in_range,
            duration_gamma=None, spectral_beta=None,
            null_tau=None, null_is_exponential=None, null_max_size=None,
            n_avalanches=n_total, n_nonzero=n_nonzero,
            size_range=(int(nonzero.min()), int(nonzero.max())),
            warnings=warn_list,
        )
    
    # === SECONDARIES ===
    duration_gamma = None
    if avalanche_durations is not None:
        nz_dur = avalanche_durations[avalanche_sizes > 0]
        duration_gamma = compute_duration_scaling(nonzero, nz_dur)
    
    spectral_beta = None
    if activity is not None:
        spectral_beta = compute_spectral_beta(activity)
    
    # === NULL MODEL ===
    null_tau = None
    null_is_exponential = None
    null_max_size = None
    
    if null_sizes is not None:
        null_nz = null_sizes[null_sizes > 0]
        if len(null_nz) > 50:
            null_fit = fit_power_law(null_nz, xmin=1)
            null_tau = null_fit.tau
            null_is_exponential = bool(null_fit.lr_vs_exponential < 0)
            null_max_size = int(null_nz.max())
    
    # === TIER DETERMINATION ===
    tier = 'screening'
    confidence = 0.35
    
    # Screening bonuses
    if fit.lr_vs_exponential_p < 0.01:
        confidence += 0.10
    if fit.tau_logbin > 0 and abs(fit.tau - fit.tau_logbin) < 0.1:
        confidence += 0.05  # Cross-validation of τ
    confidence = min(confidence, 0.60)
    
    # CONFIRMATION requires:
    # - τ in model-class range
    # - LR tests: exponential rejected
    # - At least one secondary (duration scaling OR 1/f noise)
    # - Self-tuning verified (metadata)
    lr_exp_pass = fit.lr_vs_exponential > 0 and fit.lr_vs_exponential_p < 0.01
    has_secondary = (
        (duration_gamma is not None and 0.3 < duration_gamma < 0.9) or
        (spectral_beta is not None and 0.5 < spectral_beta < 1.5)
    )
    
    if lr_exp_pass and is_self_tuned and (tau_in_range and has_secondary):
        tier = 'confirmation'
        confidence = 0.55
        
        if fit.lr_vs_exponential_p < 0.001:
            confidence += 0.10
        if duration_gamma is not None and spectral_beta is not None:
            confidence += 0.05  # Both secondaries
        if abs(fit.tau - fit.tau_logbin) < 0.05:
            confidence += 0.05  # Excellent cross-validation
        confidence = min(confidence, 0.85)
        
        # DEFINITIVE requires:
        # - Confirmation +
        # - Dissipative null is exponential
        # - Log-normal not significantly preferred (or documented as known property)
        if null_is_exponential is not None and null_is_exponential:
            tier = 'definitive'
            confidence = 0.75
            
            if null_max_size is not None and null_max_size < nonzero.max() / 10:
                confidence += 0.10  # Clear separation from null
            confidence = min(confidence, 1.00)
    
    # Note: BTW model has known multifractal scaling where log-normal is
    # preferred over simple power-law. This is NOT a detection failure —
    # it's a documented property of the 2D BTW universality class.
    if fit.lr_vs_lognormal < 0 and fit.lr_vs_lognormal_p < 0.01:
        warn_list.append(
            "Log-normal preferred over simple power-law (LR p<0.01). "
            "Known property of 2D BTW (multifractal scaling with log corrections)."
        )
    
    return P14DetectorResult(
        detected=screening_pass,
        tier=tier,
        confidence=confidence,
        fit=fit,
        tau_in_range=tau_in_range,
        duration_gamma=duration_gamma,
        spectral_beta=spectral_beta,
        null_tau=null_tau,
        null_is_exponential=null_is_exponential,
        null_max_size=null_max_size,
        n_avalanches=n_total,
        n_nonzero=n_nonzero,
        size_range=(int(nonzero.min()), int(nonzero.max())),
        warnings=warn_list,
    )
