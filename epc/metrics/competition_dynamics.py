"""Metrics for field-mediated resource competition (P37).

Operates on a trajectory of frames each carrying 'abundances' (species, shape (n,))
and 'resources' (shape (k,)). The signature of Huisman-Weissing competition is
SUSTAINED non-equilibrium dynamics with multi-species coexistence: abundances keep
oscillating (do not reach a fixed point) while several species persist together.
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional

import numpy as np


def coexistence_oscillation(history: List[Dict[str, Any]], persist_frac: float = 0.01,
                            late_frac: float = 2.0 / 3.0) -> Optional[Dict[str, float]]:
    """Sustained-oscillation + coexistence summary over the late window.

    Returns None if the frames do not carry both 'abundances' and 'resources'
    (substrate mismatch). Otherwise:
      osc_cv      : median per-species temporal coefficient of variation over the
                    late window (0 = settled fixed point; >0 = still oscillating)
      n_persist   : number of species above persist_frac * global-peak in the late mean
      n_resources : number of limiting resources (the competition 'field')
      n_species   : total species
    """
    ab = [np.asarray(f["abundances"], dtype=float) for f in history
          if isinstance(f, dict) and "abundances" in f]
    res = [np.asarray(f["resources"], dtype=float) for f in history
           if isinstance(f, dict) and "resources" in f]
    if len(ab) < 6 or not res:
        return None
    A = np.array(ab)                                   # (T, n)
    n_res = int(len(res[0]))
    late = A[int(len(A) * late_frac):]
    peak = float(A.max()) or 1.0
    persist_mask = late.mean(axis=0) > persist_frac * peak
    n_persist = int(persist_mask.sum())
    cvs = []
    for i in np.where(persist_mask)[0]:
        s = late[:, i]
        mu = float(s.mean())
        if mu > 1e-9:
            cvs.append(float(s.std() / mu))
    osc = float(np.median(cvs)) if cvs else 0.0
    return {"osc_cv": osc, "n_persist": float(n_persist),
            "n_resources": float(n_res), "n_species": float(A.shape[1])}


# ---------------------------------------------------------------------------
# Sustained-oscillation test (P37 detector v2, post-panel fix 2026-11-15).
#
# coexistence_oscillation() above measures "several species persist and their
# abundances vary". That is shared by white noise, drift and alternation, so it is
# not the Huisman-Weissing signature, which is SUSTAINED OSCILLATION of coexisting
# species competing through shared depletable resources (Huisman & Weissing 1999,
# Nature 402:407). The functions below test that signature directly:
#   * oscillation   = a dominant periodic band in the detrended abundance series,
#                     with at least `min_cycles` cycles in the window and below the
#                     Nyquist band (period-2 alternation excluded), significant
#                     against surrogates that keep the series' amplitude
#                     distribution but destroy its temporal order (permutations);
#   * coexistence   = >= 3 species persisting (unchanged criterion);
#   * coupling      = the oscillation is tied to the resource field: abundance-
#                     resource cross-correlation beyond what independent series
#                     with the same resource spectra give (multivariate
#                     phase-randomised resource surrogates, Theiler et al. 1992).
# coexistence_oscillation() is kept unchanged (the generic emergence channel and
# the standalone validation script use it).
# ---------------------------------------------------------------------------

def _detrend(x: np.ndarray) -> np.ndarray:
    t = np.arange(len(x), dtype=float)
    return x - np.polyval(np.polyfit(t, x, 1), t)


def _band_bins(n: int, min_cycles: int, nyquist_band: float) -> np.ndarray:
    k = np.arange(n // 2 + 1)
    return (k >= min_cycles) & (k / n < nyquist_band)


def band_peak_fraction(resid: np.ndarray, valid: np.ndarray) -> float:
    """Largest share of the detrended variance carried by any 3 adjacent
    periodogram bins that all lie inside the oscillation band."""
    P = np.abs(np.fft.rfft(resid)) ** 2
    total = float(P[1:].sum())
    if total <= 0.0:
        return 0.0
    best = 0.0
    for c in range(1, len(P) - 1):
        if valid[c - 1] and valid[c] and valid[c + 1]:
            best = max(best, float(P[c - 1:c + 2].sum()) / total)
    return best


def _max_lagged_corr(A: np.ndarray, R: np.ndarray, max_lag: int) -> float:
    """max |corr(A[:, i](t), R[:, j](t + lag))| over columns i, j and |lag| <= max_lag."""
    best = 0.0
    n = A.shape[0]
    for lag in range(-max_lag, max_lag + 1):
        a = A[max(0, -lag): n - max(0, lag)]
        r = R[max(0, lag): n - max(0, -lag)]
        if len(a) < 4:
            continue
        a = (a - a.mean(0)) / (a.std(0) + 1e-300)
        r = (r - r.mean(0)) / (r.std(0) + 1e-300)
        c = np.abs(a.T @ r) / len(a)
        best = max(best, float(c.max()))
    return best


def _phase_randomise(R: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """Multivariate phase-randomised surrogate: the same random phases for every
    column, so each column keeps its power spectrum and the columns keep their
    mutual relations, but the alignment with any other series is destroyed."""
    n = R.shape[0]
    F = np.fft.rfft(R, axis=0)
    ph = rng.uniform(0.0, 2.0 * np.pi, size=F.shape[0])
    ph[0] = 0.0
    if n % 2 == 0:
        ph[-1] = 0.0
    return np.fft.irfft(F * np.exp(1j * ph)[:, None], n=n, axis=0)


def sustained_oscillation_test(history: List[Dict[str, Any]], n_surrogates: int = 199,
                               seed: int = 42, persist_frac: float = 0.01,
                               window_frac: float = 0.5, min_cycles: int = 3,
                               nyquist_band: float = 0.45, max_lag: int = 2,
                               amp_floor: float = 0.05) -> Optional[Dict[str, Any]]:
    """Community oscillation significance + resource coupling over the late window.

    Window = the last ``window_frac`` of the frames (the sustained regime). Returns
    None on substrate mismatch (frames without both 'abundances' and 'resources',
    or fewer than 24 such frames).

    Persisting species: late-window mean > persist_frac x global peak. For each
    persisting species i, on the linearly detrended window series r_i:
      amp_i  = std(r_i) / mean_i                       (oscillation amplitude)
      S_i    = band_peak_fraction(r_i)                 (variance share of the dominant
               3-bin band, >= min_cycles cycles per window, frequency < nyquist_band)
    Community oscillation: S_pool = mean_i S_i over fluctuating species
    (amp_i > 1e-9). Null: n_surrogates draws in which every species series is
    independently permuted (amplitude distribution kept, temporal order destroyed);
    p_osc = (1 + #{S_pool_null >= S_pool}) / (n_surrogates + 1), z_osc = standardised
    S_pool against the same null. Per-species p_i come from the same draws.
    Coupling: C = max |corr(r_i(t), rho_j(t + lag))| over fluctuating species i,
    non-constant detrended resources rho_j and |lag| <= max_lag; null: multivariate
    phase-randomised resource surrogates; p_coupling (None if the resource field is
    constant over the window).
    """
    ab = [np.asarray(f["abundances"], dtype=float) for f in history
          if isinstance(f, dict) and "abundances" in f and "resources" in f]
    rs = [np.asarray(f["resources"], dtype=float) for f in history
          if isinstance(f, dict) and "abundances" in f and "resources" in f]
    if len(ab) < 24:
        return None
    A = np.array(ab)
    R = np.array(rs)
    T = len(A)
    n_win = int(round(T * window_frac))
    W = A[T - n_win:]
    RW = R[T - n_win:]
    rng = np.random.default_rng(seed)
    peak = float(A.max()) or 1.0
    persist = W.mean(axis=0) > persist_frac * peak
    valid = _band_bins(n_win, min_cycles, nyquist_band)
    h = n_win // 2
    species: List[Dict[str, Any]] = []
    resid_cols: List[np.ndarray] = []
    for i in np.where(persist)[0]:
        x = W[:, i]
        mu = float(x.mean())
        r = _detrend(x)
        amp = float(r.std() / mu) if mu > 1e-12 else 0.0
        a1, a2 = float(r[:h].std()), float(r[h:].std())
        species.append({"species": int(i), "amp": amp, "band_fraction": band_peak_fraction(r, valid),
                        "sustained": bool(a1 > 0 and a2 >= 0.5 * a1), "p": 1.0,
                        "fluctuating": bool(amp > 1e-9 and r.std() > 1e-12),
                        "material": bool(amp > amp_floor)})
        if species[-1]["fluctuating"]:
            resid_cols.append(r)
    out: Dict[str, Any] = {"n_window": n_win, "n_persist": int(persist.sum()),
                           "n_species": int(A.shape[1]), "n_resources": int(R.shape[1]),
                           "species": species, "S_pool": None, "p_osc": None, "z_osc": None,
                           "coupling_corr": None, "p_coupling": None}
    fl = [sp for sp in species if sp["fluctuating"]]
    if not fl:
        return out
    S_obs = np.array([sp["band_fraction"] for sp in fl])
    null = np.array([[band_peak_fraction(rng.permutation(r), valid) for r in resid_cols]
                     for _ in range(n_surrogates)])                # (n_surrogates, n_fluct)
    for j, sp in enumerate(fl):
        sp["p"] = float((1 + np.sum(null[:, j] >= S_obs[j])) / (n_surrogates + 1))
    S_pool = float(S_obs.mean())
    pool_null = null.mean(axis=1)
    sd = float(pool_null.std())
    out["S_pool"] = S_pool
    out["p_osc"] = float((1 + np.sum(pool_null >= S_pool)) / (n_surrogates + 1))
    out["z_osc"] = float((S_pool - pool_null.mean()) / sd) if sd > 1e-12 else 0.0
    Rres = np.column_stack([_detrend(RW[:, j]) for j in range(RW.shape[1])])
    keep = Rres.std(axis=0) > 1e-9 * (np.abs(RW).mean(axis=0) + 1e-12)
    if keep.any():
        Ares = np.column_stack(resid_cols)
        Rk = Rres[:, keep]
        C_obs = _max_lagged_corr(Ares, Rk, max_lag)
        cnull = np.array([_max_lagged_corr(Ares, _phase_randomise(Rk, rng), max_lag)
                          for _ in range(n_surrogates)])
        out["coupling_corr"] = C_obs
        out["p_coupling"] = float((1 + np.sum(cnull >= C_obs)) / (n_surrogates + 1))
    return out
