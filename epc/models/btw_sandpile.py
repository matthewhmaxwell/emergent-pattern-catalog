"""
Bak-Tang-Wiesenfeld (BTW) Sandpile Model.

Canonical model for P14 (self-organized criticality). 2D square lattice with
slow drive (add one grain to random site) and fast relaxation (topple when
z >= z_c = 4, redistribute to 4 neighbors). Open boundary: grains at edges
dissipate.

The system self-tunes to a critical state with power-law distributed
avalanche sizes P(s) ~ s^(-τ) where τ ≈ 1.20 for the 2D BTW model.

Reference: Bak, P., Tang, C. & Wiesenfeld, K. (1987). Self-organized
criticality: An explanation of the 1/f noise. Physical Review Letters,
59(4), 381-384.

Exponent reference: τ ≈ 1.20 (2D lattice), see:
- Lübeck & Usadel (1997), Phys. Rev. E 55, 4095
- Scientific Reports (2021), doi:10.1038/s41598-021-97592-x
"""

from __future__ import annotations

import numpy as np
from dataclasses import dataclass, field
from typing import List, Dict, Any, Optional


@dataclass
class BTWSandpileParams:
    """Parameters for the BTW sandpile model."""
    L: int = 64                    # Lattice side length
    z_c: int = 4                   # Critical height (= number of neighbors)
    n_drive: int = 100_000         # Number of grain additions (driving events)
    n_burn: int = 10_000           # Burn-in driving events (reach SOC state)
    seed: int = 42                 # RNG seed
    record_grid_every: int = 0     # Record full grid every N drive events (0=never)


@dataclass
class AvalancheRecord:
    """Record of a single avalanche."""
    size: int           # Total number of topplings
    duration: int       # Number of parallel update rounds
    area: int           # Number of distinct cells that toppled
    dissipated: int     # Grains lost at boundary


@dataclass
class SandpileResult:
    """Complete result of a sandpile simulation."""
    avalanches: List[AvalancheRecord]
    avalanche_sizes: np.ndarray       # Array of sizes for easy analysis
    avalanche_durations: np.ndarray   # Array of durations
    activity: np.ndarray              # Topplings per drive event
    final_grid: np.ndarray            # Final height configuration
    params: BTWSandpileParams
    n_total_drive: int                # Total driving events run
    n_measured: int                   # Driving events after burn-in
    grid_history: List[Dict]          # Sparse grid snapshots (if recorded)


def run_sandpile(params: Optional[BTWSandpileParams] = None) -> SandpileResult:
    """Run BTW sandpile simulation.
    
    Parameters
    ----------
    params : BTWSandpileParams, optional
        Simulation parameters. Uses defaults if None.
    
    Returns
    -------
    SandpileResult
        Complete simulation results including avalanche statistics.
    """
    if params is None:
        params = BTWSandpileParams()
    
    L = params.L
    z_c = params.z_c
    rng = np.random.default_rng(params.seed)
    
    # Initialize grid — random subcritical heights
    grid = rng.integers(0, z_c, size=(L, L), dtype=np.int32)
    
    avalanches = []
    grid_history = []
    
    n_total = params.n_burn + params.n_drive
    
    for drive_idx in range(n_total):
        # DRIVE: add one grain to random site
        r, c = rng.integers(0, L), rng.integers(0, L)
        grid[r, c] += 1
        
        # RELAX: topple until stable
        size = 0
        duration = 0
        area_mask = np.zeros((L, L), dtype=bool)
        dissipated = 0
        
        while True:
            unstable = grid >= z_c
            if not np.any(unstable):
                break
            
            duration += 1
            n_topple = np.sum(unstable)
            size += n_topple
            area_mask |= unstable
            
            # Topple all unstable sites simultaneously
            topple_count = np.zeros((L, L), dtype=np.int32)
            topple_count[unstable] = grid[unstable] // z_c
            
            # Redistribute: each toppling sends 1 grain to each of 4 neighbors
            redistribute = topple_count * z_c
            grid -= redistribute
            
            # Add to neighbors (with boundary dissipation)
            # North
            grid[1:, :] += topple_count[:-1, :]
            # South
            grid[:-1, :] += topple_count[1:, :]
            # East
            grid[:, 1:] += topple_count[:, :-1]
            # West
            grid[:, :-1] += topple_count[:, 1:]
            
            # Count dissipated grains (edges)
            dissipated += int(np.sum(topple_count[0, :]))   # top row → north
            dissipated += int(np.sum(topple_count[-1, :]))  # bottom → south
            dissipated += int(np.sum(topple_count[:, 0]))   # left col → west
            dissipated += int(np.sum(topple_count[:, -1]))  # right col → east
        
        # Record avalanche (after burn-in)
        if drive_idx >= params.n_burn:
            avalanches.append(AvalancheRecord(
                size=size,
                duration=duration,
                area=int(np.sum(area_mask)),
                dissipated=dissipated,
            ))
            
            # Optional grid recording
            if (params.record_grid_every > 0 and
                (drive_idx - params.n_burn) % params.record_grid_every == 0):
                grid_history.append({
                    'grid': grid.copy(),
                    'drive_idx': drive_idx - params.n_burn,
                })
    
    # Build arrays for analysis
    sizes = np.array([a.size for a in avalanches], dtype=np.int64)
    durations = np.array([a.duration for a in avalanches], dtype=np.int64)
    activity = sizes.copy()  # topplings per drive event
    
    return SandpileResult(
        avalanches=avalanches,
        avalanche_sizes=sizes,
        avalanche_durations=durations,
        activity=activity,
        final_grid=grid.copy(),
        params=params,
        n_total_drive=n_total,
        n_measured=params.n_drive,
        grid_history=grid_history,
    )


def run_dissipative_sandpile(
    params: Optional[BTWSandpileParams] = None,
    p_diss: float = 0.2,
) -> SandpileResult:
    """Run sandpile with bulk dissipation (subcritical null model).
    
    Each toppling event has probability p_diss of losing a grain
    (in addition to normal boundary dissipation). This drives the
    system subcritical, producing exponentially distributed avalanches.
    
    Parameters
    ----------
    params : BTWSandpileParams
        Simulation parameters.
    p_diss : float
        Probability of grain loss per toppling (0 = standard BTW).
    """
    if params is None:
        params = BTWSandpileParams()
    
    L = params.L
    z_c = params.z_c
    rng = np.random.default_rng(params.seed + 1000)
    
    grid = rng.integers(0, z_c, size=(L, L), dtype=np.int32)
    
    avalanches = []
    n_total = params.n_burn + params.n_drive
    
    for drive_idx in range(n_total):
        r, c = rng.integers(0, L), rng.integers(0, L)
        grid[r, c] += 1
        
        size = 0
        duration = 0
        area_mask = np.zeros((L, L), dtype=bool)
        dissipated = 0
        
        while True:
            unstable = grid >= z_c
            if not np.any(unstable):
                break
            
            duration += 1
            n_topple = np.sum(unstable)
            size += n_topple
            area_mask |= unstable
            
            topple_count = np.zeros((L, L), dtype=np.int32)
            topple_count[unstable] = grid[unstable] // z_c
            
            redistribute = topple_count * z_c
            grid -= redistribute
            
            # Bulk dissipation: randomly remove grains
            if p_diss > 0:
                diss_mask = rng.random((L, L)) < p_diss
                diss_grains = topple_count * diss_mask.astype(np.int32)
                dissipated += int(np.sum(diss_grains))
                # Reduce grains sent to neighbors
                topple_count_eff = topple_count - diss_grains
                topple_count_eff = np.maximum(topple_count_eff, 0)
            else:
                topple_count_eff = topple_count
            
            grid[1:, :] += topple_count_eff[:-1, :]
            grid[:-1, :] += topple_count_eff[1:, :]
            grid[:, 1:] += topple_count_eff[:, :-1]
            grid[:, :-1] += topple_count_eff[:, 1:]
            
            dissipated += int(np.sum(topple_count_eff[0, :]))
            dissipated += int(np.sum(topple_count_eff[-1, :]))
            dissipated += int(np.sum(topple_count_eff[:, 0]))
            dissipated += int(np.sum(topple_count_eff[:, -1]))
        
        if drive_idx >= params.n_burn:
            avalanches.append(AvalancheRecord(
                size=size, duration=duration,
                area=int(np.sum(area_mask)), dissipated=dissipated,
            ))
    
    sizes = np.array([a.size for a in avalanches], dtype=np.int64)
    durations = np.array([a.duration for a in avalanches], dtype=np.int64)
    
    return SandpileResult(
        avalanches=avalanches,
        avalanche_sizes=sizes,
        avalanche_durations=durations,
        activity=sizes.copy(),
        final_grid=grid.copy(),
        params=params,
        n_total_drive=n_total,
        n_measured=params.n_drive,
        grid_history=[],
    )
