"""
Minimal Schelling Segregation Model.

Canonical model for P1 (similarity-driven aggregation). Agents on a 2D grid
with two types move to random empty cells if fewer than `threshold` fraction
of their neighbors are same-type.

This provides a REAL positive control for the P1 temporal convergence guard:
- Starts with random placement (low Moran's I)
- Converges to clustered state (high Moran's I)
- Temporal signature: monotonic increase then plateau

Reference: Schelling, T.C. (1971). Dynamic models of segregation.
Journal of Mathematical Sociology, 1(2), 143-186.
"""

from __future__ import annotations

import numpy as np
from typing import List, Dict


def run_schelling(
    grid_size: int = 50,
    density: float = 0.9,
    threshold: float = 0.375,
    n_steps: int = 200,
    seed: int = 42,
) -> List[Dict]:
    """Run Schelling segregation model.
    
    Parameters
    ----------
    grid_size : int
        Side length of square grid.
    density : float
        Fraction of cells occupied.
    threshold : float
        Minimum same-type neighbor fraction for satisfaction.
        Schelling's original: 3/8 = 0.375 (3 of 8 neighbors).
    n_steps : int
        Number of time steps.
    seed : int
        RNG seed.
    
    Returns
    -------
    List of state dicts with 'grid' key (2D array: 0=empty, 1=type A, 2=type B).
    """
    rng = np.random.default_rng(seed)
    
    N = grid_size * grid_size
    n_occupied = int(N * density)
    n_type_a = n_occupied // 2
    n_type_b = n_occupied - n_type_a
    
    # Initialize random placement
    cells = np.zeros(N, dtype=np.int8)
    cells[:n_type_a] = 1
    cells[n_type_a:n_type_a + n_type_b] = 2
    rng.shuffle(cells)
    grid = cells.reshape((grid_size, grid_size))
    
    history = [{'grid': grid.copy()}]
    
    for step in range(n_steps):
        # Find unhappy agents
        padded = np.pad(grid, 1, mode='constant', constant_values=0)
        
        unhappy = []
        empty = list(zip(*np.where(grid == 0)))
        
        for r in range(grid_size):
            for c in range(grid_size):
                if grid[r, c] == 0:
                    continue
                
                my_type = grid[r, c]
                neighbors = []
                for dr in [-1, 0, 1]:
                    for dc in [-1, 0, 1]:
                        if dr == 0 and dc == 0:
                            continue
                        val = padded[r+1+dr, c+1+dc]
                        if val > 0:
                            neighbors.append(val)
                
                if len(neighbors) == 0:
                    continue
                
                same_frac = sum(1 for n in neighbors if n == my_type) / len(neighbors)
                if same_frac < threshold:
                    unhappy.append((r, c))
        
        if not unhappy or not empty:
            history.append({'grid': grid.copy()})
            continue
        
        # Move unhappy agents to random empty cells
        rng.shuffle(unhappy)
        rng.shuffle(empty)
        
        n_moves = min(len(unhappy), len(empty))
        for i in range(n_moves):
            r_from, c_from = unhappy[i]
            r_to, c_to = empty[i]
            grid[r_to, c_to] = grid[r_from, c_from]
            grid[r_from, c_from] = 0
            empty[i] = (r_from, c_from)
        
        history.append({'grid': grid.copy()})
    
    return history
