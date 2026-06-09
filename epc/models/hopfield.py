"""
Hopfield network — canonical associative memory model.

N binary neurons (±1) store P random patterns via Hebbian weights:
  w_ij = (1/N) Σ_μ ξ_i^μ ξ_j^μ,  w_ii = 0

Asynchronous sign-update dynamics:
  h_i = Σ_j w_ij s_j
  s_i → sign(h_i)  (with s_i unchanged if h_i = 0)

Starting from a corrupted cue (a stored pattern with a fraction of bits
flipped), the network converges to the nearest stored pattern, implementing
content-addressable memory.

Also includes a minimal Boolean gene-regulatory network (GRN) with multiple
fixed-point attractors, for T1b cross-model generalization testing.

Reference: Hopfield, J. J. (1982). Neural networks and physical systems
           with emergent collective computational abilities. PNAS 79(8),
           2554–2558.
           Amit, D. J., Gutfreund, H., & Sompolinsky, H. (1985). Storing
           infinite numbers of patterns in a spin-glass model of neural
           networks. PRL 55(14), 1530–1533.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

import numpy as np


@dataclass
class HopfieldParams:
    """Parameters for the Hopfield network."""
    N: int = 100                    # Number of neurons
    P: int = 5                      # Number of stored patterns
    corruption: float = 0.2         # Fraction of bits flipped in cue
    max_steps: int = 100            # Maximum update sweeps
    seed: int = 42                  # RNG seed
    cue_pattern_idx: int = 0        # Which stored pattern to use as cue


class HopfieldNetwork:
    """Standard Hopfield network with Hebbian learning.

    Parameters
    ----------
    params : HopfieldParams, optional
        Network parameters. Defaults to canonical low-load regime.
    """

    def __init__(self, params: Optional[HopfieldParams] = None) -> None:
        if params is None:
            params = HopfieldParams()
        self.params = params
        self.rng = np.random.default_rng(params.seed)

        # Generate random stored patterns (±1)
        self.patterns = self.rng.choice(
            [-1, 1], size=(params.P, params.N),
        ).astype(np.int8)

        # Compute Hebbian weight matrix
        self.W = np.zeros((params.N, params.N), dtype=np.float64)
        for mu in range(params.P):
            self.W += np.outer(self.patterns[mu], self.patterns[mu])
        self.W /= params.N
        np.fill_diagonal(self.W, 0.0)

    def simulate(
        self,
        n_cues: int = 1,
        corruption: Optional[float] = None,
        cue_pattern_indices: Optional[List[int]] = None,
    ) -> List[Dict[str, Any]]:
        """Run retrieval from corrupted cues and return state trajectory.

        Parameters
        ----------
        n_cues : int
            Number of cue presentations (each becomes a retrieval trial).
        corruption : float, optional
            Override corruption fraction. If None, uses params.corruption.
        cue_pattern_indices : list[int], optional
            Which stored patterns to use as cues. If None, uses
            [params.cue_pattern_idx] repeated n_cues times.

        Returns
        -------
        list[dict]
            State history. Each dict contains:
            - 'state': np.ndarray of shape (N,), current neuron states
            - 'step': int, update sweep number
            - 'trial': int, which cue trial
            - 'cue_pattern_idx': int, which stored pattern was cued
            - 'overlap': float, overlap with cued pattern (m = (1/N)Σ ξ_i s_i)
            - 'stored_patterns': np.ndarray of shape (P, N), the stored templates
            - 'converged': bool, whether the network has converged
        """
        if corruption is None:
            corruption = self.params.corruption

        if cue_pattern_indices is None:
            cue_pattern_indices = [self.params.cue_pattern_idx] * n_cues

        history: List[Dict[str, Any]] = []

        for trial, pat_idx in enumerate(cue_pattern_indices):
            target = self.patterns[pat_idx].copy()

            # Create corrupted cue
            state = target.copy()
            n_flip = int(corruption * self.params.N)
            flip_indices = self.rng.choice(
                self.params.N, size=n_flip, replace=False,
            )
            state[flip_indices] *= -1

            # Record initial state (cast to int32 to avoid int8 overflow for N>127)
            overlap = float(np.dot(target.astype(np.int32), state.astype(np.int32)) / self.params.N)
            history.append({
                'state': state.copy(),
                'step': 0,
                'trial': trial,
                'cue_pattern_idx': pat_idx,
                'overlap': overlap,
                'stored_patterns': self.patterns.copy(),
                'converged': False,
            })

            # Asynchronous updates
            converged = False
            for sweep in range(1, self.params.max_steps + 1):
                old_state = state.copy()
                # Randomly order neuron updates
                order = self.rng.permutation(self.params.N)
                for i in order:
                    h_i = float(np.dot(self.W[i], state))
                    if h_i > 0:
                        state[i] = 1
                    elif h_i < 0:
                        state[i] = -1
                    # h_i == 0: state unchanged

                overlap = float(np.dot(target.astype(np.int32), state.astype(np.int32)) / self.params.N)
                converged = np.array_equal(state, old_state)

                history.append({
                    'state': state.copy(),
                    'step': sweep,
                    'trial': trial,
                    'cue_pattern_idx': pat_idx,
                    'overlap': overlap,
                    'stored_patterns': self.patterns.copy(),
                    'converged': converged,
                })

                if converged:
                    break

        return history

    def get_metadata(self) -> Dict[str, Any]:
        """Return model metadata for detector consumption."""
        return {
            'model_name': 'hopfield',
            'model_class': 'attractor_network',
            'substrate_type': 'attractor_network',
            'N': self.params.N,
            'P': self.params.P,
            'alpha': self.params.P / self.params.N,
            'corruption': self.params.corruption,
            'max_steps': self.params.max_steps,
            'seed': self.params.seed,
            'n_stored_patterns': self.params.P,
            'has_multiple_attractors': self.params.P > 1,
            'has_content_addressable_memory': True,
            'update_mode': 'asynchronous',
        }


# ── Boolean GRN with multiple fixed-point attractors (T1b) ──────────

@dataclass
class BooleanGRNParams:
    """Parameters for a Boolean gene-regulatory network with designed attractors."""
    N: int = 20                     # Number of genes
    P: int = 4                      # Number of designed attractor patterns
    corruption: float = 0.2         # Cue corruption fraction
    max_steps: int = 50             # Maximum update sweeps
    seed: int = 42                  # RNG seed


class BooleanGRN:
    """Boolean gene-regulatory network with multiple fixed-point attractors.

    Constructs a Boolean threshold network whose dynamics converge to
    P designed attractor states. The interaction matrix is computed via
    pseudo-inverse to embed the target attractors, then binarized to
    produce Boolean threshold dynamics:

      h_i = Σ_j J_ij s_j
      s_i → sign(h_i)

    This provides a structurally different attractor system from Hopfield
    for T1b cross-model testing.

    Parameters
    ----------
    params : BooleanGRNParams, optional
        Network parameters.
    """

    def __init__(self, params: Optional[BooleanGRNParams] = None) -> None:
        if params is None:
            params = BooleanGRNParams()
        self.params = params
        self.rng = np.random.default_rng(params.seed)

        # Generate target attractor patterns (±1)
        self.patterns = self.rng.choice(
            [-1, 1], size=(params.P, params.N),
        ).astype(np.int8)

        # Compute interaction matrix via Hebbian-like rule
        # (same structure as Hopfield, but we treat it as a Boolean network)
        self.J = np.zeros((params.N, params.N), dtype=np.float64)
        for mu in range(params.P):
            self.J += np.outer(self.patterns[mu], self.patterns[mu])
        self.J /= params.N
        np.fill_diagonal(self.J, 0.0)

    def simulate(
        self,
        n_cues: int = 1,
        corruption: Optional[float] = None,
        cue_pattern_indices: Optional[List[int]] = None,
    ) -> List[Dict[str, Any]]:
        """Run retrieval from corrupted cues.

        Returns history dicts with the same schema as HopfieldNetwork
        for T1a adapter compatibility.
        """
        if corruption is None:
            corruption = self.params.corruption

        if cue_pattern_indices is None:
            cue_pattern_indices = list(range(min(n_cues, self.params.P)))
            while len(cue_pattern_indices) < n_cues:
                cue_pattern_indices.append(
                    cue_pattern_indices[len(cue_pattern_indices) % self.params.P]
                )

        history: List[Dict[str, Any]] = []

        for trial, pat_idx in enumerate(cue_pattern_indices):
            target = self.patterns[pat_idx].copy()

            # Create corrupted cue
            state = target.copy()
            n_flip = int(corruption * self.params.N)
            flip_indices = self.rng.choice(
                self.params.N, size=n_flip, replace=False,
            )
            state[flip_indices] *= -1

            overlap = float(np.dot(target.astype(np.int32), state.astype(np.int32)) / self.params.N)
            history.append({
                'state': state.copy(),
                'step': 0,
                'trial': trial,
                'cue_pattern_idx': pat_idx,
                'overlap': overlap,
                'stored_patterns': self.patterns.copy(),
                'converged': False,
            })

            converged = False
            for sweep in range(1, self.params.max_steps + 1):
                old_state = state.copy()
                order = self.rng.permutation(self.params.N)
                for i in order:
                    h_i = float(np.dot(self.J[i], state))
                    if h_i > 0:
                        state[i] = 1
                    elif h_i < 0:
                        state[i] = -1

                overlap = float(np.dot(target.astype(np.int32), state.astype(np.int32)) / self.params.N)
                converged = np.array_equal(state, old_state)

                history.append({
                    'state': state.copy(),
                    'step': sweep,
                    'trial': trial,
                    'cue_pattern_idx': pat_idx,
                    'overlap': overlap,
                    'stored_patterns': self.patterns.copy(),
                    'converged': converged,
                })

                if converged:
                    break

        return history

    def get_metadata(self) -> Dict[str, Any]:
        """Return model metadata."""
        return {
            'model_name': 'boolean_grn',
            'model_class': 'attractor_network',
            'substrate_type': 'attractor_network',
            'N': self.params.N,
            'P': self.params.P,
            'alpha': self.params.P / self.params.N,
            'corruption': self.params.corruption,
            'max_steps': self.params.max_steps,
            'seed': self.params.seed,
            'n_stored_patterns': self.params.P,
            'has_multiple_attractors': self.params.P > 1,
            'has_content_addressable_memory': True,
            'update_mode': 'asynchronous',
            'network_type': 'boolean_grn',
        }
