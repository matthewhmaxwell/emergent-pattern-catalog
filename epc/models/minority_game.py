"""Minority Game and El Farol Bar models for P23 anti-coordination.

Implements the standard Minority Game (Challet & Zhang, 1997) and the
El Farol Bar problem (Arthur, 1994). Both exhibit anti-coordination:
agents dynamically distribute across options to avoid overcrowding,
producing attendance variance below the random-choice baseline.

References:
    Challet, D. & Zhang, Y.-C. (1997). Emergence of cooperation and
        organization in an evolutionary game. Physica A, 246, 407–418.
    Arthur, W. B. (1994). Inductive reasoning and bounded rationality
        (the El Farol Bar problem). American Economic Review, 84(2), 406–411.
    Savit, R., Manuca, R. & Riolo, R. (1999). Adaptive competition,
        market efficiency, and phase transitions. Physical Review Letters,
        82(10), 2203–2206.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

import numpy as np


@dataclass
class MinorityGameParams:
    """Parameters for the standard Minority Game."""
    n_agents: int = 101          # Must be odd
    memory: int = 3              # m: agents use last m winning sides
    n_strategies: int = 2        # S: strategies per agent
    n_rounds: int = 500
    seed: int = 42

    def __post_init__(self) -> None:
        if self.n_agents % 2 == 0:
            raise ValueError(f"n_agents must be odd, got {self.n_agents}")
        if self.memory < 1:
            raise ValueError(f"memory must be >= 1, got {self.memory}")


@dataclass
class ElFarolParams:
    """Parameters for the El Farol Bar variant."""
    n_agents: int = 100
    capacity: int = 60           # Bar capacity (threshold for 'go')
    n_predictors: int = 5        # Number of predictors per agent
    history_length: int = 10     # How many past attendance values agents see
    n_rounds: int = 500
    seed: int = 42


class MinorityGame:
    """Standard Minority Game (Challet & Zhang, 1997).

    N agents (odd N) each hold S strategies. A strategy maps the last
    m winning sides (encoded as a binary history of length m → 2^m
    possible histories) to a binary choice {0, 1}. At each round the
    minority side wins; agents update strategy scores. Over time agents
    learn to anti-coordinate, reducing attendance variance below the
    random-choice baseline σ²/N = 1/4.

    The control parameter is α = 2^m / N. For α < α_c ≈ 0.34, agents
    enter the symmetric (maladaptive) phase; for α >> 1 they approach
    random behavior. The minimum variance (efficient phase) occurs near
    α_c.
    """

    def __init__(self, params: MinorityGameParams) -> None:
        self.params = params
        self.N = params.n_agents
        self.m = params.memory
        self.S = params.n_strategies
        self.P = 2 ** self.m  # Number of possible history patterns
        self.alpha = self.P / self.N
        self.rng = np.random.default_rng(params.seed)

        # Initialize strategies: each agent has S strategies
        # Each strategy maps P history patterns to {0, 1}
        self.strategies = self.rng.integers(
            0, 2, size=(self.N, self.S, self.P), dtype=np.int8,
        )

        # Strategy scores (virtual points)
        self.scores = np.zeros((self.N, self.S), dtype=float)

        # Initial history pattern (random m-bit string)
        self.history_bits = self.rng.integers(0, 2, size=self.m, dtype=np.int8)

    def _history_index(self) -> int:
        """Convert binary history to integer index."""
        idx = 0
        for bit in self.history_bits:
            idx = 2 * idx + int(bit)
        return idx

    def simulate(self) -> List[Dict[str, Any]]:
        """Run the Minority Game and return attendance time series.

        Returns
        -------
        list[dict]
            Each dict contains:
            - 'round': int, round number
            - 'attendance': int, number of agents choosing side 1
            - 'n_agents': int, total N
            - 'choices': np.ndarray of shape (N,), agents' choices {0, 1}
            - 'winning_side': int, 0 or 1
        """
        records: List[Dict[str, Any]] = []

        for t in range(self.params.n_rounds):
            mu = self._history_index()

            # Each agent plays their best strategy
            best_strat = np.argmax(self.scores, axis=1)
            choices = self.strategies[
                np.arange(self.N), best_strat, mu
            ]

            attendance = int(np.sum(choices))
            winning_side = 0 if attendance > self.N // 2 else 1

            records.append({
                'round': t,
                'attendance': attendance,
                'n_agents': self.N,
                'choices': choices.copy(),
                'winning_side': winning_side,
            })

            # Update strategy scores: +1 for each strategy that would have
            # predicted the winning side for this history pattern
            for s in range(self.S):
                predictions = self.strategies[:, s, mu]
                self.scores[:, s] += np.where(predictions == winning_side, 1, -1)

            # Update history bits
            self.history_bits = np.roll(self.history_bits, -1)
            self.history_bits[-1] = winning_side

        return records

    def get_metadata(self) -> Dict[str, Any]:
        """Return model metadata."""
        return {
            'model_name': 'minority_game',
            'model_class': 'minority_game',
            'n_agents': self.N,
            'memory': self.m,
            'n_strategies': self.S,
            'n_rounds': self.params.n_rounds,
            'alpha': self.alpha,
            'seed': self.params.seed,
            'substrate': 'choice_timeseries',
            'has_strategy_adaptation': True,
            'has_binary_choice': True,
        }


class ElFarolBar:
    """El Farol Bar problem (Arthur, 1994).

    N agents decide each round whether to attend a bar with capacity C.
    Each agent has K predictors that forecast next attendance from recent
    history. Agents go if their best predictor predicts attendance < C.
    Predictor scores are updated based on prediction accuracy.

    Anti-coordination emerges: attendance fluctuates around capacity with
    variance lower than the random-choice baseline.
    """

    def __init__(self, params: ElFarolParams) -> None:
        self.params = params
        self.N = params.n_agents
        self.C = params.capacity
        self.K = params.n_predictors
        self.H = params.history_length
        self.rng = np.random.default_rng(params.seed)

        # Initialize attendance history with random plausible values
        self.attendance_history = self.rng.integers(
            self.N // 4, 3 * self.N // 4, size=self.H,
        ).tolist()

        # Predictors: linear combinations of recent attendance
        # Each predictor is weights for a linear forecast
        self.predictor_weights = self.rng.standard_normal(
            size=(self.N, self.K, self.H),
        ) * 0.3  # Scale to avoid extreme predictions

        # Predictor offsets
        self.predictor_offsets = self.rng.uniform(
            self.N * 0.3, self.N * 0.7, size=(self.N, self.K),
        )

        # Predictor scores
        self.predictor_scores = np.zeros((self.N, self.K), dtype=float)

    def _predict(self) -> np.ndarray:
        """Each agent's best predictor forecasts next attendance.

        Returns
        -------
        np.ndarray of shape (N,)
            Predicted attendance for each agent.
        """
        recent = np.array(self.attendance_history[-self.H:], dtype=float)
        recent_norm = recent / self.N  # Normalize to [0, 1] range

        # All predictors' forecasts: (N, K)
        forecasts = np.einsum('nkh,h->nk', self.predictor_weights, recent_norm)
        forecasts = forecasts * self.N + self.predictor_offsets

        # Each agent uses their best predictor
        best_k = np.argmax(self.predictor_scores, axis=1)
        predictions = forecasts[np.arange(self.N), best_k]

        return predictions

    def simulate(self) -> List[Dict[str, Any]]:
        """Run the El Farol simulation and return attendance time series.

        Returns
        -------
        list[dict]
            Each dict contains:
            - 'round': int
            - 'attendance': int, number of agents who went
            - 'n_agents': int
            - 'capacity': int
            - 'choices': np.ndarray of shape (N,), 1=go, 0=stay
        """
        records: List[Dict[str, Any]] = []

        for t in range(self.params.n_rounds):
            predictions = self._predict()

            # Go if predicted attendance < capacity
            choices = (predictions < self.C).astype(np.int8)

            attendance = int(np.sum(choices))

            records.append({
                'round': t,
                'attendance': attendance,
                'n_agents': self.N,
                'capacity': self.C,
                'choices': choices.copy(),
            })

            # Update predictor scores based on accuracy
            actual = attendance
            recent = np.array(self.attendance_history[-self.H:], dtype=float)
            recent_norm = recent / self.N

            all_forecasts = np.einsum(
                'nkh,h->nk', self.predictor_weights, recent_norm,
            )
            all_forecasts = all_forecasts * self.N + self.predictor_offsets

            # Score: negative absolute error
            errors = np.abs(all_forecasts - actual)
            self.predictor_scores -= errors  # Accumulate negative error

            self.attendance_history.append(attendance)

        return records

    def get_metadata(self) -> Dict[str, Any]:
        """Return model metadata."""
        return {
            'model_name': 'el_farol',
            'model_class': 'minority_game',
            'n_agents': self.N,
            'capacity': self.C,
            'n_predictors': self.K,
            'history_length': self.H,
            'n_rounds': self.params.n_rounds,
            'seed': self.params.seed,
            'substrate': 'choice_timeseries',
            'has_strategy_adaptation': True,
            'has_capacity_threshold': True,
        }


class RandomChoiceBaseline:
    """Random-choice agents (no strategy adaptation) — negative control.

    Each agent independently chooses side 1 with probability 0.5 each round.
    Attendance follows Binomial(N, 0.5) → variance = N/4.
    """

    def __init__(self, n_agents: int = 101, n_rounds: int = 500,
                 seed: int = 42) -> None:
        if n_agents % 2 == 0:
            raise ValueError(f"n_agents must be odd, got {n_agents}")
        self.N = n_agents
        self.n_rounds = n_rounds
        self.rng = np.random.default_rng(seed)

    def simulate(self) -> List[Dict[str, Any]]:
        """Run random-choice baseline."""
        records: List[Dict[str, Any]] = []
        for t in range(self.n_rounds):
            choices = self.rng.integers(0, 2, size=self.N, dtype=np.int8)
            attendance = int(np.sum(choices))
            records.append({
                'round': t,
                'attendance': attendance,
                'n_agents': self.N,
                'choices': choices.copy(),
                'winning_side': 0 if attendance > self.N // 2 else 1,
            })
        return records

    def get_metadata(self) -> Dict[str, Any]:
        return {
            'model_name': 'random_choice_baseline',
            'model_class': 'random_baseline',
            'n_agents': self.N,
            'n_rounds': self.n_rounds,
            'substrate': 'choice_timeseries',
            'has_strategy_adaptation': False,
        }
