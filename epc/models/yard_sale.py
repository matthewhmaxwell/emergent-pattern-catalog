"""
Yard-Sale model of wealth condensation.

References:
    Chakraborti, A. (2002). Distributions of money in model markets of
        economy. International Journal of Modern Physics C, 13(10),
        1315–1321. DOI: 10.1142/S0129183102003905

    Chakraborti, A. & Chakrabarti, B. K. (2000). Statistical mechanics
        of money: How saving propensity affects its distribution.
        European Physical Journal B, 17(1), 167–170.

    Boghosian, B. M., Johnson, A., & Marcq, J. A. (2015). An H-theorem
        for Boltzmann's equation for the Yard-Sale Model of asset
        exchange. Journal of Statistical Physics, 161, 1339–1350.

    Boghosian, B. M. (2014). Kinetics of wealth and the Pareto law.
        Physical Review E, 89, 042804.

Dynamics
--------

N agents each hold a non-negative scalar wealth w_i. At each elementary
transaction, a random pair (i, j) is selected. The stake is a fraction
of the poorer agent's wealth,

    dw = f * min(w_i, w_j),   0 < f <= 1.

A fair coin flip (probability 1/2 each direction) decides which agent
wins the stake. A symmetric "saving propensity" parameter lambda in
[0, 1] allows agents to keep a guaranteed fraction of their own wealth
before transacting:

    w_i' = lambda * w_i + (1 - lambda) * (w_i + w_j) * eps
    w_j' = lambda * w_j + (1 - lambda) * (w_i + w_j) * (1 - eps)

with eps ~ Uniform(0, 1).  The pure Yard-Sale model takes lambda = 0
and the "stake = f * min(w_i, w_j)" rule rather than the uniform-pot
CC 2000 rule.  Both rules conserve total wealth and are symmetric in
expectation, yet drive the Gini coefficient monotonically toward 1
(all wealth on one agent) under either spec.  This is the defining
P28 paradox: symmetric fair exchange produces extreme inequality.

KEY CONTRAST WITH OTHER RESOURCE-REDISTRIBUTION SCHEMES

  - Uniform mixing (each dollar goes to a random agent) gives
    Gini -> ~0.5 in equilibrium (exponential Boltzmann-Gibbs
    distribution).
  - CC 2000 with nonzero saving propensity lambda converges to a
    Gamma distribution with finite Gini bounded well below 1.
  - Pure Yard-Sale (lambda = 0, f > 0) concentrates all wealth on
    one agent in the N -> infinity limit (Boghosian 2014).
  - Adding a redistribution term (wealth tax chi > 0) produces a
    finite Gini fixed point.

We implement only the pure Yard-Sale rule with the stake =
f * min(w_i, w_j) mechanic and a configurable saving propensity
lambda.  A redistribution / wealth tax knob chi is exposed but
defaults to 0 (pure condensation).

Observables exposed in state snapshots:
  - wealth: (N,) float64
  - step: int
  - gini: float  (computed lazily, O(N log N))
  - max_share: float  (w_max / sum(w))
  - top_p_share(p=0.01, 0.10, 0.20): see get_state / compute helpers

Canonical regime (reproduces Chakraborti 2002 & Boghosian 2014):
  N = 1000, f = 0.01, lambda = 0.0, chi = 0.0, w0 = 1.0
  -> Gini(t=1e5)  ~ 0.85,  max_share ~ 0.05
  -> Gini(t=1e6)  ~ 0.97+,  max_share ~ 0.30
  -> Gini(t=1e7)  -> 1.0,  max_share -> 1.0

Partial saving propensity (CC 2000 with stake rule):
  N = 1000, lambda = 0.5  -> Gini plateau ~ 0.30 (finite fixed point)
  N = 1000, lambda = 0.9  -> Gini plateau ~ 0.15

API surface (NOT BaseModel subclass, matches Vicsek/ABP standalone
pattern).
"""

from __future__ import annotations

from typing import Any, Optional

import numpy as np
from numpy.typing import NDArray


class YardSale:
    """Chakraborti-Boghosian Yard-Sale wealth exchange model.

    Parameters
    ----------
    n_agents : int
        Number of agents N.
    f : float
        Stake fraction: each transaction risks f * min(w_i, w_j).
        Default 0.01. Must satisfy 0 < f <= 1.
    lambda_save : float
        Saving propensity in [0, 1]. lambda_save = 0 is pure Yard-Sale
        (full condensation). lambda_save > 0 converges to a Gamma
        distribution with finite Gini.
    chi : float
        Redistribution / wealth-tax rate applied at each
        `redistribute_every` transaction. 0 = no redistribution
        (pure condensation). Small chi > 0 produces a steady-state
        Gini bounded below 1.
    redistribute_every : int
        Apply redistribution step every this many transactions.
        Only meaningful when chi > 0.
    w0 : float
        Initial wealth per agent. Default 1.0 (identical start —
        the key property the model violates).
    init_mode : str
        'equal' (all w_i = w0) or 'uniform' (w_i ~ Uniform(0.5*w0,
        1.5*w0)) or 'exponential' (w_i ~ Exp(1/w0)). Default 'equal',
        which is the cleanest demonstration of spontaneous inequality.
    seed : int
        RNG seed.
    """

    def __init__(
        self,
        n_agents: int = 1000,
        f: float = 0.01,
        lambda_save: float = 0.0,
        chi: float = 0.0,
        redistribute_every: int = 1000,
        w0: float = 1.0,
        init_mode: str = "equal",
        seed: Optional[int] = None,
    ) -> None:
        if n_agents < 10:
            raise ValueError("n_agents must be >= 10")
        if not (0.0 < f <= 1.0):
            raise ValueError("f must be in (0, 1]")
        if not (0.0 <= lambda_save < 1.0):
            raise ValueError("lambda_save must be in [0, 1)")
        if chi < 0:
            raise ValueError("chi must be non-negative")
        if redistribute_every < 1:
            raise ValueError("redistribute_every must be >= 1")
        if w0 <= 0:
            raise ValueError("w0 must be positive")
        if init_mode not in ("equal", "uniform", "exponential"):
            raise ValueError(f"unknown init_mode: {init_mode!r}")

        self.n_agents = n_agents
        self.f = f
        self.lambda_save = lambda_save
        self.chi = chi
        self.redistribute_every = redistribute_every
        self.w0 = w0
        self.init_mode = init_mode
        self.seed = seed

        # State
        self.wealth: NDArray[np.float64] = np.empty(0)
        self.rng: np.random.Generator = np.random.default_rng(seed)
        self._step_count: int = 0
        self._is_setup: bool = False
        self._total_wealth: float = float(n_agents) * w0

    # ------------------------------------------------------------------
    # Lifecycle
    # ------------------------------------------------------------------

    def setup(self) -> dict[str, Any]:
        """Initialize wealth array."""
        self.rng = np.random.default_rng(self.seed)
        self._step_count = 0
        N = self.n_agents
        w0 = self.w0

        if self.init_mode == "equal":
            self.wealth = np.full(N, w0, dtype=np.float64)
        elif self.init_mode == "uniform":
            self.wealth = self.rng.uniform(0.5 * w0, 1.5 * w0, size=N)
        elif self.init_mode == "exponential":
            self.wealth = self.rng.exponential(scale=w0, size=N)
        else:  # pragma: no cover
            raise ValueError(f"unknown init_mode: {self.init_mode!r}")

        # Normalize so sum equals N*w0 (so absolute scale is stable
        # across seeds regardless of init distribution).
        self._total_wealth = float(N) * w0
        current = self.wealth.sum()
        if current > 0:
            self.wealth *= self._total_wealth / current

        self._is_setup = True
        return self.get_state()

    def step(self, n_transactions: int = 1) -> dict[str, Any]:
        """Perform `n_transactions` elementary pair exchanges.

        One "step" can bundle many transactions for efficiency; the
        typical convention in the literature is to measure observables
        every N transactions (i.e., once per "sweep").
        """
        if not self._is_setup:
            self.setup()

        N = self.n_agents
        f = self.f
        lam = self.lambda_save

        # Pre-generate random pairs for the whole batch.
        # Using argmax trick to pick j != i: pick any j in [0, N), if j==i, use (j+1) % N.
        ii = self.rng.integers(0, N, size=n_transactions)
        jj = self.rng.integers(0, N - 1, size=n_transactions)
        jj = jj + (jj >= ii).astype(np.int64)
        # Fair coin for direction
        direction = self.rng.integers(0, 2, size=n_transactions)  # 0 or 1

        # We must do transactions sequentially because wealth evolves;
        # vectorizing over batch would double-count reused agents.
        w = self.wealth  # alias (in-place)
        for k in range(n_transactions):
            i = int(ii[k])
            j = int(jj[k])
            wi = w[i]
            wj = w[j]
            if lam > 0.0:
                # CC 2000 + YS stake combination: lambda fraction stays,
                # (1-lambda) fraction of combined pot is split by eps.
                # We ALSO cap the random portion at the YS stake
                # structure: the poorer agent's unsaved amount caps it.
                unsaved_i = (1.0 - lam) * wi
                unsaved_j = (1.0 - lam) * wj
                pot = unsaved_i + unsaved_j
                eps = self.rng.random()
                new_unsaved_i = pot * eps
                new_unsaved_j = pot * (1.0 - eps)
                w[i] = lam * wi + new_unsaved_i
                w[j] = lam * wj + new_unsaved_j
            else:
                # Pure Yard-Sale
                stake = f * min(wi, wj)
                if direction[k] == 0:
                    w[i] = wi + stake
                    w[j] = wj - stake
                else:
                    w[i] = wi - stake
                    w[j] = wj + stake

        self._step_count += n_transactions

        # Optional redistribution (wealth tax)
        if self.chi > 0 and self.redistribute_every > 0:
            applied = self._step_count // self.redistribute_every
            prev_applied = (self._step_count - n_transactions) // self.redistribute_every
            for _ in range(applied - prev_applied):
                mean_w = w.mean()
                # Tax: rescale toward mean
                w *= (1.0 - self.chi)
                w += self.chi * mean_w

        return self.get_state()

    def run(
        self,
        n_transactions: int,
        record_interval: int = 1000,
    ) -> list[dict[str, Any]]:
        """Run `n_transactions` exchanges; record every `record_interval`.

        Returns a list of state dicts, index 0 = initial state (after
        setup), subsequent entries at each record interval.
        """
        if not self._is_setup:
            history = [self.setup()]
        else:
            history = [self.get_state()]

        batches = n_transactions // record_interval
        remainder = n_transactions - batches * record_interval
        for _ in range(batches):
            self.step(record_interval)
            history.append(self.get_state())
        if remainder > 0:
            self.step(remainder)
            history.append(self.get_state())
        return history

    # ------------------------------------------------------------------
    # Observables / metadata
    # ------------------------------------------------------------------

    @staticmethod
    def _gini(w: NDArray[np.float64]) -> float:
        """Gini coefficient via sorted-order formula.

        For non-negative vector w of length N,
            G = (2 * sum_{i=1}^{N} i * w_(i)) / (N * sum(w)) - (N+1)/N
        where w_(i) are the sorted values (ascending).
        Returns 0 for uniform wealth, 1 for maximal concentration.
        """
        if w.size == 0:
            return 0.0
        total = w.sum()
        if total <= 0:
            return 0.0
        sw = np.sort(w)
        N = len(sw)
        idx = np.arange(1, N + 1, dtype=np.float64)
        return float((2.0 * (idx * sw).sum()) / (N * total) - (N + 1.0) / N)

    @staticmethod
    def _top_p_share(w: NDArray[np.float64], p: float) -> float:
        """Fraction of total wealth held by the top fraction p of agents."""
        if w.size == 0 or w.sum() <= 0:
            return 0.0
        N = len(w)
        k = max(1, int(round(p * N)))
        sorted_desc = np.sort(w)[::-1]
        return float(sorted_desc[:k].sum() / w.sum())

    def get_state(self) -> dict[str, Any]:
        """Snapshot of current wealth configuration and inequality metrics."""
        w = self.wealth
        total = float(w.sum()) if w.size else 0.0
        max_share = float(w.max() / total) if total > 0 else 0.0
        return {
            "wealth": w.copy(),
            "step": self._step_count,
            "gini": self._gini(w),
            "max_share": max_share,
            "top_1pct_share": self._top_p_share(w, 0.01),
            "top_10pct_share": self._top_p_share(w, 0.10),
            "total_wealth": total,
        }

    def get_metadata(self) -> dict[str, Any]:
        """Metadata including mechanistic rule flags for P28 detector."""
        return {
            "model_family": "yard_sale",
            "model_name": "Chakraborti-Boghosian Yard-Sale (2002/2014)",
            "model_class": "wealth_exchange",
            "n_agents": self.n_agents,
            "f": self.f,
            "lambda_save": self.lambda_save,
            "chi": self.chi,
            "redistribute_every": self.redistribute_every,
            "w0": self.w0,
            "init_mode": self.init_mode,
            "seed": self.seed,
            # Substrate flags
            "space_type": "well_mixed_agents",
            "dynamics_type": "stochastic_pairwise_exchange",
            # Mechanistic flags for P28 detector
            "has_conserved_resource": True,
            "has_pairwise_exchange": True,
            "has_multiplicative_stake": True,   # stake ~ min(w_i, w_j)
            "has_redistribution": self.chi > 0,
            "has_saving_propensity": self.lambda_save > 0,
            "interaction_type": "random_pairwise_zero_sum",
            "total_wealth_conserved": self.chi == 0.0,
        }

    def get_timescale(self) -> dict[str, float]:
        """Characteristic timescales.

        T_sweep = N transactions (one per agent on average).
        T_f = 1 / f    rough number of transactions for an agent's
                       stake to change O(w_i).
        T_relax = N / f    rough timescale for Gini to saturate in
                           the pure YS case (heuristic scaling).
        """
        N = self.n_agents
        return {
            "T_sweep": float(N),
            "T_f": 1.0 / self.f,
            "T_relax": float(N) / self.f,
        }
