"""Cross-model detection tests — beginning of the transfer matrix.

Tests that detectors correctly identify their home patterns and correctly
reject non-applicable models.
"""

import numpy as np
import pytest


class TestP13OnModels:
    """P13 excitable wave detector across model families."""

    def test_p13_on_greenberg_hastings_detected(self):
        """P13 should detect on its canonical model (broken wave → spiral)."""
        from epc.models.greenberg_hastings import GreenbergHastings
        from epc.detectors.p13_excitable_waves import P13ExcitableWaveDetector

        m = GreenbergHastings(
            rows=50, cols=50, n_states=5, threshold=1,
            init_mode="broken_wave", seed=0,
        )
        h = m.run(500)
        d = P13ExcitableWaveDetector(n_null_runs=99)
        r = d.detect(h, m.get_metadata())
        assert r.detected, f"P13 not detected on GH CA: {r.summary()}"

    def test_p13_on_vicsek_not_detected(self):
        """P13 should not detect on Vicsek (no grid, no excitable states)."""
        from epc.models.vicsek import VicsekModel
        from epc.detectors.p13_excitable_waves import P13ExcitableWaveDetector

        m = VicsekModel(n_particles=100, seed=0)
        h = m.run(100)
        d = P13ExcitableWaveDetector()
        r = d.detect(h, m.get_metadata())
        assert not r.detected, f"P13 false positive on Vicsek: {r.summary()}"

    def test_p13_on_quiescent_gh_not_detected(self):
        """P13 should not detect on quiescent GH (high threshold)."""
        from epc.models.greenberg_hastings import GreenbergHastings
        from epc.detectors.p13_excitable_waves import P13ExcitableWaveDetector

        m = GreenbergHastings(rows=30, cols=30, threshold=9, seed=0)
        h = m.run(50)
        d = P13ExcitableWaveDetector()
        r = d.detect(h, m.get_metadata())
        assert not r.detected, f"P13 false positive on quiescent GH: {r.summary()}"


class TestP5OnModels:
    """P5 flocking detector across model families."""

    def test_p5_on_vicsek_low_noise_detected(self):
        """P5 should detect flocking in low-noise Vicsek."""
        from epc.models.vicsek import VicsekModel
        from epc.detectors.p5_flocking import P5FlockingDetector

        m = VicsekModel(
            n_particles=100, box_size=5.0, speed=0.03,
            noise=0.05, interaction_radius=1.0, seed=42,
        )
        h = m.run(3000)
        d = P5FlockingDetector(n_permutations=99)
        r = d.detect(h, m.get_metadata())
        assert r.detected, f"P5 not detected on low-noise Vicsek: {r.summary()}"

    def test_p5_on_vicsek_high_noise_not_detected(self):
        """P5 should not detect flocking in high-noise Vicsek."""
        from epc.models.vicsek import VicsekModel
        from epc.detectors.p5_flocking import P5FlockingDetector

        m = VicsekModel(
            n_particles=100, box_size=10.0, speed=0.03,
            noise=2.0, interaction_radius=1.0, seed=42,
        )
        h = m.run(500)
        d = P5FlockingDetector()
        r = d.detect(h, m.get_metadata())
        assert not r.detected, f"P5 false positive on high-noise Vicsek: {r.summary()}"

    def test_p5_on_greenberg_hastings_not_detected(self):
        """P5 should not detect on GH CA (no velocities)."""
        from epc.models.greenberg_hastings import GreenbergHastings
        from epc.detectors.p5_flocking import P5FlockingDetector

        m = GreenbergHastings(rows=30, cols=30, seed=0)
        h = m.run(50)
        d = P5FlockingDetector()
        r = d.detect(h, m.get_metadata())
        assert not r.detected, f"P5 false positive on GH CA: {r.summary()}"


class TestP6OnModels:
    """P6 milling detector across model families."""

    def test_p6_on_greenberg_hastings_not_detected(self):
        """P6 should not detect on GH CA."""
        from epc.models.greenberg_hastings import GreenbergHastings
        from epc.detectors.p6_milling import P6MillingDetector

        m = GreenbergHastings(rows=30, cols=30, seed=0)
        h = m.run(50)
        d = P6MillingDetector()
        r = d.detect(h, m.get_metadata())
        assert not r.detected, f"P6 false positive on GH CA: {r.summary()}"

    def test_p6_on_standard_vicsek_not_detected(self):
        """P6 should not detect milling in standard Vicsek (no attraction)."""
        from epc.models.vicsek import VicsekModel
        from epc.detectors.p6_milling import P6MillingDetector

        m = VicsekModel(
            n_particles=100, box_size=10.0, speed=0.03,
            noise=0.1, seed=42,
        )
        h = m.run(500)
        d = P6MillingDetector()
        r = d.detect(h, m.get_metadata())
        # Standard Vicsek typically doesn't produce milling
        # At most screening-level is acceptable
        if r.detected:
            assert r.tier.value == "screening", (
                f"P6 unexpectedly confirmed on standard Vicsek: {r.summary()}"
            )


class TestP24CrossModel:
    """P24 homeostasis detector — T1b cross-model generalization test.

    The detector must fire on a SECOND, independent homeostat implementation
    (integral controller) that it was NOT tuned on. This is the minimum
    evidence that the detector recognizes the *phenomenon* (homeostatic
    regulation), not its *implementation* (proportional control).
    """

    def test_p24_on_integral_homeostat_detected(self):
        """P24 must detect on integral controller (independent implementation)."""
        from epc.models.homeostasis import IntegralHomeostat, HomeostatParams, PerturbationSchedule
        from epc.detectors.p24_homeostasis import P24HomeostasisDetector

        params = HomeostatParams(gain=2.0, setpoint=10.0, seed=42, dt=0.1)
        schedule = PerturbationSchedule(onset=50.0, amplitude=5.0)
        model = IntegralHomeostat(params)
        history = model.simulate(2000, schedule=schedule)
        det = P24HomeostasisDetector(n_permutations=199, seed=42)
        result = det.detect(history, model.get_metadata())
        assert result.detected, (
            f"P24 not detected on integral homeostat: {result.summary()}"
        )
        assert result.tier.value in ("confirmation", "definitive"), (
            f"Expected confirmation+, got {result.tier.value}: {result.summary()}"
        )

    def test_p24_integral_vs_proportional_both_definitive(self):
        """Both controller types should reach definitive tier."""
        from epc.models.homeostasis import (
            ProportionalHomeostat, IntegralHomeostat,
            HomeostatParams, PerturbationSchedule,
        )
        from epc.detectors.p24_homeostasis import P24HomeostasisDetector

        schedule = PerturbationSchedule(onset=50.0, amplitude=5.0)
        det = P24HomeostasisDetector(n_permutations=199, seed=42)

        # Proportional
        p_params = HomeostatParams(gain=5.0, setpoint=10.0, seed=42, dt=0.1)
        p_model = ProportionalHomeostat(p_params)
        p_hist = p_model.simulate(2000, schedule=schedule)
        p_result = det.detect(p_hist, p_model.get_metadata())

        # Integral
        i_params = HomeostatParams(gain=2.0, setpoint=10.0, seed=42, dt=0.1)
        i_model = IntegralHomeostat(i_params)
        i_hist = i_model.simulate(2000, schedule=schedule)
        i_result = det.detect(i_hist, i_model.get_metadata())

        assert p_result.tier.value == "definitive", (
            f"Proportional not definitive: {p_result.summary()}"
        )
        assert i_result.tier.value == "definitive", (
            f"Integral not definitive: {i_result.summary()}"
        )


class TestP23CrossModel:
    """P23 anti-coordination detector — T1b cross-model generalization test.

    The detector must fire on the El Farol Bar variant (Arthur 1994), which
    it was NOT tuned on. This is the minimum evidence that the detector
    recognizes the *phenomenon* (anti-coordination / load balancing), not
    its *implementation* (Minority Game strategy lookup tables).
    """

    def test_p23_on_el_farol_detected(self):
        """P23 must detect on El Farol (independent implementation)."""
        from epc.models.minority_game import ElFarolBar, ElFarolParams
        from epc.detectors.p23_anticoordination import P23AnticoordinationDetector

        ef = ElFarolBar(ElFarolParams(
            n_agents=100, capacity=60, n_rounds=1000, seed=42,
        ))
        history = ef.simulate()
        det = P23AnticoordinationDetector(n_permutations=199, seed=42)
        result = det.detect(history, ef.get_metadata())
        assert result.detected, (
            f"P23 not detected on El Farol: {result.notes}"
        )
        assert result.tier.value in ("confirmation", "definitive"), (
            f"Expected confirmation+, got {result.tier.value}: {result.notes}"
        )

    def test_p23_both_models_detected(self):
        """Both MG and El Farol should be detected."""
        from epc.models.minority_game import (
            MinorityGame, MinorityGameParams,
            ElFarolBar, ElFarolParams,
        )
        from epc.detectors.p23_anticoordination import P23AnticoordinationDetector

        det = P23AnticoordinationDetector(n_permutations=199, seed=42)

        # Minority Game
        mg = MinorityGame(MinorityGameParams(
            n_agents=101, memory=6, n_rounds=1000, seed=42,
        ))
        mg_result = det.detect(mg.simulate(), mg.get_metadata())

        # El Farol
        ef = ElFarolBar(ElFarolParams(
            n_agents=100, capacity=60, n_rounds=1000, seed=42,
        ))
        ef_result = det.detect(ef.simulate(), ef.get_metadata())

        assert mg_result.tier.value == "definitive", (
            f"MG not definitive: {mg_result.summary()}"
        )
        assert ef_result.detected, (
            f"El Farol not detected: {ef_result.summary()}"
        )


class TestP26CrossModel:
    """P26 stochastic resonance detector — T1b cross-model generalization test.

    The detector must fire on a SECOND, independent SR implementation
    (threshold unit) that it was NOT tuned on. This is the minimum
    evidence that the detector recognizes the *phenomenon* (stochastic
    resonance), not its *implementation* (bistable double-well).
    """

    def test_p26_on_threshold_unit_detected(self):
        """P26 must detect on threshold unit (independent implementation)."""
        from epc.models.stochastic_resonance import ThresholdUnit, ThresholdUnitParams
        from epc.detectors.p26_stochastic_resonance import P26StochasticResonanceDetector

        params = ThresholdUnitParams(seed=42)
        model = ThresholdUnit(params)
        history = model.simulate()
        det = P26StochasticResonanceDetector(n_permutations=199, seed=42)
        result = det.detect(history, model.get_metadata())
        assert result.detected, (
            f"P26 not detected on threshold unit: {result.notes}"
        )
        assert result.tier.value in ("confirmation", "definitive"), (
            f"Expected confirmation+, got {result.tier.value}"
        )

    def test_p26_both_models_definitive(self):
        """Both SR implementations should reach definitive tier."""
        from epc.models.stochastic_resonance import (
            BistableDoubleWell, ThresholdUnit,
            DoubleWellParams, ThresholdUnitParams,
        )
        from epc.detectors.p26_stochastic_resonance import P26StochasticResonanceDetector

        det = P26StochasticResonanceDetector(n_permutations=199, seed=42)

        # Double-well
        dw_params = DoubleWellParams(seed=42, n_trials=20, n_steps=20000)
        dw_model = BistableDoubleWell(dw_params)
        dw_hist = dw_model.simulate()
        dw_result = det.detect(dw_hist, dw_model.get_metadata())

        # Threshold unit
        tu_params = ThresholdUnitParams(seed=42)
        tu_model = ThresholdUnit(tu_params)
        tu_hist = tu_model.simulate()
        tu_result = det.detect(tu_hist, tu_model.get_metadata())

        assert dw_result.tier.value == "definitive", (
            f"Double-well not definitive: {dw_result.tier.value}"
        )
        assert tu_result.tier.value == "definitive", (
            f"Threshold unit not definitive: {tu_result.tier.value}"
        )


class TestP16CrossModel:
    """P16 associative memory detector — T1b cross-model generalization test.

    The detector must fire on a SECOND, independent attractor system
    (Boolean GRN) that it was NOT tuned on. This is the minimum
    evidence that the detector recognizes the *phenomenon* (associative
    memory / pattern completion), not its *implementation* (Hopfield).
    """

    def test_p16_on_boolean_grn_detected(self):
        """P16 must detect on Boolean GRN (independent implementation)."""
        from epc.models.hopfield import BooleanGRN, BooleanGRNParams
        from epc.detectors.p16_associative_memory import P16AssociativeMemoryDetector

        params = BooleanGRNParams(N=50, P=4, corruption=0.2, seed=42)
        model = BooleanGRN(params)
        history = model.simulate(n_cues=4, cue_pattern_indices=[0, 1, 2, 3])
        det = P16AssociativeMemoryDetector(n_permutations=199, seed=42)
        result = det.detect(history, model.get_metadata())
        assert result.detected, (
            f"P16 not detected on Boolean GRN: {result.notes}"
        )
        assert result.tier.value in ("confirmation", "definitive"), (
            f"Expected confirmation+, got {result.tier.value}"
        )

    def test_p16_both_models_detected(self):
        """Both Hopfield and Boolean GRN should be detected."""
        from epc.models.hopfield import (
            HopfieldNetwork, HopfieldParams,
            BooleanGRN, BooleanGRNParams,
        )
        from epc.detectors.p16_associative_memory import P16AssociativeMemoryDetector

        det = P16AssociativeMemoryDetector(n_permutations=199, seed=42)

        # Hopfield
        hop = HopfieldNetwork(HopfieldParams(N=100, P=5, corruption=0.2, seed=42))
        hop_result = det.detect(
            hop.simulate(n_cues=5, cue_pattern_indices=[0, 1, 2, 3, 4]),
            hop.get_metadata(),
        )

        # Boolean GRN
        grn = BooleanGRN(BooleanGRNParams(N=50, P=4, corruption=0.2, seed=42))
        grn_result = det.detect(
            grn.simulate(n_cues=4, cue_pattern_indices=[0, 1, 2, 3]),
            grn.get_metadata(),
        )

        assert hop_result.tier.value == "definitive", (
            f"Hopfield not definitive: {hop_result.summary()}"
        )
        assert grn_result.detected, (
            f"Boolean GRN not detected: {grn_result.summary()}"
        )


class TestP25OnModels:
    """P25 equifinality detector across model families (T1b)."""

    def test_p25_on_canalized_landscape_detected(self):
        """P25 should detect on canonical CanalizedLandscape."""
        from epc.models.canalization import CanalizedLandscape, CanalizedLandscapeParams
        from epc.detectors.p25_equifinality import P25EquifinalityDetector

        params = CanalizedLandscapeParams(
            n_dims=10, n_ics=20, n_steps=200, seed=42,
        )
        model = CanalizedLandscape(params)
        h = model.simulate()
        det = P25EquifinalityDetector(n_permutations=199, seed=42)
        r = det.detect(h, model.get_metadata())
        assert r.detected, f"P25 not detected on CanalizedLandscape: {r.summary()}"

    def test_p25_on_multi_basin_grn_detected(self):
        """P25 should detect on MultiBasinGRN (T1b cross-model)."""
        from epc.models.canalization import MultiBasinGRN
        from epc.detectors.p25_equifinality import P25EquifinalityDetector

        grn = MultiBasinGRN(n_genes=10, n_ics=20, n_steps=400, seed=42)
        h = grn.simulate()
        det = P25EquifinalityDetector(n_permutations=199, seed=42)
        r = det.detect(h, grn.get_metadata())
        assert r.detected, f"P25 not detected on MultiBasinGRN: {r.summary()}"

    def test_p25_on_diffusive_not_detected(self):
        """P25 should NOT detect on diffusive dynamics."""
        from epc.models.canalization import DiffusiveDynamics
        from epc.detectors.p25_equifinality import P25EquifinalityDetector

        model = DiffusiveDynamics(n_dims=10, n_ics=20, n_steps=200, seed=42)
        h = model.simulate()
        det = P25EquifinalityDetector(n_permutations=99, seed=42)
        r = det.detect(h, model.get_metadata())
        assert not r.detected, f"P25 false positive on diffusive: {r.summary()}"


class TestP20OnModels:
    """P20 quorum sensing detector across model families (T1b)."""

    def test_p20_on_autoinducer_quorum_detected(self):
        """P20 should detect on canonical AutoinducerQuorum."""
        from epc.models.quorum_sensing import AutoinducerQuorum, AutoinducerParams
        from epc.detectors.p20_quorum_sensing import P20QuorumSensingDetector

        model = AutoinducerQuorum(AutoinducerParams(seed=42))
        h = model.simulate()
        det = P20QuorumSensingDetector(n_permutations=199, seed=42)
        r = det.detect(h, model.get_metadata())
        assert r.detected, f"P20 not detected on AutoinducerQuorum: {r.summary()}"
        assert r.tier.value == "definitive"

    def test_p20_on_fraction_threshold_detected(self):
        """P20 should detect on FractionThresholdModel (T1b cross-model)."""
        from epc.models.quorum_sensing import FractionThresholdModel, FractionThresholdParams
        from epc.detectors.p20_quorum_sensing import P20QuorumSensingDetector

        model = FractionThresholdModel(FractionThresholdParams(seed=42))
        h = model.simulate()
        det = P20QuorumSensingDetector(n_permutations=199, seed=42)
        r = det.detect(h, model.get_metadata())
        assert r.detected, f"P20 not detected on FractionThreshold: {r.summary()}"

    def test_p20_on_graded_not_detected(self):
        """P20 should NOT detect on graded response (no threshold)."""
        from epc.models.quorum_sensing import GradedResponseModel
        from epc.detectors.p20_quorum_sensing import P20QuorumSensingDetector

        model = GradedResponseModel(seed=42)
        h = model.simulate()
        det = P20QuorumSensingDetector(n_permutations=99, seed=42)
        r = det.detect(h, model.get_metadata())
        assert not r.detected, f"P20 false positive on graded: {r.summary()}"


class TestP4OnModels:
    """P4 territoriality detector across model families (T1b)."""

    def test_p4_on_scent_marking_detected(self):
        """P4 should detect on canonical ScentMarkingModel."""
        from epc.models.territoriality import ScentMarkingModel, ScentMarkingParams
        from epc.detectors.p4_territoriality import P4TerritorialityDetector

        model = ScentMarkingModel(ScentMarkingParams(
            n_steps=20000, snapshot_interval=100, seed=42,
        ))
        h = model.simulate()
        det = P4TerritorialityDetector(n_permutations=199, seed=42)
        r = det.detect(h, model.get_metadata())
        assert r.detected, f"P4 not detected on ScentMarking: {r.warnings}"
        assert r.tier.value == "definitive"

    def test_p4_on_pheromone_repulsion_detected(self):
        """P4 should detect on PheromoneRepulsionModel (T1b cross-model)."""
        from epc.models.territoriality import PheromoneRepulsionModel, PheromoneRepulsionParams
        from epc.detectors.p4_territoriality import P4TerritorialityDetector

        model = PheromoneRepulsionModel(PheromoneRepulsionParams(seed=42))
        h = model.simulate()
        det = P4TerritorialityDetector(n_permutations=199, seed=42)
        r = det.detect(h, model.get_metadata())
        assert r.detected, f"P4 not detected on PheromoneRepulsion: {r.warnings}"

    def test_p4_on_random_walk_not_detected(self):
        """P4 should NOT detect on PlainRandomWalkModel."""
        from epc.models.territoriality import PlainRandomWalkModel
        from epc.detectors.p4_territoriality import P4TerritorialityDetector

        model = PlainRandomWalkModel(
            n_agents=4, grid_size=48, n_steps=5000,
            snapshot_interval=50, seed=42,
        )
        h = model.simulate()
        det = P4TerritorialityDetector(n_permutations=199, seed=42)
        r = det.detect(h, model.get_metadata())
        assert not r.detected, f"P4 false positive on random walk: {r.warnings}"


class TestP29OnModels:
    """P29 trail network detector across model families (T1b)."""

    def test_p29_on_ant_trail_detected(self):
        """P29 should detect on canonical AntTrailModel."""
        from epc.models.trail_network import AntTrailModel, AntTrailParams
        from epc.detectors.p29_trail_network import P29TrailNetworkDetector

        model = AntTrailModel(AntTrailParams(
            n_nodes=7, n_agents=50, grid_size=100,
            alpha=1.0, beta=2.0, deposition_rate=10.0,
            evaporation_rate=0.02, n_steps=500, snapshot_interval=50,
            node_layout="grid", seed=42,
        ))
        h = model.simulate()
        det = P29TrailNetworkDetector(n_permutations=199, seed=42)
        r = det.detect(h, model.get_metadata())
        assert r.detected, f"P29 not detected on AntTrail: {r.summary()}"

    def test_p29_on_physarum_detected(self):
        """P29 should detect on PhysarumModel (T1b cross-model)."""
        from epc.models.trail_network import PhysarumModel, PhysarumParams
        from epc.detectors.p29_trail_network import P29TrailNetworkDetector

        model = PhysarumModel(PhysarumParams(
            n_nodes=5, grid_size=80,
            gamma=1.8, decay_rate=0.01,
            n_steps=2000, snapshot_interval=50,
            node_layout="grid", seed=42,
        ))
        h = model.simulate()
        det = P29TrailNetworkDetector(n_permutations=199, seed=42)
        r = det.detect(h, model.get_metadata())
        assert r.detected, f"P29 not detected on Physarum: {r.summary()}"

    def test_p29_on_no_reinforcement_not_detected(self):
        """P29 should NOT detect on NoReinforcementModel."""
        from epc.models.trail_network import NoReinforcementModel, NoReinforcementParams
        from epc.detectors.p29_trail_network import P29TrailNetworkDetector

        model = NoReinforcementModel(NoReinforcementParams(
            n_nodes=7, n_agents=50, grid_size=100,
            n_steps=500, snapshot_interval=50,
            node_layout="grid", seed=42,
        ))
        h = model.simulate()
        det = P29TrailNetworkDetector(n_permutations=99, seed=42)
        r = det.detect(h, model.get_metadata())
        assert not r.detected, f"P29 false positive on no-reinforcement: {r.summary()}"


class TestP32OnModels:
    """P32 specialization detector across model families (Sprint 90 T1b)."""

    def test_p32_on_response_threshold_detected(self):
        """P32 should detect on canonical ResponseThresholdModel."""
        from epc.models.division_of_labor import ResponseThresholdModel
        from epc.detectors.p32_specialization import P32SpecializationDetector

        m = ResponseThresholdModel(
            n_agents=20, n_tasks=3, reinforcement_rate=0.05,
            forgetting_rate=0.01, seed=42,
        )
        h = m.run(1000)
        det = P32SpecializationDetector(n_permutations=199, seed=42)
        r = det.detect(h, m.get_metadata())
        assert r.detected, f"P32 not detected on ResponseThresholdModel: {r.summary()}"

    def test_p32_on_second_dol_model_detected(self):
        """T1b: P32 fires on a second division-of-labor model.

        Uses a different reinforcement scheme: instead of additive threshold
        change, we use a multiplicative scheme where performing a task
        multiplies the threshold by (1 - rate) and not performing multiplies
        by (1 + rate). This tests OOD-readiness via the T1a adapter.
        """
        import numpy as np
        # Build a synthetic history mimicking a multiplicative-threshold model
        # where agents specialize over time
        rng = np.random.default_rng(123)
        n_agents, n_tasks, T = 20, 3, 600

        # Simulate multiplicative threshold reinforcement
        thresholds = np.full((n_agents, n_tasks), 0.5)
        stimulus = np.full(n_tasks, 0.5)
        history = []

        for t in range(T):
            # Response probability
            s2 = stimulus ** 2
            th2 = thresholds ** 2
            denom = s2[np.newaxis, :] + th2
            denom = np.where(denom > 0, denom, 1e-12)
            prob = s2[np.newaxis, :] / denom

            rand = rng.random((n_agents, n_tasks))
            responds = rand < prob

            assignments = np.full(n_agents, -1, dtype=np.int64)
            for i in range(n_agents):
                active = np.where(responds[i])[0]
                if len(active) > 0:
                    assignments[i] = active[np.argmax(prob[i, active])]

            history.append({
                "task_assignments": assignments.copy(),
                "thresholds": thresholds.copy(),
                "stimulus": stimulus.copy(),
                "n_agents": n_agents,
                "n_tasks": n_tasks,
                "step": t,
            })

            # Multiplicative threshold update (different from additive)
            for i in range(n_agents):
                if assignments[i] >= 0:
                    task = assignments[i]
                    thresholds[i, task] *= 0.95  # decrease for performed
                    thresholds[i, task] = max(0.01, thresholds[i, task])
                    for j in range(n_tasks):
                        if j != task:
                            thresholds[i, j] *= 1.02  # increase for unperformed
                            thresholds[i, j] = min(1.0, thresholds[i, j])

            workers = np.zeros(n_tasks)
            for i in range(n_agents):
                if assignments[i] >= 0:
                    workers[assignments[i]] += 1
            stimulus = np.clip(stimulus + 0.1 - 0.1 * workers, 0, 1)

        from epc.detectors.p32_specialization import P32SpecializationDetector
        det = P32SpecializationDetector(n_permutations=199, seed=42)
        r = det.detect(history)
        assert r.detected, f"P32 not detected on multiplicative-threshold model: {r.summary()}"

    def test_p32_on_no_reinforcement_not_detected(self):
        """P32 should NOT detect on NoReinforcementModel (negative control)."""
        from epc.models.division_of_labor import NoReinforcementModel
        from epc.detectors.p32_specialization import P32SpecializationDetector

        m = NoReinforcementModel(n_agents=20, n_tasks=3, seed=42)
        h = m.run(500)
        det = P32SpecializationDetector(n_permutations=99, seed=42)
        r = det.detect(h, m.get_metadata())
        assert not r.detected, f"P32 false positive on no-reinforcement: {r.summary()}"
