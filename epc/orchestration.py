"""
Substrate-aware detector dispatch and orchestration.

Maps models to compatible detectors based on substrate type, preventing
cross-substrate false positives. The transfer matrix is block-diagonal
by substrate type.

7 substrate types:
- lattice_1d: Zhang sorting (chimeric), Nagel-Schreckenberg traffic
- lattice_2d: GH, GoL, BTW sandpile, Schelling, Nowak-May, SIR, RPS, LV
- lattice_2d_continuous: Gray-Scott (continuous-valued field on 2D lattice)
- continuous_2d: Vicsek, D'Orsogna, ABP
- oscillator: Kuramoto
- opinion_space: Hegselmann-Krause
- scalar_wealth: Yard-Sale (Sprint 17, new)

Architecture decision #25 (updated Sprint 17):
  17 models × 16 detectors — compatible pairs identified by substrate.
  Gray-Scott (Sprint 13) occupies the lattice_2d_continuous substrate;
  P3 (Sprint 13) is restricted to it by registration. Nagel-Schreckenberg
  (Sprint 15) shares lattice_1d with Zhang but is the only lattice_1d
  model with a 'velocities' observable; P8 (Sprint 15) is restricted
  to lattice_1d AND requires integer 1D velocities at content-level.
  ABP (Sprint 16) shares continuous_2d with Vicsek and D'Orsogna;
  P2 (Sprint 16) is restricted to continuous_2d and uses metadata
  flags (has_alignment_rule, has_attraction_rule,
  has_density_dependent_speed) to discriminate MIPS from flocking
  and milling at the mechanistic-null (DEFINITIVE) gate.
  Yard-Sale (Sprint 17) occupies the new scalar_wealth substrate —
  the first well-mixed (non-spatial) agent population in the registry.
  P28 (Sprint 17) is restricted to scalar_wealth and uses metadata
  flags (has_conserved_resource, has_multiplicative_stake,
  has_saving_propensity, has_redistribution) to discriminate pure
  condensation from redistributive or saving-propensity regimes at
  the mechanistic-null (DEFINITIVE) gate. See Decisions 47–49.
"""

from __future__ import annotations

from typing import Dict, List, Tuple, Optional, Any
from dataclasses import dataclass


@dataclass
class ModelRegistration:
    """Registration entry for a model."""
    name: str
    substrate_type: str          # lattice_1d, lattice_2d, continuous_2d, oscillator
    observables: List[str]       # Available observable keys in state history
    primary_patterns: List[str]  # Patterns this model is canonical for
    metadata_keys: List[str]     # Available metadata keys


@dataclass
class DetectorRegistration:
    """Registration entry for a detector."""
    pattern_id: str
    required_substrate: List[str]    # Compatible substrate types
    required_observables: List[str]  # Observable keys needed
    observable_scope: str            # state_history_only, model_metadata_assisted, model_metadata_required


@dataclass
class CompatibilityResult:
    """Result of checking model-detector compatibility."""
    compatible: bool
    reason: str                  # 'compatible', 'substrate_mismatch', 'missing_observable', etc.
    model: str
    detector: str


# === Model Registry ===

MODEL_REGISTRY: Dict[str, ModelRegistration] = {
    'zhang_sequential': ModelRegistration(
        name='zhang_sequential',
        substrate_type='lattice_1d',
        observables=['cell_types', 'n', 'positions'],
        primary_patterns=['P1', 'P31'],
        metadata_keys=['algorithm'],
    ),
    'zhang_threaded': ModelRegistration(
        name='zhang_threaded',
        substrate_type='lattice_1d',
        observables=['cell_types', 'n', 'positions'],
        primary_patterns=['P1', 'P31'],
        metadata_keys=['algorithm'],
    ),
    'schelling': ModelRegistration(
        name='schelling',
        substrate_type='lattice_2d',
        observables=['grid', 'grid_dims'],
        primary_patterns=['P1'],
        metadata_keys=['threshold', 'density'],
    ),
    'greenberg_hastings': ModelRegistration(
        name='greenberg_hastings',
        substrate_type='lattice_2d',
        observables=['grid', 'grid_dims'],
        primary_patterns=['P13'],
        metadata_keys=['n_states', 'threshold', 'neighborhood'],
    ),
    'game_of_life': ModelRegistration(
        name='game_of_life',
        substrate_type='lattice_2d',
        observables=['grid', 'grid_dims'],
        primary_patterns=['P15'],
        metadata_keys=['rule'],
    ),
    'btw_sandpile': ModelRegistration(
        name='btw_sandpile',
        substrate_type='lattice_2d',
        observables=['avalanche_sizes', 'avalanche_durations', 'activity'],
        primary_patterns=['P14'],
        metadata_keys=['L', 'z_c'],
    ),
    'vicsek': ModelRegistration(
        name='vicsek',
        substrate_type='continuous_2d',
        observables=['positions', 'velocities', 'headings'],
        primary_patterns=['P5'],
        metadata_keys=['N', 'eta', 'v0', 'r_interact'],
    ),
    'dorsogna': ModelRegistration(
        name='dorsogna',
        substrate_type='continuous_2d',
        observables=['positions', 'velocities'],
        primary_patterns=['P6'],
        metadata_keys=['N', 'C_a', 'C_r', 'l_a', 'l_r', 'alpha', 'beta'],
    ),
    'kuramoto': ModelRegistration(
        name='kuramoto',
        substrate_type='oscillator',
        observables=['theta', 'r', 'psi', 'omega'],
        primary_patterns=['P9'],
        metadata_keys=['N', 'K', 'gamma', 'freq_dist',
                       'has_nonlocal_coupling', 'has_frequency_heterogeneity',
                       'coupling_kernel'],
    ),
    'kuramoto_nonlocal': ModelRegistration(
        name='kuramoto_nonlocal',
        substrate_type='oscillator',
        observables=['theta', 'r', 'psi', 'positions'],
        primary_patterns=['P10'],
        metadata_keys=['N', 'A', 'beta', 'alpha', 'dt', 'record_dt',
                       'kernel', 'init_mode',
                       'has_nonlocal_coupling', 'has_frequency_heterogeneity',
                       'coupling_kernel', 'freq_dist'],
    ),
    'nowak_may': ModelRegistration(
        name='nowak_may',
        substrate_type='lattice_2d',
        observables=['grid', 'grid_dims', 'coop_fraction', 'moran_i'],
        primary_patterns=['P27'],
        metadata_keys=['b', 'pd_structure', 'has_movement'],
    ),
    'hegselmann_krause': ModelRegistration(
        name='hegselmann_krause',
        substrate_type='opinion_space',
        observables=['opinions', 'n_clusters', 'variance'],
        primary_patterns=['P21'],
        metadata_keys=['epsilon', 'n_agents', 'init_mode'],
    ),
    'sir_epidemic': ModelRegistration(
        name='sir_epidemic',
        substrate_type='lattice_2d',
        observables=['grid', 'grid_dims'],
        primary_patterns=['P22'],
        metadata_keys=['infection_prob', 'recovery_prob', 'neighborhood', 'r0_approx'],
    ),
    'rps_spatial': ModelRegistration(
        name='rps_spatial',
        substrate_type='lattice_2d',
        observables=['grid', 'grid_dims'],
        primary_patterns=['P12'],
        metadata_keys=['mobility', 'exchange_rate', 'selection_rate',
                       'reproduction_rate', 'neighborhood', 'dominance_map'],
    ),
    'gray_scott': ModelRegistration(
        name='gray_scott',
        substrate_type='lattice_2d_continuous',
        observables=['field', 'field_u', 'grid_dims'],
        primary_patterns=['P3'],
        metadata_keys=['feed_rate', 'kill_rate', 'Du', 'Dv', 'state_dtype'],
    ),
    'nagel_schreckenberg': ModelRegistration(
        name='nagel_schreckenberg',
        substrate_type='lattice_1d',
        observables=['positions', 'velocities', 'gaps', 'n_cars', 'density',
                     'L', 'v_max', 'mean_velocity', 'flow', 'stopped_fraction'],
        primary_patterns=['P8'],
        metadata_keys=['L', 'v_max', 'p_slow', 'density', 'init_mode',
                       'boundary', 'interaction_type', 'update_mode'],
    ),
    'abp': ModelRegistration(
        name='abp',
        substrate_type='continuous_2d',
        observables=['positions', 'velocities', 'headings', 'speeds',
                     'local_density'],
        primary_patterns=['P2'],
        metadata_keys=['n_particles', 'box_size', 'v0', 'rho_star', 'D_r',
                       'r_cg', 'dt', 'peclet', 'packing_fraction',
                       'has_alignment_rule', 'has_attraction_rule',
                       'has_density_dependent_speed', 'interaction_type'],
    ),
    'yard_sale': ModelRegistration(
        name='yard_sale',
        substrate_type='scalar_wealth',
        observables=['wealth', 'gini', 'max_share', 'top_1pct_share',
                     'top_10pct_share', 'total_wealth'],
        primary_patterns=['P28'],
        metadata_keys=['n_agents', 'f', 'lambda_save', 'chi',
                       'redistribute_every', 'w0', 'init_mode',
                       'has_conserved_resource', 'has_pairwise_exchange',
                       'has_multiplicative_stake', 'has_redistribution',
                       'has_saving_propensity', 'interaction_type',
                       'total_wealth_conserved'],
    ),
}

# === Detector Registry ===

DETECTOR_REGISTRY: Dict[str, DetectorRegistration] = {
    'P1': DetectorRegistration(
        pattern_id='P1',
        required_substrate=['lattice_1d', 'lattice_2d'],
        required_observables=['cell_types', 'grid'],  # either one
        observable_scope='state_history_only',
    ),
    'P5': DetectorRegistration(
        pattern_id='P5',
        required_substrate=['continuous_2d'],
        required_observables=['positions', 'headings'],
        observable_scope='state_history_only',
    ),
    'P6': DetectorRegistration(
        pattern_id='P6',
        required_substrate=['continuous_2d'],
        required_observables=['positions', 'velocities'],
        observable_scope='state_history_only',
    ),
    'P9': DetectorRegistration(
        pattern_id='P9',
        required_substrate=['oscillator'],
        required_observables=['theta'],
        observable_scope='state_history_only',
    ),
    'P13': DetectorRegistration(
        pattern_id='P13',
        required_substrate=['lattice_2d'],
        required_observables=['grid'],
        observable_scope='state_history_only',
    ),
    'P14': DetectorRegistration(
        pattern_id='P14',
        required_substrate=['lattice_2d'],
        required_observables=['avalanche_sizes'],
        observable_scope='model_metadata_required',
    ),
    'P15': DetectorRegistration(
        pattern_id='P15',
        required_substrate=['lattice_2d'],
        required_observables=['grid'],
        observable_scope='state_history_only',
    ),
    'P31': DetectorRegistration(
        pattern_id='P31',
        required_substrate=['lattice_1d'],
        required_observables=['cell_types'],
        observable_scope='state_history_only',
    ),
    'P21': DetectorRegistration(
        pattern_id='P21',
        required_substrate=['opinion_space'],
        required_observables=['opinions'],
        observable_scope='model_metadata_assisted',
    ),
    'P27': DetectorRegistration(
        pattern_id='P27',
        required_substrate=['lattice_2d'],
        required_observables=['coop_fraction'],  # PD-specific; only Nowak-May produces this
        observable_scope='model_metadata_required',
    ),
    'P22': DetectorRegistration(
        pattern_id='P22',
        required_substrate=['lattice_2d'],
        required_observables=['grid'],
        observable_scope='state_history_only',
    ),
    'P12': DetectorRegistration(
        pattern_id='P12',
        required_substrate=['lattice_2d'],
        required_observables=['grid'],
        observable_scope='state_history_only',
    ),
    'P3': DetectorRegistration(
        pattern_id='P3',
        required_substrate=['lattice_2d_continuous'],
        required_observables=['field'],
        observable_scope='state_history_only',
    ),
    'P8': DetectorRegistration(
        pattern_id='P8',
        required_substrate=['lattice_1d'],
        required_observables=['velocities'],
        observable_scope='state_history_only',
    ),
    'P2': DetectorRegistration(
        pattern_id='P2',
        required_substrate=['continuous_2d'],
        required_observables=['positions'],
        observable_scope='model_metadata_assisted',
    ),
    'P28': DetectorRegistration(
        pattern_id='P28',
        required_substrate=['scalar_wealth'],
        required_observables=['wealth'],
        observable_scope='model_metadata_assisted',
    ),
    'P10': DetectorRegistration(
        pattern_id='P10',
        required_substrate=['oscillator'],
        required_observables=['theta'],
        observable_scope='model_metadata_assisted',
    ),
}


def check_compatibility(
    model_name: str,
    detector_id: str,
) -> CompatibilityResult:
    """Check if a model-detector pair is compatible.
    
    Parameters
    ----------
    model_name : str
        Model name from MODEL_REGISTRY.
    detector_id : str
        Detector pattern ID from DETECTOR_REGISTRY.
    
    Returns
    -------
    CompatibilityResult
    """
    if model_name not in MODEL_REGISTRY:
        return CompatibilityResult(
            compatible=False,
            reason=f'unknown_model: {model_name}',
            model=model_name,
            detector=detector_id,
        )
    
    if detector_id not in DETECTOR_REGISTRY:
        return CompatibilityResult(
            compatible=False,
            reason=f'unknown_detector: {detector_id}',
            model=model_name,
            detector=detector_id,
        )
    
    model = MODEL_REGISTRY[model_name]
    detector = DETECTOR_REGISTRY[detector_id]
    
    # Substrate check
    if model.substrate_type not in detector.required_substrate:
        return CompatibilityResult(
            compatible=False,
            reason=f'substrate_mismatch: {model.substrate_type} not in {detector.required_substrate}',
            model=model_name,
            detector=detector_id,
        )
    
    # Observable check (at least one required observable must be present)
    has_observable = any(
        obs in model.observables for obs in detector.required_observables
    )
    if not has_observable:
        return CompatibilityResult(
            compatible=False,
            reason=f'missing_observable: need one of {detector.required_observables}',
            model=model_name,
            detector=detector_id,
        )
    
    return CompatibilityResult(
        compatible=True,
        reason='compatible',
        model=model_name,
        detector=detector_id,
    )


def get_compatible_pairs() -> List[Tuple[str, str]]:
    """Return all compatible (model, detector) pairs."""
    pairs = []
    for model_name in MODEL_REGISTRY:
        for det_id in DETECTOR_REGISTRY:
            result = check_compatibility(model_name, det_id)
            if result.compatible:
                pairs.append((model_name, det_id))
    return pairs


def get_compatibility_matrix() -> Dict[str, Dict[str, CompatibilityResult]]:
    """Return full compatibility matrix as nested dict."""
    matrix = {}
    for model_name in MODEL_REGISTRY:
        matrix[model_name] = {}
        for det_id in DETECTOR_REGISTRY:
            matrix[model_name][det_id] = check_compatibility(model_name, det_id)
    return matrix


def print_compatibility_matrix() -> None:
    """Print human-readable compatibility matrix."""
    det_ids = sorted(DETECTOR_REGISTRY.keys())
    model_names = sorted(MODEL_REGISTRY.keys())
    
    # Header
    header = f"{'Model':<20s}" + "".join(f"{d:>6s}" for d in det_ids)
    print(header)
    print("-" * len(header))
    
    n_compatible = 0
    n_total = 0
    
    for model_name in model_names:
        row = f"{model_name:<20s}"
        for det_id in det_ids:
            result = check_compatibility(model_name, det_id)
            n_total += 1
            if result.compatible:
                row += f"{'✓':>6s}"
                n_compatible += 1
            else:
                row += f"{'×':>6s}"
        print(row)
    
    print(f"\n{n_compatible}/{n_total} compatible pairs, "
          f"{n_total - n_compatible} substrate mismatches")
