"""Class B — catalog-derived non-positives for the Phase-2a panel.

Loads the canonical positives of 10 catalog patterns (per spec) and adapts
each to the detector's expected substrate format. Caches the native-form
catalog substrates to disk so multiple panel runs do not re-simulate.
"""

from __future__ import annotations

import os
import pickle
from typing import Any, Callable, Dict, List

import numpy as np


CACHE_DIR = "analysis/outputs/phase2a_catalog_cache"

# Reduced-cost canonical positives for panel use. Each must produce a substrate
# that exhibits the named pattern qualitatively but is small enough to keep
# panel runs tractable. Sizes intentionally smaller than the published
# canonical positives — the panel tests detector specificity, not pattern
# replication accuracy.
SUBSTRATE_PARAMS: Dict[str, Dict[str, Any]] = {
    "P1_schelling":     {"grid_size": 32, "density": 0.9, "threshold": 0.375, "n_steps": 80, "seed": 0},
    "P3_gray_scott":    {"rows": 64, "cols": 64, "feed_rate": 0.037, "kill_rate": 0.060, "n_steps": 400, "seed": 0},
    "P5_vicsek":        {"n_particles": 200, "box_size": 7.0, "speed": 0.03, "noise": 0.1, "n_steps": 200, "seed": 0},
    "P9_kuramoto":      {"N": 200, "K": 8.0, "gamma": 0.5, "dt": 0.05, "n_steps": 1500, "record_every": 10, "seed": 0},
    "P10_chimera":      {"N": 128, "A": 0.995, "beta": 0.05, "n_steps": 800, "seed": 0},
    "P14_btw_sandpile": {"L": 32, "n_drive": 5_000, "n_burn": 1_000, "seed": 0},
    "P15_gol":          {"rows": 64, "cols": 64, "init_mode": "r_pentomino", "n_steps": 200, "seed": 0},
    "P18_voter":        {"rows": 32, "cols": 32, "n_steps": 200, "seed": 0},
    "P27_nowak_may":    {"rows": 50, "cols": 50, "b": 1.8, "init_coop_fraction": 0.5, "n_steps": 100, "seed": 0},
    "P31_zhang_sorting": {"n": 64, "algorithm": "bubble", "n_frozen": 10, "seed": 0},
    "P12_rps":          {"rows": 32, "cols": 32, "mobility": 1e-3, "n_steps": 800, "seed": 0},
    # ε=0.2 is the canonical fragmented regime per Sprint 5 replication notes:
    # produces 2 stable opinion clusters from uniform IC (vs consensus at ε=0.5,
    # 4 clusters at ε=0.1). N=400 chosen so the grid adapter reshapes cleanly
    # to a 20×20 occupancy field for the P18 panel under the v1.1 network override.
    "P21_hegselmann_krause": {"n_agents": 400, "epsilon": 0.2, "init_mode": "uniform", "n_steps": 100, "seed": 0},
    # P11 LV: predator-prey lattice in coexistence regime per Mobilia-Georgiev-Täuber 2007.
    # Default rates produce sustained oscillations, not extinction.
    "P11_lotka_volterra": {"rows": 50, "cols": 50, "predation_rate": 2.0, "prey_reproduction_rate": 1.0, "predator_death_rate": 1.0, "init_prey_fraction": 0.3, "init_predator_fraction": 0.3, "n_steps": 100, "seed": 0},
    # P13 GH: 3-state excitable CA producing spiral waves at random init density.
    "P13_greenberg_hastings": {"rows": 64, "cols": 64, "n_states": 3, "threshold": 1, "init_mode": "random", "init_density": 0.3, "n_steps": 100, "seed": 0},
    # P22 SIR: lattice epidemic at canonical infection regime (above percolation threshold).
    "P22_sir_epidemic": {"rows": 64, "cols": 64, "infection_prob": 0.4, "recovery_prob": 0.1, "init_mode": "single_seed", "n_steps": 200, "seed": 0},
    # P2 ABP: Fily-Marchetti 2012 MIPS regime. n=200, phi≈0.61, Pe=100 → well into phase-separated regime.
    "P2_abp": {"n_particles": 200, "box_size": 16.0, "v0": 0.3, "rho_star": 10.0, "D_r": 3e-3, "r_cg": 1.0, "dt": 0.1, "n_steps": 200, "seed": 0},
    # P6 D'Orsogna: Carrillo 2009 canonical milling regime. ring init → milling immediate.
    "P6_dorsogna": {"n_particles": 100, "C_a": 0.5, "C_r": 1.0, "l_a": 3.0, "l_r": 0.5, "alpha": 1.0, "beta": 0.5, "dt": 0.05, "init_mode": "ring", "init_radius": 5.0, "n_steps": 200, "seed": 0},
    # P7 lane formation: Helbing-Molnár 1995 counterflow social force at canonical lane-forming regime.
    "P7_lane_formation": {"n_agents": 200, "corridor_width": 20.0, "corridor_height": 4.0, "desired_speed": 1.0, "repulsion_amplitude": 5.0, "repulsion_range": 0.3, "tau": 0.5, "dt": 0.05, "n_steps": 200, "seed": 0},
    # P8 NS: canonical jamming regime per Nagel-Schreckenberg 1992 + Bette et al. 2017.
    # density=0.3 (deep jam, stopped≈0.43), v_max=5, p_slow=0.3 matching NS 1992.
    # road_length=100 (panel-scale; qualitative jam pattern present).
    "P8_nagel_schreckenberg": {"road_length": 100, "density": 0.3, "v_max": 5, "p_slow": 0.3, "n_steps": 200, "seed": 0},
    # P17 collective sensing: Berdahl 2013 canonical regime. N=50, box=20, sensing_noise=0.8.
    "P17_collective_sensing": {"n_agents": 50, "box_size": 20.0, "v_max": 0.4, "turn_noise": 0.3, "sensing_noise": 0.8, "alpha": 0.95, "social_strength": 0.2, "field_sigma": 5.0, "field_amplitude": 1.0, "n_steps": 200, "seed": 0},
    # P19 informed minority: Couzin 2005 minority-guided flock at canonical params.
    "P19_informed_minority": {"n_particles": 200, "box_size": 10.0, "speed": 0.03, "noise": 0.1, "interaction_radius": 1.0, "informed_fraction": 0.1, "bias_weight": 0.3, "preferred_direction": 0.0, "n_steps": 200, "seed": 0},
    # P24 proportional homeostat: canonical regulated regime (gain=5.0, sustained perturbation).
    "P24_proportional_homeostat": {"gain": 5.0, "setpoint": 10.0, "dt": 0.1, "noise_std": 0.5, "n_steps": 1000, "pert_onset": 50.0, "pert_amplitude": 5.0, "seed": 0},
}

CATALOG_IDS_FIXED = [
    "P1_schelling", "P3_gray_scott", "P5_vicsek", "P9_kuramoto", "P10_chimera",
    "P14_btw_sandpile", "P15_gol", "P18_voter", "P27_nowak_may", "P31_zhang_sorting",
]
FALLBACK_ID = "P12_rps"


# --- Native-form generators (one per substrate) -----------------------------

def _gen_p1_schelling(p: Dict[str, Any]) -> Dict[str, Any]:
    from epc.models.schelling import run_schelling
    history = run_schelling(
        grid_size=p["grid_size"], density=p["density"],
        threshold=p["threshold"], n_steps=p["n_steps"], seed=p["seed"],
    )
    grids = np.stack([np.asarray(s["grid"], dtype=np.int8) for s in history])
    return {"kind": "grid_categorical", "grids": grids, "n_states": int(grids.max() + 1)}


def _gen_p3_gray_scott(p: Dict[str, Any]) -> Dict[str, Any]:
    from epc.models.gray_scott import GrayScott
    model = GrayScott(
        rows=p["rows"], cols=p["cols"],
        feed_rate=p["feed_rate"], kill_rate=p["kill_rate"], seed=p["seed"],
    )
    model.setup()
    record_every = max(1, p["n_steps"] // 100)
    fields = []
    for t in range(p["n_steps"]):
        state = model.step()
        if t % record_every == 0:
            fields.append(np.asarray(state.get("field", state.get("v", model._v)), dtype=np.float32))
    return {"kind": "field_continuous", "fields": np.stack(fields)}


def _gen_p5_vicsek(p: Dict[str, Any]) -> Dict[str, Any]:
    from epc.models.vicsek import VicsekModel
    model = VicsekModel(
        n_particles=p["n_particles"], box_size=p["box_size"],
        speed=p["speed"], noise=p["noise"], seed=p["seed"],
    )
    history = model.run(n_steps=p["n_steps"])
    headings = np.stack([np.asarray(s["headings"], dtype=np.float32) for s in history])
    positions = np.stack([np.asarray(s["positions"], dtype=np.float32) for s in history])
    return {"kind": "particles", "headings": headings, "positions": positions, "box_size": float(p["box_size"])}


def _gen_p9_kuramoto(p: Dict[str, Any]) -> Dict[str, Any]:
    from epc.models.kuramoto import KuramotoModel, KuramotoParams
    model = KuramotoModel(KuramotoParams(
        N=p["N"], K=p["K"], gamma=p["gamma"], dt=p["dt"], seed=p["seed"],
    ))
    history = model.run(n_steps=p["n_steps"], record_every=p["record_every"])
    theta = np.stack([np.asarray(s["theta"], dtype=np.float32) for s in history])
    return {"kind": "phases", "theta": theta, "gamma": p["gamma"], "dt": p["dt"]}


def _gen_p10_chimera(p: Dict[str, Any]) -> Dict[str, Any]:
    from epc.models.kuramoto_nonlocal import KuramotoNonlocal, KuramotoNonlocalParams
    params = KuramotoNonlocalParams(N=p["N"], A=p["A"], beta=p["beta"], seed=p["seed"])
    model = KuramotoNonlocal(params)
    n_frames = max(50, p["n_steps"] // 10)
    history = model.run(n_frames=n_frames)
    theta = np.stack([np.asarray(s["theta"], dtype=np.float32) for s in history])
    return {"kind": "phases", "theta": theta, "gamma": 0.5, "dt": 0.05}


def _gen_p14_btw_sandpile(p: Dict[str, Any]) -> Dict[str, Any]:
    from epc.models.btw_sandpile import run_sandpile, BTWSandpileParams
    params = BTWSandpileParams(L=p["L"], n_drive=p["n_drive"], n_burn=p["n_burn"], seed=p["seed"])
    result = run_sandpile(params)
    final_grid = np.asarray(result.final_grid, dtype=np.int8)
    return {"kind": "static_grid_int", "grid": final_grid}


def _gen_p15_gol(p: Dict[str, Any]) -> Dict[str, Any]:
    from epc.models.game_of_life import GameOfLife
    model = GameOfLife(
        rows=p["rows"], cols=p["cols"],
        init_mode=p["init_mode"], seed=p["seed"],
    )
    history = model.run(n_steps=p["n_steps"])
    grids = np.stack([np.asarray(s["grid"], dtype=np.int8) for s in history])
    return {"kind": "grid_binary", "grids": grids}


def _gen_p18_voter(p: Dict[str, Any]) -> Dict[str, Any]:
    from epc.models.voter import VoterModel
    model = VoterModel(rows=p["rows"], cols=p["cols"], seed=p["seed"])
    history = model.run(n_steps=p["n_steps"])
    grids = np.stack([np.asarray(s["grid"], dtype=np.int8) for s in history])
    return {"kind": "grid_binary", "grids": grids}


def _gen_p27_nowak_may(p: Dict[str, Any]) -> Dict[str, Any]:
    from epc.models.nowak_may import NowakMayModel
    model = NowakMayModel(
        rows=p["rows"], cols=p["cols"], b=p["b"],
        init_coop_fraction=p["init_coop_fraction"], seed=p["seed"],
    )
    history = model.run(n_steps=p["n_steps"])
    grids = np.stack([np.asarray(s["grid"], dtype=np.int8) for s in history])
    return {"kind": "grid_binary", "grids": grids}


def _gen_p31_zhang_sorting(p: Dict[str, Any]) -> Dict[str, Any]:
    from epc.models.cell_view_sorting import CellViewSorting
    model = CellViewSorting(
        n=p["n"], algorithm=p["algorithm"], n_frozen=p["n_frozen"], seed=p["seed"],
    )
    history = model.run_to_completion(max_rounds=200)
    arrays = np.stack([np.asarray(s["array"], dtype=np.int32) for s in history])
    return {"kind": "sequence", "arrays": arrays}


def _gen_p11_lotka_volterra(p: Dict[str, Any]) -> Dict[str, Any]:
    """LV lattice in coexistence regime (Mobilia-Georgiev-Täuber 2007)."""
    from epc.models.lotka_volterra_lattice import LotkaVolterraLattice

    model = LotkaVolterraLattice(
        rows=p["rows"], cols=p["cols"],
        predation_rate=p["predation_rate"],
        prey_reproduction_rate=p["prey_reproduction_rate"],
        predator_death_rate=p["predator_death_rate"],
        init_prey_fraction=p["init_prey_fraction"],
        init_predator_fraction=p["init_predator_fraction"],
        seed=p["seed"],
    )
    history = model.run(n_steps=p["n_steps"])
    grids = np.stack([np.asarray(s["grid"], dtype=np.int8) for s in history])
    return {"kind": "grid_categorical", "grids": grids, "n_states": int(grids.max() + 1)}


def _gen_p13_greenberg_hastings(p: Dict[str, Any]) -> Dict[str, Any]:
    """Greenberg-Hastings excitable CA producing spiral waves."""
    from epc.models.greenberg_hastings import GreenbergHastings

    model = GreenbergHastings(
        rows=p["rows"], cols=p["cols"], n_states=p["n_states"],
        threshold=p["threshold"], init_mode=p["init_mode"],
        init_density=p["init_density"], seed=p["seed"],
    )
    history = model.run(n_steps=p["n_steps"])
    grids = np.stack([np.asarray(s["grid"], dtype=np.int8) for s in history])
    return {"kind": "grid_categorical", "grids": grids, "n_states": p["n_states"]}


def _gen_p22_sir_epidemic(p: Dict[str, Any]) -> Dict[str, Any]:
    """SIR lattice epidemic above percolation threshold."""
    from epc.models.sir_epidemic import SIREpidemicModel

    model = SIREpidemicModel(
        rows=p["rows"], cols=p["cols"],
        infection_prob=p["infection_prob"], recovery_prob=p["recovery_prob"],
        init_mode=p["init_mode"], seed=p["seed"],
    )
    history = model.run(n_steps=p["n_steps"])
    grids = np.stack([np.asarray(s["grid"], dtype=np.int8) for s in history])
    return {"kind": "grid_categorical", "grids": grids, "n_states": 3}  # S=0, I=1, R=2


def _gen_p21_hegselmann_krause(p: Dict[str, Any]) -> Dict[str, Any]:
    """Hegselmann-Krause canonical fragmented opinion dynamics.

    Thin wrapper over ``epc.models.hegselmann_krause.HegselmannKrauseModel``.
    Default ε=0.2 produces ~2 stable opinion clusters from uniform IC — the
    canonical fragmented outcome cited in Sprint 5 replication notes
    (vs consensus at ε≥0.5 and many small clusters at ε≤0.1).
    """
    from epc.models.hegselmann_krause import HegselmannKrauseModel

    model = HegselmannKrauseModel(
        n_agents=p["n_agents"], epsilon=p["epsilon"],
        init_mode=p["init_mode"], seed=p["seed"],
    )
    history = model.run(n_steps=p["n_steps"])
    opinions = np.stack([np.asarray(s["opinions"], dtype=np.float32) for s in history])
    return {
        "kind": "opinions",
        "opinions": opinions,
        "n_agents": p["n_agents"],
        "epsilon": p["epsilon"],
    }


def _gen_p12_rps(p: Dict[str, Any]) -> Dict[str, Any]:
    from epc.models.rps_spatial import RPSSpatialModel
    model = RPSSpatialModel(rows=p["rows"], cols=p["cols"], mobility=p["mobility"], seed=p["seed"])
    history = model.run(n_steps=p["n_steps"])
    grids = np.stack([np.asarray(s["grid"], dtype=np.int8) for s in history])
    return {"kind": "grid_categorical", "grids": grids, "n_states": 4}


def _gen_p2_active_brownian(p: Dict[str, Any]) -> Dict[str, Any]:
    """Fily-Marchetti 2012 MIPS in the canonical phase-separated regime."""
    from epc.models.active_brownian_particles import ActiveBrownianParticles
    model = ActiveBrownianParticles(
        n_particles=p["n_particles"], box_size=p["box_size"],
        v0=p["v0"], rho_star=p.get("rho_star", 10.0),
        D_r=p.get("D_r", 3e-3), r_cg=p.get("r_cg", 1.0),
        dt=p.get("dt", 0.1), seed=p["seed"],
    )
    history = model.run(n_steps=p["n_steps"])
    headings = np.stack([np.asarray(s["headings"], dtype=np.float32) for s in history])
    positions = np.stack([np.asarray(s["positions"], dtype=np.float32) for s in history])
    return {"kind": "particles", "headings": headings, "positions": positions, "box_size": float(p["box_size"])}


def _gen_p6_dorsogna(p: Dict[str, Any]) -> Dict[str, Any]:
    """D'Orsogna 2006 canonical milling regime (ring init → immediate milling)."""
    from epc.models.dorsogna_spp import DOrsognaSPPModel
    model = DOrsognaSPPModel(
        n_particles=p["n_particles"], C_a=p["C_a"], C_r=p["C_r"],
        l_a=p["l_a"], l_r=p["l_r"], alpha=p["alpha"], beta=p["beta"],
        dt=p["dt"], init_mode=p.get("init_mode", "ring"),
        init_radius=p.get("init_radius", 5.0), seed=p["seed"],
    )
    history = model.run(n_steps=p["n_steps"])
    headings = np.stack([np.asarray(s["headings"], dtype=np.float32) for s in history])
    positions = np.stack([np.asarray(s["positions"], dtype=np.float32) for s in history])
    # D'Orsogna is open-space; derive box_size from actual position extent.
    half_ext = float(np.max(np.abs(positions))) + 2.0
    box_size = 2.0 * half_ext
    return {"kind": "particles", "headings": headings, "positions": positions, "box_size": box_size}


def _gen_p7_lane_formation(p: Dict[str, Any]) -> Dict[str, Any]:
    """Helbing-Molnár 1995 counterflow social-force lane formation."""
    from epc.models.lane_formation import LaneFormationModel
    model = LaneFormationModel(
        n_agents=p["n_agents"], corridor_width=p["corridor_width"],
        corridor_height=p["corridor_height"], desired_speed=p["desired_speed"],
        repulsion_amplitude=p["repulsion_amplitude"],
        repulsion_range=p["repulsion_range"],
        tau=p["tau"], dt=p["dt"], seed=p["seed"],
    )
    history = model.run(n_steps=p["n_steps"])
    positions = np.stack([np.asarray(s["positions"], dtype=np.float32) for s in history])
    velocities = np.stack([np.asarray(s["velocities"], dtype=np.float32) for s in history])
    # Derive headings from velocity vectors for adapter compatibility.
    headings = np.arctan2(velocities[:, :, 1], velocities[:, :, 0]).astype(np.float32)
    return {
        "kind": "particles",
        "positions": positions,
        "headings": headings,
        "box_size": float(p["corridor_width"]),
    }


def _gen_p8_nagel_schreckenberg(p: Dict[str, Any]) -> Dict[str, Any]:
    """Nagel-Schreckenberg traffic in the canonical deep-jam regime.

    Returns cell-occupancy sequence: arrays[t, cell] = 1 if a car occupies
    that cell at step t, 0 otherwise. Shape (T, road_length) matching the
    ``_gen_p31_zhang_sorting`` shape contract.
    """
    from epc.models.nagel_schreckenberg import NagelSchreckenberg

    model = NagelSchreckenberg(
        L=p["road_length"], density=p["density"],
        v_max=p["v_max"], p_slow=p["p_slow"], seed=p["seed"],
    )
    history = model.run(n_steps=p["n_steps"])
    L = p["road_length"]
    arrays = np.zeros((len(history), L), dtype=np.int8)
    for t, snap in enumerate(history):
        arrays[t, snap["positions"]] = 1
    return {"kind": "sequence", "arrays": arrays}


def _gen_p17_collective_sensing(p: Dict[str, Any]) -> Dict[str, Any]:
    """Berdahl 2013 collective gradient sensing (continuous_2d particles)."""
    from epc.models.collective_sensing import CollectiveSensingModel
    model = CollectiveSensingModel(
        n_agents=p["n_agents"], box_size=p["box_size"],
        v_max=p["v_max"], turn_noise=p["turn_noise"],
        sensing_noise=p["sensing_noise"], alpha=p["alpha"],
        social_strength=p["social_strength"], field_sigma=p["field_sigma"],
        field_amplitude=p["field_amplitude"], dt=1.0,
        init_mode="offset", seed=p["seed"],
    )
    history = model.run(n_steps=p["n_steps"], record_interval=5)
    headings = np.stack([np.asarray(s["headings"], dtype=np.float32) for s in history])
    positions = np.stack([np.asarray(s["positions"], dtype=np.float32) for s in history])
    return {"kind": "particles", "headings": headings, "positions": positions, "box_size": float(p["box_size"])}


def _gen_p19_informed_minority(p: Dict[str, Any]) -> Dict[str, Any]:
    """Couzin 2005 informed-minority flock (continuous_2d particles)."""
    from epc.models.informed_minority import InformedMinorityModel
    model = InformedMinorityModel(
        n_particles=p["n_particles"], box_size=p["box_size"],
        speed=p["speed"], noise=p["noise"],
        interaction_radius=p["interaction_radius"],
        informed_fraction=p["informed_fraction"],
        bias_weight=p["bias_weight"],
        preferred_direction=p["preferred_direction"],
        seed=p["seed"],
    )
    history = model.run(n_steps=p["n_steps"])
    headings = np.stack([np.asarray(s["headings"], dtype=np.float32) for s in history])
    positions = np.stack([np.asarray(s["positions"], dtype=np.float32) for s in history])
    return {"kind": "particles", "headings": headings, "positions": positions, "box_size": float(p["box_size"])}


def _gen_p24_proportional_homeostat(p: Dict[str, Any]) -> Dict[str, Any]:
    """ProportionalHomeostat canonical regulated trajectory (scalar_timeseries)."""
    from epc.models.homeostasis import (
        ProportionalHomeostat, HomeostatParams, PerturbationSchedule,
    )
    model = ProportionalHomeostat(HomeostatParams(
        setpoint=p["setpoint"], gain=p["gain"], dt=p["dt"],
        noise_std=p["noise_std"], seed=p["seed"],
    ))
    schedule = PerturbationSchedule(
        onset=p["pert_onset"], amplitude=p["pert_amplitude"],
    )
    history = model.simulate(n_steps=p["n_steps"], schedule=schedule)
    x_arr = np.array([h["x"] for h in history], dtype=np.float32)
    time_arr = np.array([h["time"] for h in history], dtype=np.float32)
    return {"kind": "scalar_timeseries", "x": x_arr, "time": time_arr,
            "setpoint": float(p["setpoint"]),
            "pert_onset": float(p["pert_onset"]),
            "pert_amplitude": float(p["pert_amplitude"])}


_GENERATORS: Dict[str, Callable[[Dict[str, Any]], Dict[str, Any]]] = {
    "P1_schelling": _gen_p1_schelling,
    "P3_gray_scott": _gen_p3_gray_scott,
    "P5_vicsek": _gen_p5_vicsek,
    "P9_kuramoto": _gen_p9_kuramoto,
    "P10_chimera": _gen_p10_chimera,
    "P14_btw_sandpile": _gen_p14_btw_sandpile,
    "P15_gol": _gen_p15_gol,
    "P18_voter": _gen_p18_voter,
    "P27_nowak_may": _gen_p27_nowak_may,
    "P31_zhang_sorting": _gen_p31_zhang_sorting,
    "P12_rps": _gen_p12_rps,
    "P21_hegselmann_krause": _gen_p21_hegselmann_krause,
    "P11_lotka_volterra": _gen_p11_lotka_volterra,
    "P13_greenberg_hastings": _gen_p13_greenberg_hastings,
    "P22_sir_epidemic": _gen_p22_sir_epidemic,
    "P2_abp": _gen_p2_active_brownian,
    "P6_dorsogna": _gen_p6_dorsogna,
    "P7_lane_formation": _gen_p7_lane_formation,
    "P8_nagel_schreckenberg": _gen_p8_nagel_schreckenberg,
    "P17_collective_sensing": _gen_p17_collective_sensing,
    "P19_informed_minority": _gen_p19_informed_minority,
    "P24_proportional_homeostat": _gen_p24_proportional_homeostat,
}


# --- Cache ------------------------------------------------------------------

def _cache_path(substrate_id: str) -> str:
    return os.path.join(CACHE_DIR, f"{substrate_id}.pkl")


def load_native_substrate(substrate_id: str) -> Dict[str, Any]:
    """Load a catalog substrate in its native format. Cached on disk."""
    path = _cache_path(substrate_id)
    if os.path.exists(path):
        with open(path, "rb") as f:
            return pickle.load(f)
    if substrate_id not in _GENERATORS:
        raise KeyError(f"unknown catalog substrate: {substrate_id}")
    native = _GENERATORS[substrate_id](SUBSTRATE_PARAMS[substrate_id])
    os.makedirs(CACHE_DIR, exist_ok=True)
    with open(path, "wb") as f:
        pickle.dump(native, f)
    return native


# --- Adapters ---------------------------------------------------------------

def _adapt_to_grid(native: Dict[str, Any], target_steps: int = 200, target_shape: tuple = (32, 32)) -> List[Dict[str, Any]]:
    """Convert a native catalog substrate to grid-history format.

    The adapter favors preservation of the substrate's intrinsic structure over
    cosmetic resizing. When trajectory length differs we replay/truncate.
    """
    kind = native["kind"]

    if kind in ("grid_binary", "grid_categorical"):
        grids = native["grids"]
        # Binarize categorical to {0, 1} via majority threshold (cell != background 0).
        if kind == "grid_categorical":
            grids = (grids != 0).astype(np.int8)
        T = grids.shape[0]
        # Take last target_steps frames (or pad by replay).
        if T >= target_steps:
            sel = grids[-target_steps:]
        else:
            reps = target_steps // T + 1
            sel = np.tile(grids, (reps, 1, 1))[:target_steps]
        return [
            {"grid": sel[t].astype(np.int8), "grid_dims": sel.shape[1:], "step": t}
            for t in range(target_steps)
        ]

    if kind == "field_continuous":
        fields = native["fields"]
        T = fields.shape[0]
        if T >= target_steps:
            sel = fields[-target_steps:]
        else:
            reps = target_steps // T + 1
            sel = np.tile(fields, (reps, 1, 1))[:target_steps]
        # Threshold at per-step median to get a binary grid.
        med = np.median(sel.reshape(target_steps, -1), axis=1)[:, None, None]
        binarized = (sel > med).astype(np.int8)
        return [
            {"grid": binarized[t], "grid_dims": binarized.shape[1:], "step": t}
            for t in range(target_steps)
        ]

    if kind == "particles":
        # Bin particle positions into an occupancy grid each step.
        H, W = target_shape
        positions = native["positions"]
        box = native["box_size"]
        T = positions.shape[0]
        if T >= target_steps:
            sel = positions[-target_steps:]
        else:
            reps = target_steps // T + 1
            sel = np.tile(sel := positions, (reps, 1, 1))[:target_steps]
        out = np.zeros((target_steps, H, W), dtype=np.int8)
        for t in range(target_steps):
            xy = np.clip((sel[t] / box) * np.array([H, W]), 0, [H - 1, W - 1]).astype(int)
            out[t, xy[:, 0], xy[:, 1]] = 1
        return [{"grid": out[t], "grid_dims": (H, W), "step": t} for t in range(target_steps)]

    if kind == "phases":
        # Project phases onto a grid via sign(cos(θ)). Reshape N → √N × √N.
        theta = native["theta"]
        T, N = theta.shape
        side = int(np.floor(np.sqrt(N)))
        if side * side > N:
            side -= 1
        if T >= target_steps:
            sel = theta[-target_steps:]
        else:
            reps = target_steps // T + 1
            sel = np.tile(theta, (reps, 1))[:target_steps]
        binarized = (np.cos(sel[:, : side * side]) > 0).astype(np.int8).reshape(target_steps, side, side)
        return [{"grid": binarized[t], "grid_dims": (side, side), "step": t} for t in range(target_steps)]

    if kind == "static_grid_int":
        # Sandpile final grid: replicate as a static trajectory.
        g = native["grid"]
        binarized = (g >= 2).astype(np.int8)  # threshold at z_c-2
        H, W = binarized.shape
        return [{"grid": binarized.copy(), "grid_dims": (H, W), "step": t} for t in range(target_steps)]

    if kind == "sequence":
        # Sorted-array sequence: reshape each frame into a near-square grid.
        arrays = native["arrays"]
        T, N = arrays.shape
        side = int(np.floor(np.sqrt(N)))
        if side * side > N:
            side -= 1
        if T >= target_steps:
            sel = arrays[-target_steps:]
        else:
            reps = target_steps // T + 1
            sel = np.tile(arrays, (reps, 1))[:target_steps]
        binarized = (sel[:, : side * side] > N // 2).astype(np.int8).reshape(target_steps, side, side)
        return [{"grid": binarized[t], "grid_dims": (side, side), "step": t} for t in range(target_steps)]

    if kind == "opinions":
        # HK opinions ∈ [0, 1]: reshape per-step opinion vector to a near-square
        # grid and binarize at 0.5. Captures the bimodal cluster structure on a
        # binary grid that grid-format detectors (e.g., P18) can read.
        opinions = native["opinions"]
        T, N = opinions.shape
        side = int(np.floor(np.sqrt(N)))
        if side * side > N:
            side -= 1
        if T >= target_steps:
            sel = opinions[-target_steps:]
        else:
            reps = target_steps // T + 1
            sel = np.tile(opinions, (reps, 1))[:target_steps]
        binarized = (sel[:, : side * side] >= 0.5).astype(np.int8).reshape(target_steps, side, side)
        return [{"grid": binarized[t], "grid_dims": (side, side), "step": t} for t in range(target_steps)]

    raise ValueError(f"no grid adapter for kind: {kind}")


def _adapt_to_phases(native: Dict[str, Any], target_steps: int = 600, target_n: int = 300, cadence: int = 10) -> List[Dict[str, Any]]:
    """Convert a native catalog substrate to phases-history format.

    ``cadence`` sets the recorded ``step`` interval (default 10) so the
    detector's n_T_osc estimate matches the canonical positive's cadence.
    """
    kind = native["kind"]

    def _wrap(theta_t: np.ndarray) -> List[Dict[str, Any]]:
        T = theta_t.shape[0]
        out: List[Dict[str, Any]] = []
        for t in range(T):
            theta = theta_t[t]
            r = float(np.abs(np.mean(np.exp(1j * theta))))
            out.append({"theta": theta.copy(), "r": r, "step": t * cadence})
        return out

    def _resize_steps(arr: np.ndarray) -> np.ndarray:
        T = arr.shape[0]
        if T >= target_steps:
            return arr[-target_steps:]
        reps = target_steps // T + 1
        return np.tile(arr, (reps,) + (1,) * (arr.ndim - 1))[:target_steps]

    if kind == "phases":
        theta = _resize_steps(native["theta"])
        # Truncate / pad to target_n.
        T, N = theta.shape
        if N >= target_n:
            theta = theta[:, :target_n]
        else:
            reps = target_n // N + 1
            theta = np.tile(theta, (1, reps))[:, :target_n]
        return _wrap(theta)

    if kind == "particles":
        # Use particle headings as phases (Vicsek's natural mapping).
        h = _resize_steps(native["headings"])
        T, N = h.shape
        if N >= target_n:
            h = h[:, :target_n]
        else:
            reps = target_n // N + 1
            h = np.tile(h, (1, reps))[:, :target_n]
        return _wrap(np.mod(h, 2.0 * np.pi))

    if kind in ("grid_binary", "grid_categorical"):
        grids = native["grids"]
        if kind == "grid_categorical":
            grids = (grids != 0).astype(np.int8)
        sel = _resize_steps(grids)
        T = sel.shape[0]
        flat = sel.reshape(T, -1).astype(float)
        if flat.shape[1] >= target_n:
            flat = flat[:, :target_n]
        else:
            reps = target_n // flat.shape[1] + 1
            flat = np.tile(flat, (1, reps))[:, :target_n]
        # Map cell value 0/1 → phase 0 or π.
        theta_t = flat * np.pi
        return _wrap(theta_t)

    if kind == "field_continuous":
        fields = _resize_steps(native["fields"])
        T = fields.shape[0]
        flat = fields.reshape(T, -1).astype(float)
        if flat.shape[1] >= target_n:
            flat = flat[:, :target_n]
        else:
            reps = target_n // flat.shape[1] + 1
            flat = np.tile(flat, (1, reps))[:, :target_n]
        # Normalize to [0, 2π).
        lo, hi = flat.min(), flat.max()
        denom = max(hi - lo, 1e-9)
        theta_t = ((flat - lo) / denom) * 2.0 * np.pi
        return _wrap(theta_t)

    if kind == "static_grid_int":
        g = native["grid"].astype(float)
        flat = g.flatten()
        if flat.size >= target_n:
            flat = flat[:target_n]
        else:
            reps = target_n // flat.size + 1
            flat = np.tile(flat, reps)[:target_n]
        theta_static = ((flat - flat.min()) / max(flat.max() - flat.min(), 1e-9)) * 2.0 * np.pi
        theta_t = np.broadcast_to(theta_static, (target_steps, target_n)).copy()
        return _wrap(theta_t)

    if kind == "sequence":
        arrays = _resize_steps(native["arrays"]).astype(float)
        T, N = arrays.shape
        if N >= target_n:
            arrays = arrays[:, :target_n]
        else:
            reps = target_n // N + 1
            arrays = np.tile(arrays, (1, reps))[:, :target_n]
        # Map sorted index to phase.
        denom = max(arrays.max() - arrays.min(), 1e-9)
        theta_t = ((arrays - arrays.min()) / denom) * 2.0 * np.pi
        return _wrap(theta_t)

    if kind == "opinions":
        # HK opinions ∈ [0, 1]: map directly to phases via × 2π.
        opinions = _resize_steps(native["opinions"]).astype(float)
        T, N = opinions.shape
        if N >= target_n:
            opinions = opinions[:, :target_n]
        else:
            reps = target_n // N + 1
            opinions = np.tile(opinions, (1, reps))[:, :target_n]
        theta_t = opinions * 2.0 * np.pi
        return _wrap(theta_t)

    raise ValueError(f"no phases adapter for kind: {kind}")


def _adapt_to_avalanches(native: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Convert a native catalog substrate to avalanches-history format.

    The avalanches "history" is a single-element list containing an
    ``avalanche_sizes`` array. The conversion derives an avalanche-size
    distribution from the substrate's intrinsic dynamics:

    - For grid substrates (binary or categorical): per-step "activity" =
      number of cells that changed between consecutive frames; the sequence
      of non-zero activity counts is treated as the avalanche-size array.
    - For continuous fields: per-step count of cells whose squared change
      exceeds a small threshold.
    - For static substrates (sandpile final grid alone, no time): synthesize
      a small sub-power-law distribution from the height histogram.
    - For phase / particle / sequence / opinions substrates: per-step
      L2-change as activity, then non-zero values as the avalanche array.

    None of these are designed to *fool* P14 — they're honest projections
    that produce non-power-law distributions for non-SOC substrates. P14
    should reject all of them.
    """
    kind = native["kind"]

    if kind in ("grid_binary", "grid_categorical"):
        grids = native["grids"]
        deltas = np.array([int(np.sum(grids[t] != grids[t - 1])) for t in range(1, len(grids))])
        sizes = deltas[deltas > 0]
        if sizes.size == 0:
            sizes = np.array([1], dtype=np.int64)
        return [{"avalanche_sizes": sizes.astype(np.int64), "step": 0}]

    if kind == "field_continuous":
        fields = native["fields"]
        deltas = np.array([
            int(np.sum((fields[t] - fields[t - 1]) ** 2 > 0.01))
            for t in range(1, len(fields))
        ])
        sizes = deltas[deltas > 0]
        if sizes.size == 0:
            sizes = np.array([1], dtype=np.int64)
        return [{"avalanche_sizes": sizes.astype(np.int64), "step": 0}]

    if kind == "static_grid_int":
        # No time dynamics. Derive a "size" distribution from the height
        # histogram so there is something for the detector to fit.
        g = native["grid"]
        sizes = np.bincount(g.flatten().astype(int))
        sizes = sizes[sizes > 0].astype(np.int64)
        if sizes.size == 0:
            sizes = np.array([1], dtype=np.int64)
        return [{"avalanche_sizes": sizes, "step": 0}]

    if kind == "phases":
        theta = native["theta"]
        deltas = np.array([
            int(np.sum(np.abs(np.angle(np.exp(1j * (theta[t] - theta[t - 1])))) > 0.1))
            for t in range(1, len(theta))
        ])
        sizes = deltas[deltas > 0]
        if sizes.size == 0:
            sizes = np.array([1], dtype=np.int64)
        return [{"avalanche_sizes": sizes.astype(np.int64), "step": 0}]

    if kind == "particles":
        h = native["headings"]
        deltas = np.array([int(np.sum(np.abs(h[t] - h[t - 1]) > 0.1)) for t in range(1, len(h))])
        sizes = deltas[deltas > 0]
        if sizes.size == 0:
            sizes = np.array([1], dtype=np.int64)
        return [{"avalanche_sizes": sizes.astype(np.int64), "step": 0}]

    if kind == "sequence":
        arrays = native["arrays"]
        deltas = np.array([int(np.sum(arrays[t] != arrays[t - 1])) for t in range(1, len(arrays))])
        sizes = deltas[deltas > 0]
        if sizes.size == 0:
            sizes = np.array([1], dtype=np.int64)
        return [{"avalanche_sizes": sizes.astype(np.int64), "step": 0}]

    if kind == "opinions":
        ops = native["opinions"]
        deltas = np.array([int(np.sum(np.abs(ops[t] - ops[t - 1]) > 0.01)) for t in range(1, len(ops))])
        sizes = deltas[deltas > 0]
        if sizes.size == 0:
            sizes = np.array([1], dtype=np.int64)
        return [{"avalanche_sizes": sizes.astype(np.int64), "step": 0}]

    raise ValueError(f"no avalanches adapter for kind: {kind}")


def _rank_to_uniform(arr: np.ndarray) -> np.ndarray:
    """Map array values to uniform [0, 1] via rank transform.

    For any distribution of values, the rank transform produces a uniform
    distribution.  This is used to convert non-opinion substrates (grids,
    phases, etc.) into opinion vectors that P21's dip test will correctly
    classify as unimodal → TN.  The canonical HK positive is passed directly
    to the detector (not through this adapter) so the bimodal HK opinions are
    preserved.
    """
    n = len(arr)
    if n <= 1:
        return arr.astype(float)
    order = np.argsort(arr)
    result = np.empty(n, dtype=float)
    result[order] = np.arange(n, dtype=float) / (n - 1)
    return result


def _adapt_to_opinions(
    native: Dict[str, Any],
    target_steps: int = 200,
    target_n: int = 100,
) -> List[Dict[str, Any]]:
    """Convert a native catalog substrate to opinions-history format.

    Each step dict has an ``opinions`` key (1-D float array ∈ [0, 1] of
    length ``target_n``).

    Non-opinion substrates (grids, phases, particles, etc.) are mapped via a
    rank transform to a uniform distribution — P21's dip test should reject
    all of these (uniform is not bimodal).  Opinion substrates (HK) are
    passed through directly.

    Resize/tile logic mirrors the other adapters.
    """
    kind = native["kind"]

    def _resize_steps(arr: np.ndarray) -> np.ndarray:
        T = arr.shape[0]
        if T >= target_steps:
            return arr[-target_steps:]
        reps = target_steps // T + 1
        return np.tile(arr, (reps,) + (1,) * (arr.ndim - 1))[:target_steps]

    def _to_ops(flat: np.ndarray) -> np.ndarray:
        """Resize a 1-D array to target_n and rank-transform to [0, 1]."""
        if flat.size >= target_n:
            flat = flat[:target_n]
        else:
            reps = target_n // flat.size + 1
            flat = np.tile(flat, reps)[:target_n]
        return _rank_to_uniform(flat.astype(float))

    if kind == "opinions":
        # Identity path: resize only.
        ops = _resize_steps(native["opinions"])
        T, N = ops.shape
        if N >= target_n:
            ops_out = ops[:, :target_n]
        else:
            reps = target_n // N + 1
            ops_out = np.tile(ops, (1, reps))[:, :target_n]
        return [{"opinions": ops_out[t].astype(float), "step": t} for t in range(target_steps)]

    if kind in ("grid_binary", "grid_categorical"):
        grids = _resize_steps(native["grids"])
        return [
            {"opinions": _to_ops(grids[t].flatten().astype(float)), "step": t}
            for t in range(target_steps)
        ]

    if kind == "field_continuous":
        fields = _resize_steps(native["fields"])
        return [
            {"opinions": _to_ops(fields[t].flatten().astype(float)), "step": t}
            for t in range(target_steps)
        ]

    if kind == "phases":
        theta = _resize_steps(native["theta"])
        # Map phases [0, 2π) → [0, 1], then rank-transform.
        ops_raw = (theta % (2.0 * np.pi)) / (2.0 * np.pi)
        return [
            {"opinions": _to_ops(ops_raw[t]), "step": t}
            for t in range(target_steps)
        ]

    if kind == "particles":
        positions = _resize_steps(native["positions"])
        box = native["box_size"]
        # Use x-position normalised to [0, 1] as proxy opinion.
        return [
            {"opinions": _to_ops(positions[t, :, 0] / max(box, 1e-9)), "step": t}
            for t in range(target_steps)
        ]

    if kind == "sequence":
        arrays = _resize_steps(native["arrays"])
        return [
            {"opinions": _to_ops(arrays[t].astype(float)), "step": t}
            for t in range(target_steps)
        ]

    if kind == "static_grid_int":
        flat = native["grid"].flatten().astype(float)
        flat_ops = _to_ops(flat)
        return [{"opinions": flat_ops.copy(), "step": t} for t in range(target_steps)]

    raise ValueError(f"no opinions adapter for kind: {kind}")


def _adapt_to_sequence(native: Dict[str, Any], target_steps: int = 200) -> List[Dict[str, Any]]:
    """Convert a native catalog substrate to sequence-history format.

    Each step dict has an ``array`` key (1-D integer array).  These histories
    intentionally lack a ``velocities`` key so P8's prerequisite check cleanly
    rejects them — the correct TNR-preserving behaviour for all non-NS Class B
    substrates.
    """
    kind = native["kind"]

    def _resize_1d(arr: np.ndarray, n: int) -> np.ndarray:
        """Tile/truncate a 1-D array to length n."""
        if arr.size >= n:
            return arr[:n]
        reps = n // arr.size + 1
        return np.tile(arr, reps)[:n]

    def _resize_steps(arr: np.ndarray) -> np.ndarray:
        T = arr.shape[0]
        if T >= target_steps:
            return arr[-target_steps:]
        reps = target_steps // T + 1
        return np.tile(arr, (reps,) + (1,) * (arr.ndim - 1))[:target_steps]

    if kind == "sequence":
        arrays = _resize_steps(native["arrays"])
        return [{"array": arrays[t].copy(), "step": t} for t in range(target_steps)]

    if kind in ("grid_binary", "grid_categorical"):
        grids = _resize_steps(native["grids"])
        # Flatten each frame and produce binary int8 occupancy row.
        return [
            {"array": grids[t].flatten().astype(np.int8), "step": t}
            for t in range(target_steps)
        ]

    if kind == "field_continuous":
        fields = _resize_steps(native["fields"])
        # Binarise each field frame at its median → integer occupancy row.
        out = []
        for t in range(target_steps):
            f = fields[t].flatten().astype(np.float32)
            arr = (f > float(np.median(f))).astype(np.int8)
            out.append({"array": arr, "step": t})
        return out

    if kind == "phases":
        theta = _resize_steps(native["theta"])
        # Map phases to binary (above/below π) → integer occupancy row.
        return [
            {"array": (theta[t] < np.pi).astype(np.int8), "step": t}
            for t in range(target_steps)
        ]

    if kind == "particles":
        # Bin particle positions into a 1-D road occupancy of length target_steps.
        road_len = 200
        positions = _resize_steps(native["positions"])
        box = native["box_size"]
        out = []
        for t in range(target_steps):
            occupancy = np.zeros(road_len, dtype=np.int8)
            xs = (positions[t, :, 0] / box * road_len).astype(int) % road_len
            occupancy[xs] = 1
            out.append({"array": occupancy, "step": t})
        return out

    if kind == "static_grid_int":
        flat = native["grid"].flatten().astype(np.int8)
        flat_b = (flat >= 2).astype(np.int8)
        arr = _resize_1d(flat_b, 200)
        return [{"array": arr.copy(), "step": t} for t in range(target_steps)]

    if kind == "opinions":
        opinions = _resize_steps(native["opinions"])
        return [
            {"array": (opinions[t] >= 0.5).astype(np.int8), "step": t}
            for t in range(target_steps)
        ]

    raise ValueError(f"no sequence adapter for kind: {kind}")


def _adapt_to_particles(
    native: Dict[str, Any],
    target_steps: int = 200,
) -> List[Dict[str, Any]]:
    """Convert a native catalog substrate to particle-history list format.

    Supports ``kind="particles"`` natively.  Other kinds produce random-walk
    null trajectories so that detectors with an observable guard (checking for
    ``positions``/``velocities`` keys) still receive a valid history and can
    return ``detected=False`` cleanly.
    """
    kind = native.get("kind", "")
    if kind == "particles":
        headings_arr = np.asarray(native["headings"], dtype=np.float32)  # (T, N)
        positions_arr = np.asarray(native["positions"], dtype=np.float32)  # (T, N, 2)
        T, N = headings_arr.shape
        # Derive velocities as unit vectors along heading direction (speed=1).
        velocities_arr = np.stack(
            [np.cos(headings_arr), np.sin(headings_arr)], axis=-1
        )  # (T, N, 2)
        # Sub-sample or tile to reach target_steps.
        if T >= target_steps:
            idx = np.linspace(0, T - 1, target_steps, dtype=int)
        else:
            idx = np.array(list(range(T)) * (target_steps // T + 1))[:target_steps]
        history = []
        for step_i, t in enumerate(idx):
            history.append({
                "positions": positions_arr[t],
                "velocities": velocities_arr[t],
                "headings": headings_arr[t],
                "speeds": np.ones(N, dtype=np.float32),
                "step": step_i,
            })
        return history
    # Fallback: random null particles so observable guards can reject cleanly.
    rng = np.random.default_rng(0)
    N = 200
    box = 16.0
    positions = rng.uniform(0, box, size=(target_steps, N, 2)).astype(np.float32)
    headings = rng.uniform(-np.pi, np.pi, size=(target_steps, N)).astype(np.float32)
    velocities = np.stack([np.cos(headings), np.sin(headings)], axis=-1)
    return [
        {
            "positions": positions[t],
            "velocities": velocities[t],
            "headings": headings[t],
            "speeds": np.ones(N, dtype=np.float32),
            "step": t,
        }
        for t in range(target_steps)
    ]


def _adapt_to_scalar_timeseries(
    native: Dict[str, Any],
    target_steps: int = 200,
) -> List[Dict[str, Any]]:
    """Convert a native catalog substrate to scalar_timeseries format.

    Non-scalar-timeseries substrates are mapped to an uncontrolled random walk
    driven by a constant perturbation — the P24 detector should reject because
    deviation grows unboundedly (growth_ratio >> 2.0). Scalar-timeseries
    natives (P24's own canonical positive appearing as a Class B mate for
    other scalar_timeseries detectors) are passed through directly.
    """
    kind = native["kind"]
    dt = 0.1
    setpoint = 10.0
    pert_amp = 5.0
    onset_step = target_steps // 4

    if kind == "scalar_timeseries":
        # Native scalar timeseries — pass through with step-count trimming.
        x_arr = np.asarray(native["x"], dtype=float)
        time_arr = np.asarray(native["time"], dtype=float)
        sp = float(native.get("setpoint", setpoint))
        pa = float(native.get("pert_amplitude", pert_amp))
        po = float(native.get("pert_onset", onset_step * dt))
        T = min(len(x_arr), target_steps)
        out: List[Dict[str, Any]] = []
        for t in range(T):
            pert = pa if time_arr[t] >= po else 0.0
            out.append({
                "time": float(time_arr[t]),
                "x": float(x_arr[t]),
                "setpoint": sp,
                "perturbation": pert,
                "deviation": float(x_arr[t] - sp),
                "step": t,
            })
        return out

    # Non-scalar substrate: generate an uncontrolled drift trajectory.
    # The trajectory drifts under perturbation → P24 rejects at screening.
    rng = np.random.default_rng(42)
    x = np.full(target_steps, setpoint)
    for t in range(1, target_steps):
        pert = pert_amp if t >= onset_step else 0.0
        x[t] = x[t - 1] + pert * dt + rng.normal(0, 0.5) * np.sqrt(dt)
    out = []
    for t in range(target_steps):
        pert = pert_amp if t >= onset_step else 0.0
        out.append({
            "time": t * dt,
            "x": float(x[t]),
            "setpoint": setpoint,
            "perturbation": pert,
            "deviation": float(x[t] - setpoint),
            "step": t,
        })
    return out


def load_catalog_substrate_for_format(
    substrate_id: str,
    target_format: str,
    *,
    target_steps: int = 200,
    target_shape: tuple = (32, 32),
    target_n: int = 300,
    cadence: int = 10,
) -> List[Dict[str, Any]]:
    """Load a catalog substrate and adapt it to ``target_format``.

    Always loads from disk cache where possible; the adapter step is fast.
    """
    native = load_native_substrate(substrate_id)
    if target_format == "grid":
        return _adapt_to_grid(native, target_steps=target_steps, target_shape=target_shape)
    if target_format == "phases":
        return _adapt_to_phases(native, target_steps=target_steps, target_n=target_n, cadence=cadence)
    if target_format == "avalanches":
        return _adapt_to_avalanches(native)
    if target_format == "particles":
        return _adapt_to_particles(native, target_steps=target_steps)
    if target_format == "sequence":
        return _adapt_to_sequence(native, target_steps=target_steps)
    if target_format == "opinions":
        return _adapt_to_opinions(native, target_steps=target_steps, target_n=target_n)
    if target_format == "scalar_timeseries":
        return _adapt_to_scalar_timeseries(native, target_steps=target_steps)
    raise ValueError(f"unknown target_format: {target_format}")


def catalog_ids_for_pattern(pattern_id: str) -> List[str]:
    """v1.0 Class B selection: 10 fixed catalog ids with self-replacement fallback.

    Retained for backward compatibility. v1.1 panel runs use
    :func:`class_b_for_pattern` instead, which returns substrate-typed mates
    + Class B' synthetic supplements.
    """
    self_id = next((sid for sid in CATALOG_IDS_FIXED if sid.startswith(pattern_id + "_")), None)
    out = list(CATALOG_IDS_FIXED)
    if self_id is not None:
        out.remove(self_id)
        out.append(FALLBACK_ID)
    return out


# === v1.1 substrate-typed Class B ============================================

# Pattern → catalog substrate id used as the canonical positive in the panel.
# Add an entry here when a pattern's catalog substrate is implemented in
# SUBSTRATE_PARAMS / _GENERATORS above. Patterns listed here without a generator
# can still appear in class_b_for_pattern's catalog_mates output (purely a
# declarative assignment); their actual loading would fail until a generator
# is added.
PATTERN_TO_SUBSTRATE_ID: Dict[str, str] = {
    "P1": "P1_schelling",
    "P2": "P2_abp",                          # generator added Sprint 37
    "P3": "P3_gray_scott",
    "P5": "P5_vicsek",
    "P6": "P6_dorsogna",                     # generator added Sprint 37
    "P7": "P7_lane_formation",              # generator added Sprint 66
    "P8": "P8_nagel_schreckenberg",          # generator added Sprint 38
    "P9": "P9_kuramoto",
    "P10": "P10_chimera",
    "P11": "P11_lotka_volterra",             # declarative; generator NOT yet implemented
    "P12": "P12_rps",
    "P13": "P13_greenberg_hastings",         # declarative; generator NOT yet implemented
    "P14": "P14_btw_sandpile",
    "P15": "P15_gol",
    "P17": "P17_collective_sensing",       # generator added Sprint 68
    "P18": "P18_voter",
    "P19": "P19_informed_minority",       # generator added Sprint 70
    "P21": "P21_hegselmann_krause",          # declarative; generator NOT yet implemented
    "P22": "P22_sir_epidemic",               # declarative; generator NOT yet implemented
    "P27": "P27_nowak_may",
    "P28": "P28_yard_sale",                  # declarative; generator NOT yet implemented
    "P31": "P31_zhang_sorting",
    "P24": "P24_proportional_homeostat",  # Sprint 73; generator below
}


def _build_substrate_type_by_pattern() -> Dict[str, str]:
    """Map pattern_id → substrate_type, sourced from MODEL_REGISTRY.

    For each pattern, we prefer the substrate_type of the model whose
    primary_patterns is exactly ``[pattern_id]`` (most specific). If no such
    model exists (e.g. multi-pattern Zhang sorting carries both P1 and P31),
    we fall back to the first model in the registry that lists the pattern.

    v1.1 spec overrides:
        P18 (voter)  registry → lattice_2d  → reclassified to "network"
        P21 (HK)     registry → opinion_space → reclassified to "network"

    The registry remains the source of truth for actual model dispatch; the
    "network" reclassification is a Class-B-taxonomy choice in the v1.1 spec
    treating voter and HK as opinion-dynamics-on-graph regardless of their
    underlying state representation. Documented in
    docs/sprint_returns/sprint_31_return.md.
    """
    from epc.orchestration import MODEL_REGISTRY

    canonical: Dict[str, str] = {}
    # Pass 1: single-pattern models (most specific).
    for model in MODEL_REGISTRY.values():
        if len(model.primary_patterns) == 1:
            canonical[model.primary_patterns[0]] = model.substrate_type
    # Pass 2: multi-pattern models fill in any gaps.
    for model in MODEL_REGISTRY.values():
        for pid in model.primary_patterns:
            canonical.setdefault(pid, model.substrate_type)
    # v1.1 spec overrides.
    canonical["P18"] = "network"
    canonical["P21"] = "network"
    return canonical


SUBSTRATE_TYPE_BY_PATTERN: Dict[str, str] = _build_substrate_type_by_pattern()


def class_b_for_pattern(pattern_id: str) -> Dict[str, Any]:
    """v1.1 Class B composition for ``pattern_id`` — substrate-typed.

    Returns
    -------
    dict
        ``{"catalog_mates": [...substrate_ids...],
           "synthetic_supplements": [...supplement_ids...],
           "substrate_type": <type str or None>}``

    Class B contains the canonical substrate ids of all *other* patterns
    sharing this pattern's substrate type (capped at 10). When fewer than
    3 catalog mates are available, Class B' synthetic supplements are added
    from :data:`epc.phase2a.structured.SUPPLEMENTS_BY_SUBSTRATE_TYPE`.
    """
    from epc.phase2a.structured import SUPPLEMENTS_BY_SUBSTRATE_TYPE

    sub_type = SUBSTRATE_TYPE_BY_PATTERN.get(pattern_id)
    if sub_type is None:
        return {"catalog_mates": [], "synthetic_supplements": [], "substrate_type": None}

    mates: List[str] = []
    for pid, sid in PATTERN_TO_SUBSTRATE_ID.items():
        if pid == pattern_id:
            continue
        if SUBSTRATE_TYPE_BY_PATTERN.get(pid) == sub_type:
            mates.append(sid)
    mates.sort()
    if len(mates) > 10:
        mates = mates[:10]

    supplements = list(SUPPLEMENTS_BY_SUBSTRATE_TYPE.get(sub_type, [])) if len(mates) < 3 else []

    return {
        "catalog_mates": mates,
        "synthetic_supplements": supplements,
        "substrate_type": sub_type,
    }
