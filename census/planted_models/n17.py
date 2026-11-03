"""N17 (planted card): two-dimensional Vicsek model in its coexistence window. It is a census grammar-1 program as it stands
(PLANTED_PIPELINE_PLAN.md): `A kA=1 v=0.5 eta=2.0 | A.ALIGN(r=1.0)`, 16 bits, census index 2358. The census interpreter
adds heading noise uniform in [-eta/2, eta/2] = [-1, 1] rad, i.e. the paper's [-eta_V*pi, eta_V*pi] with eta_V = 1/pi =
0.318, inside the paper's window 0.28 to 0.35; speed 0.5 and alignment radius 1 match the card. Census size: 400 agents
in a 20 x 20 box (density 1, the paper's density); area4: 1,600 agents in 40 x 40 (set by census.verify._set_size).
The frozen interpreter, its seed (census.runner.seed_of(bits)) and its knock-out are used unchanged.
RULE_TEXT is the frozen card text with the box, number, noise value, start and run length of the census run; the
sentence about the stripes seen in the paper's long runs is replaced by the census run length.
"""
from census import grammar_g1 as G

CID = "planted:N17"
BITS = "0010011011000000"
LAYERS = "A"
K = {"A": 1}
SPEED = 0.5
RULES = [("ALIGN", True)]
HEADLINE = {"speed": 0.5, "radius": 1.0, "noise_rad_full_width": 2.0, "eta_vicsek": 0.318, "N": 400, "box": 20}
STEP = "one census step = one Vicsek update (discrete time), as in the census"
RULE_TEXT = ("World: 400 self-propelled point particles in a periodic 2D box (square of side 20), density 1 per unit area. "
             "Rules, discrete time steps: (1) each particle sets its heading to the average heading of all particles within "
             "distance 1 (itself included), plus a uniformly random angle in [-eta*pi, eta*pi], with noise eta = 0.318 (a "
             "random angle between -1 and 1 radian); (2) it moves 0.5 along its new heading. Start from random positions and "
             "random headings along the four lattice directions. Run 1,000 steps.")


def prog(keep=None):
    from census import knockout as KO
    p = G.decode(BITS)
    return KO.knockout_program(p) if keep is not None and 1 not in keep else p


def simulate(seed0, steps=1000, S=3, size="census", keep=None):
    """the census interpreter; size is set globally by census.verify._set_size (census settings by default)."""
    from census import sim
    return sim.run(prog(keep), seed0=seed0, steps=steps, S=S)


def fidelity(budget_s=7200):
    return {"check": "not testable within the budget (plan): the paper's lanes of band trains need boxes of side ~1,200",
            "pass": None}


def presence(out):
    """travelling bands: a density modulation along the mean heading. In a 20 x 20 box at most one band fits; report the
    share of power in the lowest non-zero mode of the density profile along the mean heading (late frames), and the
    polarization (a band state has intermediate polarization)."""
    import numpy as np
    res = []
    B = float(out["meta"]["WA"])
    for s in range(out["pos"].shape[0]):
        P = out["pos"][s, -50:]; H = out["head"][s, -50:]; shares, pols = [], []
        for f in range(len(P)):
            z = np.exp(1j * H[f]).mean(); pols.append(abs(z)); u = np.array([np.cos(np.angle(z)), np.sin(np.angle(z))])
            x = (P[f] @ u) % B; prof = np.histogram(x, bins=20, range=(0, B))[0].astype(float); prof -= prof.mean()
            pw = np.abs(np.fft.rfft(prof)) ** 2; shares.append(float(pw[1] / (pw[1:].sum() + 1e-12)))
        res.append({"band_mode_share": round(float(np.mean(shares)), 3), "polarization": round(float(np.mean(pols)), 3)})
    return {"per_seed": res, "note": "lanes of counter-propagating band trains cannot fit a 20 x 20 box"}
