"""TEST ONLY (not a planted model, not in the plan's scope): agents steering toward the centroid of neighbours within
radius 2, used to exercise census/planted_pipeline.py end to end before the code freeze."""
import numpy as np
from scipy.spatial import cKDTree
from census.planted_models import base

CID = "toy:attract"; LAYERS = "A"; K = {"A": 1}; SPEED = 0.3
RULES = [("propel", False), ("attract", True), ("noise", False)]
HEADLINE = {}; STEP = "one update"; RULE_TEXT = "toy"


def simulate(seed0, steps=1000, S=3, size="census", keep=None):
    on = base.kept(keep, RULES); NA = base.SIZES[size]["NA"]; B = 20.0 * np.sqrt(NA / 400)
    P = np.zeros((S, base.n_frames(steps, "pos"), NA, 2), np.float32); H = np.zeros((S, P.shape[1], NA), np.float32)
    for s, r in enumerate(base.rngs(seed0, S)):
        x = r.uniform(0, B, (NA, 2)); h = r.uniform(0, 2 * np.pi, NA); fi = 0
        for t in range(steps + 1):
            if t % 2 == 0: P[s, fi] = x; H[s, fi] = np.mod(h, 2 * np.pi); fi += 1
            if t == steps: break
            if 2 in on:
                tr = cKDTree(x, boxsize=B); pairs = tr.query_pairs(2.0, output_type="ndarray")
                d = np.zeros((NA, 2)); c = np.zeros(NA)
                if len(pairs):
                    dv = x[pairs[:, 1]] - x[pairs[:, 0]]; dv -= B * np.round(dv / B)
                    np.add.at(d, pairs[:, 0], dv); np.add.at(d, pairs[:, 1], -dv); np.add.at(c, pairs.ravel(), 1)
                m = c > 0; h[m] = np.arctan2(d[m, 1], d[m, 0])
            h = h + r.uniform(-0.15, 0.15, NA); x = (x + SPEED * np.stack([np.cos(h), np.sin(h)], 1)) % B
            x[x >= B] = 0.0
    return {"pos": P, "head": H, "at": np.zeros(H.shape, np.int8), "meta": base.meta(steps, 64, B)}


def presence(out): return {}
