"""Shared plumbing for the planted-model plug-ins (PLANTED_PIPELINE_PLAN.md, commit 85acbe1).

A plug-in module census/planted_models/nXX.py simulates one planted model and returns arrays in the census interpreter's
output format (census.sim.run), so the frozen screen, fingerprint, knock-out decision, textbook measures, namer and
verification classifier run on them unchanged. Every module defines:

  CID        "planted:NXX"
  LAYERS     census world type: "A" particles, "N" graph, "F" fields, "C" lattice (phases ride on C or N, as in census)
  K          entity types per discrete layer, as census Program.k, e.g. {"A": 1}, {"N": 2}, {"C": 1}, {} for fields only
  SPEED      propulsion speed for the particle view's velocities (census: prog.agent[0]); 0.0 for passive particles
  RULES      [(name, is_interaction), ...] in the plan's numbering (rule 1 first)
  HEADLINE   dict of the parameters actually simulated (the plan's headline setting), with paper sources in comments
  STEP       one line: the model time in one census step, and the integrator sub-steps
  RULE_TEXT  plain-words rule text for the census card (the frozen planted card text, edited only for size/length)
  simulate(seed0, steps=1000, S=3, size="census", keep=None) -> out
      size: "census" or "area4" (SIZES below); keep: plan numbers of the INTERACTION rules kept (None = all of them,
      () = knock-out). One-way rules and the header are always kept.
  fidelity(budget_s=7200) -> dict   the plan's fixed check, the paper's own measure, at the largest size that fits
  presence(out) -> dict             the same measure on census-scale output of the real program (no census code)

Never call census.filter, census.fingerprint, census.knockout, census.namer or census.triage from a plug-in or while
developing one: pipeline outputs are not computed before the code freeze. check_out() below validates format only.
"""
import hashlib
import numpy as np

REC = {"C": 2, "pos": 2, "head": 2, "at": 2, "nt": 2, "Cph": 5, "Nph": 5, "adj": 10, "F": 10}   # = census.sim.REC
SIZES = {"census": dict(W=64, NA=400, NN=200), "area4": dict(W=128, NA=1600, NN=800)}            # = verify.CONDITIONS
BLOWUP = 50.0                                                                                    # population x start


def seed_of(cid):
    """stable seed per planted card (census.runner.seed_of hashes the program's bit string the same way)."""
    return int(hashlib.sha1(cid.encode()).hexdigest()[:8], 16)


def rngs(seed0, S):
    """one independent generator per seed."""
    return [np.random.default_rng([seed0, s, 20261103]) for s in range(S)]


def n_frames(steps, key):
    return steps // REC[key] + 1


def meta(steps_run, W, WA=None, unstable=False):
    return {"null_shuffle": False, "steps_run": int(steps_run), "frozen_at": None, "absorbed_at": None,
            "unstable": bool(unstable), "rec": dict(REC), "W": int(W), "WA": float(WA if WA is not None else W)}


def kept(keep, rules):
    """set of interaction-rule numbers (1-based, plan numbering) that are switched on."""
    inter = [i + 1 for i, (_, it) in enumerate(rules) if it]
    if keep is None: return set(inter)
    bad = set(keep) - set(inter)
    if bad: raise ValueError(f"keep names non-interaction or unknown rules: {sorted(bad)}")
    return set(keep)


class Rule:
    __slots__ = ("tmpl",)
    def __init__(self, name): self.tmpl = "PLANTED." + name


class PlugProg:
    """the attributes of a census Program that census.filter.views and census.reflib read."""
    def __init__(self, mod, keep=None):
        on = kept(keep, mod.RULES)
        self.layers = mod.LAYERS; self.k = dict(mod.K); self.agent = (float(getattr(mod, "SPEED", 0.0)), 0.0)
        self.nf = 0; self.D = ()
        self.rules = tuple(Rule(n) for i, (n, it) in enumerate(mod.RULES) if (not it) or (i + 1) in on)


def check_out(out, mod, steps, S, size="census"):
    """format-only validation of a plug-in's output against the census interpreter's conventions."""
    z = SIZES[size]; W = z["W"]; errs = []
    shapes = {"C": (S, n_frames(steps, "C"), W, W), "Cph": (S, n_frames(steps, "Cph"), W, W),
              "pos": (S, n_frames(steps, "pos"), None, 2), "head": (S, n_frames(steps, "head"), None),
              "at": (S, n_frames(steps, "at"), None), "nt": (S, n_frames(steps, "nt"), z["NN"]),
              "Nph": (S, n_frames(steps, "Nph"), z["NN"]), "adj0": (S, z["NN"], z["NN"])}
    for k, sh in shapes.items():
        if k not in out: continue
        a = out[k]
        if not hasattr(a, "shape") or len(a.shape) != len(sh) or any(x is not None and x != y for x, y in zip(sh, a.shape)):
            if not (out["meta"].get("unstable") and k in ("pos", "head", "at", "F")):
                errs.append(f"{k}: shape {getattr(a, 'shape', None)} != {sh}")
        if k in ("C", "nt", "at") and a.dtype != np.int8: errs.append(f"{k}: dtype {a.dtype} (want int8)")
        if k in ("Cph", "Nph", "pos", "head") and a.dtype != np.float32: errs.append(f"{k}: dtype {a.dtype} (want float32)")
        if k == "adj0" and a.dtype != bool: errs.append("adj0: dtype must be bool")
    if "F" in out:
        a = out["F"]
        if a.ndim != 5 or a.shape[0] != S or a.shape[3:] != (W, W): errs.append(f"F: shape {a.shape} != (S, frames, nf, {W}, {W})")
        elif not out["meta"].get("unstable") and a.shape[1] != n_frames(steps, "F"): errs.append(f"F: {a.shape[1]} frames")
        if a.dtype != np.float32: errs.append(f"F: dtype {a.dtype} (want float32)")
    if "pos" in out:
        B = out["meta"]["WA"]
        if np.nanmin(out["pos"]) < 0 or np.nanmax(out["pos"]) >= B: errs.append("pos outside [0, box)")
    if ("C" in out) != ("C" in mod.LAYERS): errs.append("C layer mismatch")
    for key in ("rec", "W", "WA", "steps_run", "unstable", "frozen_at", "absorbed_at"):
        if key not in out.get("meta", {}): errs.append(f"meta.{key} missing")
    return errs
