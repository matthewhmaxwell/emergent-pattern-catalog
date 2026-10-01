"""Naming for the census — the method that passed the validation gate (confirm3, 2026-09-28), plus BEHAVIOUR names.

Decision logic (unchanged from the gate; census.validate4 imports it from here):
  NearestNamer: a view is given library class c only if (1) all 5 nearest library examples are class c, (2) the
  nearest class-c example is within R_c (90th percentile of class-c nearest-other-variant distances), (3) the nearest
  example of any other class is >= 1.5x farther, (4) c has >= 3 parameter variants, and (5) the run passes class c's
  textbook measure (name-then-verify).

Display layer (owner decision 2026-09-28 / 2026-10-01): library classes that share one textbook measure are the same
BEHAVIOUR reached by different mechanisms, so names are reported at behaviour level — e.g. "voter coarsening" and
"majority-rule coarsening" are both "domain coarsening". The mechanism is not part of the name: it is the program's
own rule set (its ROUTE). The mapping below is presentation only; it changes no naming decision.
"""
import collections, json
import numpy as np

MIN_VARIANTS = 3

BEHAVIOUR = {  # library class -> (behaviour name, catalog entry or None)
    "voter coarsening": ("domain coarsening", "P18"),
    "majority-rule coarsening": ("domain coarsening", "P18"),
    "Schelling segregation": ("segregation with vacancies", "P1"),
    "cyclic-CA waves/spirals": ("cyclic waves and spirals", None),
    "spatial PD chaos (Nowak-May)": ("spatial game chaos", "P27"),
    "network synchronization": ("network synchronization", "P9"),
    "lattice phase locking": ("lattice phase locking", "P9"),
    "co-evolving network fragmentation": ("network fragmentation", "P34"),
    "network voter consensus": ("network consensus", "P18"),
    "flocking": ("flocking", "P5"),
}
FAMILY = {"C": "grid", "C.transient": "grid", "A": "agents", "N": "network", "C.phase": "lattice-phase",
          "N.phase": "network-phase", "F0": "field", "F1": "field", "C.aval": "avalanche"}


def behaviour(cls): return BEHAVIOUR.get(cls, (cls, None))[0]


def catalog(cls): return BEHAVIOUR.get(cls, (cls, None))[1]


def robust(X):
    """robust centre/scale with a floor: a feature that is near-constant in most rows (MAD ~ 0) is scaled by its
    standard deviation instead, so it cannot blow up into +/-8 z-units (round-7 degeneracy)."""
    med = np.median(X, 0); mad = np.median(np.abs(X - med), 0) * 1.4826
    scale = np.maximum(mad, 0.5 * X.std(0)); scale[scale < 1e-6] = 1.0
    return med, scale


class NearestNamer:
    K, VOTES, RPCT, MARGIN = 5, 5, 90, 1.5

    def __init__(self, rows, keys, med, mad):
        self.keys, self.med, self.mad = keys, med, mad
        self.Z = self._z([r["fp"] for r in rows]); self.cls = np.array([r["class"] for r in rows])
        self.var = np.array([f"{r['class']}#{r['variant']}" for r in rows])
        self.classes = sorted(c for c in set(self.cls) if len(set(self.var[self.cls == c])) >= MIN_VARIANTS)

    def _z(self, fps):
        X = np.array([[f.get(k, 0.0) for k in self.keys] for f in fps], float).reshape(len(fps), len(self.keys))
        return np.clip((X - self.med) / self.mad, -8, 8)

    def name(self, fp, checks, exclude_variant=None, exclude_class=None):
        z = self._z([fp])[0]; m = np.ones(len(self.cls), bool)
        if exclude_variant: m &= self.var != exclude_variant
        if exclude_class: m &= self.cls != exclude_class
        if m.sum() < self.K: return None
        Zr, cr, vr = self.Z[m], self.cls[m], self.var[m]
        d = np.sqrt(((Zr - z) ** 2).sum(1)); o = np.argsort(d)[:self.K]
        top, c = collections.Counter(cr[o]).most_common(1)[0]
        if top not in self.classes or top == exclude_class or c < self.VOTES: return None
        idx = np.flatnonzero(cr == top); nd = []
        for j in idx:
            oo = idx[vr[idx] != vr[j]]
            if len(oo): nd.append(np.sqrt(((Zr[oo] - Zr[j]) ** 2).sum(1)).min())
        R = np.percentile(nd, self.RPCT) if nd else 0.0
        dc = d[cr == top].min(); dother = d[cr != top].min() if (cr != top).any() else np.inf
        if dc > R or dother < self.MARGIN * dc: return None
        return top if checks.get(top, False) else None                       # name-then-verify


def family_namers(lib_examples):
    """{family: NearestNamer} from the textbook-verified library examples (same construction as the gate)."""
    by = collections.defaultdict(list)
    for e in lib_examples:
        if e["verified"]: by[FAMILY.get(e["view"], e["view"])].append(e)
    out = {}
    for fam, rows in by.items():
        keys = sorted(set().union(*[r["fp"] for r in rows]))
        X = np.array([[r["fp"].get(k, 0.0) for k in keys] for r in rows]); med, mad = robust(X)
        out[fam] = NearestNamer(rows, keys, med, mad)
    return out


class CensusNamer:
    """Names one census view from its per-seed fingerprints + per-seed textbook checks.

    primary  : library class named in >= 2 emergent seeds (each seed named independently by NearestNamer)
    also     : behaviours whose textbook measure passes in >= 2 emergent seeds, other than the primary's behaviour
    status   : NAMED | KNOWN-BEHAVIOUR-ATYPICAL-LOOK (-> review) | UNNAMED (-> literature check)
    """
    def __init__(self, library_path):
        self.namers = family_namers(json.load(open(library_path))["examples"])

    def name_view(self, view, seeds):
        em = [s for s in seeds if s.get("emergent") and s.get("fp")]
        nm = self.namers.get(FAMILY.get(view, view))
        votes = collections.Counter()
        for s in em:
            c = nm.name(s["fp"], s.get("checks", {})) if nm else None
            if c: votes[c] += 1
        primary = next((c for c, n in votes.most_common(1) if n >= 2), None)
        passed = collections.Counter(c for s in em for c, ok in s.get("checks", {}).items() if ok)
        shown = {behaviour(c) for c, n in passed.items() if n >= 2}
        pb = behaviour(primary) if primary else None
        also = sorted(shown - {pb})
        status = "NAMED" if primary else ("KNOWN-BEHAVIOUR-ATYPICAL-LOOK" if shown else "UNNAMED")
        return {"primary_class": primary, "behaviour": pb, "catalog": catalog(primary) if primary else None,
                "also": also, "shown": sorted(shown), "status": status}
