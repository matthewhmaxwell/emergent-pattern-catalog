"""Leave-one-behaviour-out sensitivity test (POST-REGISTRATION, EXPLORATORY; added 2026-10-26 after peer review).

Question: if a behaviour the census found had NOT been in the reference library (no library examples and no
textbook measure), would the pipeline have sent it to the literature check as an unnamed phenomenon of its own?

For each of the 9 behaviours named by the nearest-example namer (power-law avalanches are named by the catalog's SOC
detector, not the namer, and are reported separately):
  1. remove that behaviour's library classes from the namer of its substrate family, and remove its textbook
     checks from every seed's check dictionary (a truly unknown behaviour has no textbook measure);
  2. re-name every flagged view of that family with the frozen CensusNamer decision rule;
  3. re-run the frozen fixed-radius grouping (radius 7, census order) on the family's UNNAMED views;
  4. score the views originally named as the behaviour:
       -> UNNAMED in a group whose representative (shortest member, the only one the literature check sees) is itself
          an originally-held-out view: the behaviour REACHES the literature check;
       -> UNNAMED but in a group led by some other phenomenon: hidden behind another representative;
       -> named as another behaviour, or routed to review as another behaviour's atypical look: absorbed.
Frozen modules (namer, triage, filter, knockout) are imported, not changed. NearestNamer's per-class radius R is
cached here (it depends only on the library); equivalence with the frozen namer is asserted on a sample first.

  python -m census.loco census/runs/census_v1_g1 --lib census/reflib/census_v3/library.json
"""
import argparse, collections, json, os, random, time
import numpy as np
from census.namer import FAMILY, BEHAVIOUR, CensusNamer, NearestNamer, family_namers, behaviour, robust
from census.triage import load, mean_fp, cluster


class CachedNamer(NearestNamer):
    """NearestNamer.name with exclude_* unused and R cached per class. Identical decisions (asserted below)."""

    def _R(self, top):
        if not hasattr(self, "_rc"): self._rc = {}
        if top not in self._rc:
            idx = np.flatnonzero(self.cls == top); nd = []
            for j in idx:
                oo = idx[self.var[idx] != self.var[j]]
                if len(oo): nd.append(np.sqrt(((self.Z[oo] - self.Z[j]) ** 2).sum(1)).min())
            self._rc[top] = np.percentile(nd, self.RPCT) if nd else 0.0
        return self._rc[top]

    def name(self, fp, checks, exclude_variant=None, exclude_class=None):
        assert exclude_variant is None and exclude_class is None
        if len(self.cls) < self.K: return None
        z = self._z([fp])[0]; d = np.sqrt(((self.Z - z) ** 2).sum(1)); o = np.argsort(d)[:self.K]
        top, c = collections.Counter(self.cls[o]).most_common(1)[0]
        if top not in self.classes or c < self.VOTES: return None
        dc = d[self.cls == top].min(); dother = d[self.cls != top].min() if (self.cls != top).any() else np.inf
        if dc > self._R(top) or dother < self.MARGIN * dc: return None
        return top if checks.get(top, False) else None


def cached_family_namers(examples):
    out = {}
    for fam, nm in family_namers(examples).items():
        c = CachedNamer.__new__(CachedNamer); c.__dict__.update(nm.__dict__); out[fam] = c
    return out


def make_census_namer(examples):
    cn = CensusNamer.__new__(CensusNamer); cn.namers = cached_family_namers(examples); return cn


def strip(seeds, held):
    return [{**s, "checks": {c: ok for c, ok in (s.get("checks") or {}).items() if c not in held}} for s in seeds]


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("run"); ap.add_argument("--lib", required=True)
    ap.add_argument("--radius", type=float, default=7.0); ap.add_argument("--check", type=int, default=3000)
    a = ap.parse_args(); t0 = time.time()
    lib = json.load(open(a.lib))["examples"]; rows = load(a.run)
    views = [(r, v, e["seeds"]) for r in rows for v, e in (r.get("views") or {}).items()
             if e.get("flagged") and v != "C.aval"]
    print(f"{len(rows):,} programs, {len(views):,} flagged namer views ({time.time()-t0:.0f}s)", flush=True)

    frozen = CensusNamer(a.lib); fast = make_census_namer(lib)
    random.seed(0); sample = random.sample(views, min(a.check, len(views)))
    mism = sum(frozen.name_view(v, s) != fast.name_view(v, s) for _, v, s in sample)
    print(f"equivalence check: {mism} mismatches in {len(sample)} views", flush=True)
    assert mism == 0, "cached namer differs from the frozen namer"

    base = [fast.name_view(v, s) for _, v, s in views]
    print(f"baseline named ({time.time()-t0:.0f}s): {collections.Counter(b['status'] for b in base)}", flush=True)
    lib_fam = {e["class"]: FAMILY.get(e["view"], e["view"]) for e in lib}
    results = []
    for b in sorted({bn for bn, _ in BEHAVIOUR.values()}):
        held = {c for c, (bn, _) in BEHAVIOUR.items() if bn == b}
        fams = {lib_fam[c] for c in held if c in lib_fam}
        nm_b = make_census_namer([e for e in lib if e["class"] not in held])
        # every view the census credited to b: named as b (primary), or b shown by its textbook measure (also / atypical)
        target = [i for i, x in enumerate(base) if x["behaviour"] == b or b in x["shown"]]
        n_primary = sum(1 for i in target if base[i]["behaviour"] == b)
        fam_idx = [i for i, (_, v, _) in enumerate(views) if FAMILY.get(v, v) in fams]
        new = {i: nm_b.name_view(views[i][1], strip(views[i][2], held)) for i in fam_idx}
        out = collections.Counter()
        for i in target:
            st = new[i]["status"]
            out["unnamed" if st == "UNNAMED" else ("review_as_other" if st.startswith("KNOWN") else "named_as_other")] += 1
        # grouping of the family's unnamed views, exactly as triage does it (per family, census order)
        tset = set(target); led_groups = 0; in_led = 0; groups_total = 0; rep_lens = []
        for fam in fams:
            pool = [(views[i][0], views[i][1], mean_fp(views[i][2]), i) for i in fam_idx
                    if FAMILY.get(views[i][1], views[i][1]) == fam and new[i]["status"] == "UNNAMED"]
            if not pool: continue
            pool.sort(key=lambda x: (x[0]["len"], x[0]["idx"]))
            lab = cluster([(r, v, f) for r, v, f, _ in pool], a.radius)
            for c in sorted(set(lab)):
                mem = [pool[k] for k in np.flatnonzero(lab == c)]
                mem.sort(key=lambda x: (x[0]["len"], x[0]["idx"])); groups_total += 1
                if mem[0][3] in tset:
                    led_groups += 1; rep_lens.append(mem[0][0]["len"])
                    in_led += sum(1 for m in mem if m[3] in tset)
        row = {"behaviour": b, "held_classes": sorted(held), "families": sorted(fams), "views": len(target), "views_named_primary": n_primary,
               "unnamed": out["unnamed"], "review_as_other": out["review_as_other"], "named_as_other": out["named_as_other"],
               "groups_in_family": groups_total, "groups_led_by_behaviour": led_groups,
               "views_in_led_groups": in_led, "shortest_led_rep_bits": min(rep_lens) if rep_lens else None,
               "reaches_literature_check": led_groups > 0}
        results.append(row)
        print(json.dumps(row), f"({time.time()-t0:.0f}s)", flush=True)
    av = sum(1 for r in rows if (r.get("views") or {}).get("C.aval", {}).get("flagged"))
    out = {"run": a.run, "library": a.lib, "generated": time.strftime("%Y-%m-%d %H:%M"), "radius": a.radius,
           "equivalence_mismatches": mism, "equivalence_sample": len(sample), "behaviours": results,
           "avalanche_views_detector_named": av,
           "note": "post-registration exploratory sensitivity analysis; frozen census modules unchanged"}
    json.dump(out, open(os.path.join(a.run, "loco.json"), "w"), indent=1)
    print("DONE", os.path.join(a.run, "loco.json"), flush=True)


if __name__ == "__main__":
    main()
