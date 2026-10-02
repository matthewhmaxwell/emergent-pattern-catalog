"""Diagnostic: does the grouping of UNNAMED views keep its meaning when the census is ~500x larger than the rehearsal?

  python -m census.t_groupscale <run_dir> [<run_dir> ...] --lib census/reflib/census_v1/library.json

Ward linkage heights grow with cluster size (height = sqrt(2 nA nB / (nA + nB)) * centroid distance), so a fixed
Ward threshold splits the same blob into more and more classes as n grows. Compared here: Ward (threshold 12) vs
leader grouping with a fixed radius (shortest program first), on the real unnamed views and on k-fold jittered copies.
"""
import argparse, collections
import numpy as np
from census.namer import FAMILY, CensusNamer, robust
from census.triage import load, mean_fp


def collect(run, namer):
    un = collections.defaultdict(list)
    for r in load(run):
        for v, e in (r.get("views") or {}).items():
            if not e.get("flagged") or v == "C.aval": continue
            if namer.name_view(v, e["seeds"])["status"] == "UNNAMED": un[FAMILY.get(v, v)].append((r, v, mean_fp(e["seeds"])))
    return un


def zmat(items):
    keys = sorted(set().union(*[f for *_, f in items]))
    X = np.array([[f.get(k, 0.0) for k in keys] for *_, f in items], float); med, sc = robust(X)
    return np.clip((X - med) / sc, -8, 8)


def ward(Z, t=12.0):
    from scipy.cluster.hierarchy import linkage, fcluster
    return fcluster(linkage(Z, "ward"), t=t, criterion="distance") if len(Z) > 1 else np.array([1])


def leader(Z, radius):
    """rows already in census order (shortest first). Pass 1: a row founds a class unless a leader is within radius.
    Pass 2: every row joins its nearest leader."""
    L = [0]
    for i in range(1, len(Z)):
        if np.sqrt(((Z[L] - Z[i]) ** 2).sum(1)).min() > radius: L.append(i)
    C = Z[L]
    return np.array([int(np.sqrt(((C - z) ** 2).sum(1)).argmin()) for z in Z])


def main():
    from sklearn.metrics import adjusted_rand_score as ari
    ap = argparse.ArgumentParser(); ap.add_argument("runs", nargs="+"); ap.add_argument("--lib", required=True)
    a = ap.parse_args(); namer = CensusNamer(a.lib); radii = [4, 5, 6, 7, 8, 10, 12]; rng = np.random.default_rng(0)
    for run in a.runs:
        un = collect(run, namer); print(f"\n=== {run}: unnamed views by family {dict((k, len(v)) for k, v in un.items())}")
        tot = collections.Counter()
        for fam, items in sorted(un.items()):
            items = sorted(items, key=lambda x: (x[0]["len"], x[0]["idx"])); Z = zmat(items); w = ward(Z)
            line = f"{fam:16s} n={len(Z):4d} d={Z.shape[1]:2d} ward12={len(set(w)):3d} |"; tot["ward"] += len(set(w))
            for R in radii:
                l = leader(Z, R); tot[R] += len(set(l)); line += f" R{R}:{len(set(l))}({ari(w, l):.2f})"
            print(line)
            if len(Z) >= 20:                                      # scale test: k jittered copies of every view
                for k in (10, 40):
                    Zk = np.repeat(Z, k, 0) + rng.normal(0, 0.5, (len(Z) * k, Z.shape[1]))
                    print(f"{'':16s} x{k:<3d} n={len(Zk):5d}: ward12={len(set(ward(Zk))):4d} | "
                          + " ".join(f"R{R}:{len(set(leader(Zk, R)))}" for R in radii))
        print("TOTAL classes: ward12=%d | " % tot["ward"] + " ".join(f"R{R}:{tot[R]}" for R in radii))


if __name__ == "__main__":
    main()
