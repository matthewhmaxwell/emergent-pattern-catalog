"""Recompute the screen-to-flag counts (Tables 1 and 2 of the pipeline paper) from the compressed per-program records.

  python -m census.paper_tables            # reads census/runs/census_v1_g{1,2}/results_*.jsonl.xz

Definitions: a view PASSED THE SCREEN if it was emergent in >= 2 of 3 seeds (em_seeds >= 2); it was REMOVED BY THE
KNOCK-OUT if it passed the screen and was not flagged; "programs with interaction rules" are records whose
interaction_free field is false. G2's results_0..4 contain 1,206 ERROR records (the sandpile operator before
Deviation 1, no views); results_fix.jsonl holds their re-runs on the original seeds.
"""
import collections
import glob
import json
import lzma
import os

HERE = os.path.dirname(os.path.abspath(__file__))
FAMILY = {"F0": "fields", "F1": "fields", "A": "agents", "C": "lattice states", "C.phase": "lattice phases",
          "N": "network node states", "N.phase": "network phases", "C.transient": "lattice transients",
          "C.aval": "avalanches"}


def tables(run):
    c = collections.Counter(); fam = collections.Counter(); status = collections.Counter()
    for f in sorted(glob.glob(os.path.join(HERE, "runs", run, "results_*.jsonl.xz"))):
        with lzma.open(f, "rt") as fh:
            for line in fh:
                r = json.loads(line); status[r.get("status", "?")] += 1
                inter = not r.get("interaction_free", False)
                for vn, v in (r.get("views") or {}).items():
                    c["views"] += 1
                    if v.get("em_seeds", 0) >= 2:
                        c["screen"] += 1; fam[(FAMILY.get(vn, vn), "screen")] += 1
                        c["screen_inter"] += inter
                        if v.get("flagged"):
                            c["flagged"] += 1; fam[(FAMILY.get(vn, vn), "flagged")] += 1
                        else:
                            c["removed"] += 1; c["removed_inter"] += inter
    return status, c, fam


if __name__ == "__main__":
    for run in ("census_v1_g1", "census_v1_g2"):
        status, c, fam = tables(run)
        print(f"== {run}: records by status {dict(status)}")
        print(f"   views {c['views']:,}; passed screen {c['screen']:,}; removed by knock-out {c['removed']:,} "
              f"({c['removed'] / c['screen']:.0%}); removed among programs with interaction rules "
              f"{c['removed_inter']:,} of {c['screen_inter']:,} ({c['removed_inter'] / c['screen_inter']:.0%}); "
              f"flagged views {c['flagged']:,}")
        for name in dict.fromkeys(FAMILY.values()):
            print(f"   {name:22s} passed screen {fam[(name, 'screen')]:>8,}  flagged {fam[(name, 'flagged')]:>8,}")
