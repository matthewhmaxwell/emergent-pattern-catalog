"""Throughput and census-size projection (pilot T7): measured cost per program -> affordable maximum length L.

  python -m census.speed <sample_run_dir> census/pilot/counts_g1_24.json [--workers 5] [--days 10.5]

Cost per program is measured on a random sample run through the VALIDATED pipeline (simulation + screen +
knock-out when something is emergent + textbook measures), by layer combination. Each candidate L is then costed as
sum over (length <= L, layers) of count x mean seconds for that layer combination, and compared with the budget
(DESIGN §5: G1 census ~50% of a ~3-week run => ~10.5 days x workers).
"""
import argparse, collections, json
from census.triage import load

ap = argparse.ArgumentParser(); ap.add_argument("run"); ap.add_argument("counts")
ap.add_argument("--workers", type=int, default=5); ap.add_argument("--days", type=float, default=10.5); a = ap.parse_args()
rows = [r for r in load(a.run) if "sim_s" in r]; counts = json.load(open(a.counts))
by = collections.defaultdict(list)
for r in rows: by[r["layers"]].append(r)
print(f"sample: {len(rows)} programs | errors {sum(1 for r in load(a.run) if r['status'] == 'ERROR')} | "
      f"fingerprint errors {sum(r.get('fp_errors', 0) for r in rows)}")
print(f"{'layers':7s} {'n':>4s} {'sec/program':>12s} {'flagged':>8s} {'knock-out run':>14s}")
mean = {}
for L, rs in sorted(by.items()):
    t = [r["sim_s"] + r["post_s"] for r in rs]; mean[L] = sum(t) / len(t)
    print(f"{L:7s} {len(rs):4d} {mean[L]:12.1f} {100 * sum(r['status'] == 'FLAGGED' for r in rs) / len(rs):7.0f}% "
          f"{100 * sum(bool(r.get('knockout_run')) for r in rs) / len(rs):13.0f}%")
allt = [r["sim_s"] + r["post_s"] for r in rows]; overall = sum(allt) / len(allt)
print(f"overall mean {overall:.1f} s/program; flagged {100 * sum(r['status'] == 'FLAGGED' for r in rows) / len(rows):.0f}%")
budget = a.days * 86400 * a.workers
print(f"\nbudget: {a.days} days x {a.workers} workers = {budget / 3600:,.0f} worker-hours")
print(f"{'L':>3s} {'programs':>10s} {'worker-hours':>13s} {'days on ' + str(a.workers) + ' workers':>20s}  fits?")
for Lmax in range(18, int(counts['max_bits']) + 1):
    tot = n = 0
    for key, c in counts["by_length_layers"].items():
        ln, lay = key.split(":")
        if int(ln) <= Lmax: n += c; tot += c * mean.get(lay, overall)
    print(f"{Lmax:3d} {n:10,d} {tot / 3600:13,.0f} {tot / 86400 / a.workers:20.1f}  {'yes' if tot <= budget else 'no'}")
