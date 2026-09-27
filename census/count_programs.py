"""Pilot T7: exact canonical program counts by length (G1, G2) -> census/pilot/counts_<g>.json."""
import sys, time, json, collections
from census import grammar_g1 as G1, grammar_g2 as G2
g, L = sys.argv[1], int(sys.argv[2]); M = G1 if g == "g1" else G2
t = time.time(); cnt = collections.Counter(); by_layers = collections.Counter()
for bits, p in M.enumerate_programs(L):
    cnt[len(bits)] += 1; by_layers[(len(bits), p.layers)] += 1
out = {"grammar": g, "max_bits": L, "by_length": dict(sorted(cnt.items())),
       "cumulative": {n: sum(v for k, v in cnt.items() if k <= n) for n in sorted(cnt)},
       "by_length_layers": {f"{k[0]}:{k[1]}": v for k, v in sorted(by_layers.items())}, "seconds": time.time() - t}
json.dump(out, open(f"census/pilot/counts_{g}_{L}.json", "w"), indent=1); print(json.dumps(out["cumulative"]))
