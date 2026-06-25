"""Does the graph-structure lens earn its place? Two independent discriminators,
each with its own positive family and the SAME random-graph null:
  degree_cv  : scale-free / hub families HIGH, random LOW
  modularity : community (block) families HIGH, random LOW
ADMIT only if BOTH axes show a clean gap (structured > random)."""
import numpy as np
import networkx as nx
from analysis.blind_spot_probes import community_network, scale_free_network, null_random_graph
from epc.metrics.graph_structure import graph_structure


def H(G):
    return [{"adjacency": nx.to_numpy_array(G)}]


rows = []
# --- synthetic controls (networkx generators), several seeds each ---
for s in range(3):
    rows.append((f"ER_random[{s}]", "null",
                 graph_structure(H(nx.erdos_renyi_graph(90, 0.06, seed=s)))))
    rows.append((f"BA_scalefree[{s}]", "scale-free",
                 graph_structure(H(nx.barabasi_albert_graph(90, 2, seed=s)))))
    rows.append((f"WS_smallworld[{s}]", "small-world",
                 graph_structure(H(nx.watts_strogatz_graph(90, 6, 0.1, seed=s)))))
    sizes = [30, 30, 30]
    pin, pout = 0.18, 0.008
    P = [[pin if i == j else pout for j in range(3)] for i in range(3)]
    rows.append((f"SBM_community[{s}]", "community",
                 graph_structure(H(nx.stochastic_block_model(sizes, P, seed=s)))))
# --- EPC blind-spot probes (the real test) ---
for fn, tr in [(community_network, "community(EPC)"), (scale_free_network, "scale-free(EPC)"),
               (null_random_graph, "null(EPC)")]:
    p = fn(0)
    rows.append((fn.__name__, tr, graph_structure(p["history"])))

print(f"{'graph':<22}{'truth':<16}{'deg_cv':>8}{'modul':>8}{'clust':>8}{'assort':>8}{'<k>':>6}")
for nm, t, r in sorted(rows, key=lambda x: -((x[2] or {}).get('degree_cv', -1))):
    if r:
        print(f"{nm:<22}{t:<16}{r['degree_cv']:>8.3f}{r['modularity']:>8.3f}"
              f"{r['clustering']:>8.3f}{r['assortativity']:>8.3f}{r['mean_degree']:>6.1f}")
    else:
        print(f"{nm:<22}{t:<16}   (n/a)")

vals = [(t, r) for nm, t, r in rows if r]
randoms = [r['degree_cv'] for t, r in vals if 'null' in t]
hubs = [r['degree_cv'] for t, r in vals if 'scale-free' in t]
rand_mod = [r['modularity'] for t, r in vals if 'null' in t]
comm_mod = [r['modularity'] for t, r in vals if 'community' in t]
print(f"\ndegree_cv  hubs(scale-free): {[round(v,2) for v in sorted(hubs)]}  "
      f"random: {[round(v,2) for v in sorted(randoms)]}")
gap_cv = min(hubs) - max(randoms)
print(f"           gap = {gap_cv:+.3f}  (thr ~{(min(hubs)+max(randoms))/2:.2f})")
print(f"modularity communities:      {[round(v,2) for v in sorted(comm_mod)]}  "
      f"random: {[round(v,2) for v in sorted(rand_mod)]}")
gap_mod = min(comm_mod) - max(rand_mod)
print(f"           gap = {gap_mod:+.3f}  (thr ~{(min(comm_mod)+max(rand_mod))/2:.2f})")
ok = gap_cv > 0.1 and gap_mod > 0.05
print(f"\nVERDICT: {'ADMIT — both axes separate' if ok else 'PARTIAL/DEFER — review gaps'}")
