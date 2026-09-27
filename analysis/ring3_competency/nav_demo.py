"""#1 navigation: re-evolve the domain-randomized FSM navigator (import runs evolution -> bestrule), save the
rule, then render it re-routing around a barrier to the goal (competency = reach the goal a different way when
the straight path is blocked)."""
import sys, json, numpy as np
ROOT = "/home/matthewhmaxwell/emergent-pattern-catalog"
sys.path.insert(0, ROOT); sys.path.insert(0, ROOT + "/analysis/ring3_competency")
from gallery.mk_ring3_demo import ASSETS
import render_extra as RX
import evolve_navigation as EN                    # <- evolution runs here; EN.bestrule after
from stateful_funnel import run_fsm
N, RES, rule = EN.N, EN.RES, EN.bestrule
json.dump([list(t) for t in rule], open(ROOT + "/analysis/ring3_competency/nav_rule.json", "w"))

scenarios = [("barrier", EN.barrier()), ("serpentine", EN.SERP), ("trap", EN.trap())]
starts = [(0, 0), (12, 0), (0, N - 1), (12, 6)]
chosen = None
for name, w in scenarios:
    for s in starts:
        tr, ok = run_fsm(rule, s, RES, w, 250)
        if ok and len(tr) >= 8:
            chosen = (name, s, w, tr); break
    if chosen:
        break
if chosen is None:
    tr, ok = run_fsm(rule, (0, 0), RES, set(), 250); chosen = ("open", (0, 0), set(), tr)
name, s, w, tr = chosen
idx = np.linspace(0, len(tr) - 1, min(len(tr), 24)).astype(int); tr = [tr[i] for i in idx]


def frame(pos):
    mk = [(float(x), float(y), "#4a5568", "s", True, 95) for (x, y) in w] + \
         [(float(RES[0]), float(RES[1]), "#38a169", "s", True, 120)]
    return {"markers": mk, "agents": [(float(pos[0]), float(pos[1]), "#2b6cb0")]}


leg = [("agent", "o", "#2b6cb0", "white"), ("barrier", "s", "#4a5568", "#2d3748"), ("goal", "s", "#38a169", "#38a169")]
spr = RX.render_markers([frame(p) for p in tr], ASSETS + "/ring3_navigation",
                        "Navigation (#1): re-route around a barrier to the goal", leg, (-.5, N - .5, -.5, N - .5))
print("navigation scenario", name, "frames", len(tr), spr)
