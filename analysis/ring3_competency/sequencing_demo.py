"""#4 sequencing: visit station A BEFORE goal G (order matters; once you leave A a reactive agent oscillates, so
it needs memory). RNN+ES trains at import. Render: agent goes to A (turns green when visited), then to G."""
import sys, numpy as np
ROOT = "/home/matthewhmaxwell/emergent-pattern-catalog"
sys.path.insert(0, ROOT); sys.path.insert(0, ROOT + "/analysis/ring3_competency")
from gallery.mk_ring3_demo import ASSETS
import render_extra as RX
import rnn_es as RE                              # ES trains at import -> RE.th
th = RE.th; N = RE.H and RE.N
Wx, Wh, b, Wo, bo = RE.unpack(th)


def rollout(start, A, G, T=50):
    p = list(start); h = np.zeros(RE.H); done = 0; traj = [tuple(p)]; a_step = None
    for t in range(T):
        x = np.concatenate([RE.sgn(p, A), RE.sgn(p, G)])
        h = np.tanh(Wx @ x + Wh @ h + b); a = int(np.argmax(Wo @ h + bo))
        np_ = (p[0] + RE.DIRS[a][0], p[1] + RE.DIRS[a][1])
        if 0 <= np_[0] < N and 0 <= np_[1] < N: p = list(np_)
        traj.append(tuple(p))
        if done == 0 and tuple(p) == A: done = 1; a_step = len(traj) - 1
        if done == 1 and tuple(p) == G: return traj, a_step, True
    return traj, a_step, False


chosen = None
for k in range(400):
    r = np.random.default_rng(k); start, Ax, Gx = RE.cfg(r)
    traj, a_step, ok = rollout(start, Ax, Gx)
    if ok and a_step and 6 <= len(traj) <= 34:
        chosen = (start, Ax, Gx, traj, a_step); break
if chosen is None:
    r = np.random.default_rng(0); start, Ax, Gx = RE.cfg(r); traj, a_step, ok = rollout(start, Ax, Gx)
    chosen = (start, Ax, Gx, traj, a_step or 1)
start, Ax, Gx, traj, a_step = chosen


def frame(t):
    av = t >= a_step
    mk = [(float(Ax[0]), float(Ax[1]), "#38a169" if av else "#805ad5", "s", True, 130),
          (float(Gx[0]), float(Gx[1]), "#2f855a", "s", (t == len(traj) - 1), 130)]
    return {"markers": mk, "agents": [(float(traj[t][0]), float(traj[t][1]), "#2b6cb0")]}


leg = [("agent", "o", "#2b6cb0", "white"), ("A: visit first", "s", "#805ad5", "#805ad5"), ("goal G", "s", "none", "#2f855a")]
print("sequencing len", len(traj), "a_step", a_step,
      RX.render_markers([frame(t) for t in range(len(traj))], ASSETS + "/ring3_sequencing",
                        "Sequencing (#4): visit A first, then reach goal G", leg, (-.5, N - .5, -.5, N - .5)))
