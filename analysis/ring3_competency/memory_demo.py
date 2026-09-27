"""#2 memory: a cue at t=0, a long delay of distractors (cue GONE from the input), then a decision step where the
agent must report the cue. Memoryless = chance; competency = hold the cue in internal state across distractors.

QA fix (2026-09-26): the old single-row timeline kept the cue cell visible the whole time, which reads as "the
cue is still shown." Now two rows make the mechanism explicit: the INPUT row shows the cue only at t=0 then
distractors (cue vanishes), while the MEMORY row shows the cue being HELD (a faded tint of the cue colour) across
the whole delay -> and the held value drives the correct report at the decision step."""
import sys, random
ROOT = "/home/matthewhmaxwell/emergent-pattern-catalog"
sys.path.insert(0, ROOT); sys.path.insert(0, ROOT + "/analysis/ring3_competency")
from gallery.mk_ring3_demo import ASSETS
import render_diagram as RD
import evolve_memory as EM                      # evolves on import -> EM.best

rule = EM.best
CUEC = {0: "#4363d8", 1: "#dd6b20"}            # cue colours (blue / orange)
FADE = {0: "#a3bffa", 1: "#fbd38d"}            # faded tint = "held in memory"
cue, delay = 1, 10
r = random.Random(3)
seq = [cue] + [r.choice(EM.DISTRACT) for _ in range(delay)] + [3]     # cue, distractors, decision-prompt(3)
s = 0; outs = []
for inp in seq:
    out, s = rule[inp * EM.NSTATES + s]; s %= EM.NSTATES; outs.append(out)

steps = []
for k in range(1, len(seq) + 1):
    fr = []
    for t in range(k):
        # INPUT row (row 1, top): cue only at t0, then distractors, then the decision prompt
        if t == 0:
            fr.append((t, 1, CUEC[cue], "cue", "white"))
        elif seq[t] == 3:
            fr.append((t, 1, "#d6bcfa", "?", "#44337a"))
        else:
            fr.append((t, 1, "#e2e8f0", "·", "#333"))
        # MEMORY row (row 0, bottom): cue held (faded) across the delay; report at the decision step
        if seq[t] == 3:
            ok = outs[-1] == cue
            fr.append((t, 0, "#38a169" if ok else "#e53e3e", "report %d" % outs[-1], "white"))
        else:
            fr.append((t, 0, FADE[cue], "hold", "#5b4b1f" if cue == 1 else "#2a3f7a"))
    steps.append(fr)

leg = [("cue (t0 only)", CUEC[cue]), ("distractor", "#e2e8f0"), ("held in memory", FADE[cue]), ("reported", "#38a169")]
print("memory cue", cue, "delay", delay, "out", outs[-1],
      RD.timeline(steps, "Memory (#2): cue vanishes, held in memory, then reported",
                  ASSETS + "/ring3_memory", leg))
