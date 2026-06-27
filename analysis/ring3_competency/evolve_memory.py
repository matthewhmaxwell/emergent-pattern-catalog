"""Discovery loop for a SECOND competency face: MEMORY (cue-dependent action).

Orthogonal to navigation (temporal information-use, not spatial goal-reaching). The agent sees a
cue (0/1) at t=0, then a VARIABLE delay full of distractor inputs (cue gone), then a decision step
where it must output the cue. A memoryless reactive agent has no information at decision time ->
pinned at 50% chance BY CONSTRUCTION. Competency = it stores the cue in internal state, ignores
distractors, and reports it. Domain randomization over cue, delay length, and distractor sequence
stops timing/counting hacks. Held-out = LONGER delays + more distractors than ever trained on, plus
a cue-flip test (flip the cue -> the output must flip = proof it reads memory, not a fixed policy).
"""
import random

# inputs: 0=cue0, 1=cue1, 2=neutral, 3=decision-prompt, 4=distractorA, 5=distractorB
NIN, NSTATES = 6, 4
DISTRACT = [2, 4, 5]


def random_rule(r, nstates=NSTATES):
    return [(r.randrange(2), r.randrange(nstates)) for _ in range(NIN * nstates)]


def mutate(rule, r, rate=0.06):
    c = rule[:]
    for i in range(len(c)):
        if r.random() < rate:
            c[i] = (r.randrange(2), r.randrange(NSTATES))
    return c


def episode(rule, cue, delay, r, nstates=NSTATES):
    seq = [cue] + [r.choice(DISTRACT) for _ in range(delay)] + [3]   # cue, distractors, decision
    s = 0; out = 0
    for inp in seq:
        out, s = rule[inp * nstates + s]; s %= nstates
    return out == cue                                                 # out from the decision step


def fitness(rule, r, trials=24, dmin=1, dmax=8, nstates=NSTATES):
    return sum(episode(rule, r.randrange(2), r.randrange(dmin, dmax + 1), r, nstates)
               for _ in range(trials)) / trials


def heldout(rule, r, trials=200, dmin=12, dmax=25):                  # delays FAR longer than trained
    return sum(episode(rule, r.randrange(2), r.randrange(dmin, dmax + 1), r) for _ in range(trials)) / trials


def cue_flip_test(rule, r, trials=100):
    ok = 0
    for _ in range(trials):
        d = r.randrange(12, 25)
        st = r.getstate()
        a = episode(rule, 0, d, r); r.setstate(st)              # same distractor sequence, flipped cue
        b = episode(rule, 1, d, r)
        ok += (a and b)
    return ok / trials


# --- memoryless baseline: must be ~chance (proves the task requires memory) ---
rb = random.Random(7)
base = [fitness(random_rule(rb, nstates=1), rb, nstates=1) for _ in range(50)]
print(f"memoryless (1-state) baseline fitness: mean {sum(base)/len(base):.2f} (chance=0.50)")

# --- directed search + domain randomization ---
POP, MU, GENS = 120, 18, 70
r = random.Random(123)
pop = [random_rule(random.Random(i)) for i in range(POP)]
parents = pop[:MU]
for g in range(GENS):
    rg = random.Random(1000 + g)                               # fresh random trials each generation
    scored = sorted(((fitness(ru, rg), ru) for ru in pop), key=lambda x: -x[0])
    if g % 15 == 0 or g == GENS - 1:
        print(f"gen {g:>2}: best train fitness {scored[0][0]:.3f}")
    parents = [ru for _, ru in scored[:MU]]
    pop = parents + [mutate(r.choice(parents), r) for _ in range(POP - MU - 5)] + [random_rule(r) for _ in range(5)]

best = parents[0]
ho = heldout(best, random.Random(55))
flip = cue_flip_test(best, random.Random(77))
print(f"\nbest evolved rule:")
print(f"  HELD-OUT (delays 12-25, never trained; distractor-heavy): {ho:.2f}")
print(f"  cue-flip test (correct for BOTH cues at long delay): {flip:.2f}")
discovered = ho >= 0.9 and flip >= 0.9
print("VERDICT:", "DISCOVERED COMPETENCY: MEMORY (holds + reports a cue across long unseen delays; reads memory not a fixed policy)"
      if discovered else "no robust memory (over-fit short delays / timing hack)")
if discovered:
    print("\nCATALOG: form #2 = 'memory'. Functional signature {cue-dependent output, generalises to")
    print("  unseen long delays, cue-flip->output-flip} is DISTINCT from form #1 'navigation'.")
    print("  Governing dynamic: requires INTERNAL STATE USED AS STORAGE (memoryless = chance) -- a")
    print("  second form whose ingression needs state, toward the law 'competency needs memory'.")
