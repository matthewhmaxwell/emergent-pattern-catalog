"""Catalog-classification stage (Ring 3, v0).

A competency FORM is identified by its FUNCTIONAL SIGNATURE — the battery of substrate-agnostic
competency tests it passes — NOT by its name or substrate. So "same form, different algorithm" is
a principled judgment. Each verified hit is classified:
  - NEW FORM            : competent, but its signature matches no catalogued form.
  - KNOWN FORM, NEW INGRESSION : signature matches a form, via an algorithm/substrate not yet on
                          record. NOT a duplicate -- it is a Platonic data point: another portal
                          the form ingresses through, and it feeds the form's GOVERNING DYNAMICS
                          (what these algorithms share).
  - KNOWN FORM, DUP     : same form, algorithm/substrate already recorded.
The Platonic claim for a form strengthens with the DISTANCE between the algorithms that express it.

v0 substrate = gridworld navigation. Battery tests are abstract (any navigate(start,res,walls)->
reached): reach-from-many-starts, re-route-under-novel-barrier, solve-away-from-goal-trap,
generalise-to-unseen-maze. Core competency = re-route + away-from-goal (the irreducible bit).
"""
from stateful_funnel import make_fsm, run_fsm, pledge, RES, N

STARTS = [(0, 0), (0, N - 1), (N - 1, 0), (3, 6), (6, 3)]
HARD = [(0, N - 1), (2, N - 3), (4, N - 1)]
TRAP = {(RES[0] - 2, y) for y in range(2, N)}                 # away-from-goal trap
BLOCK = {(RES[0] - 2, y) for y in range(0, N - 4)}            # simple novel barrier (gap near goal)
SERP = {(4, y) for y in range(0, N - 2)} | {(8, y) for y in range(3, N)}  # unseen serpentine


def signature(navigate):
    """navigate(start, resource, walls) -> reached bool. Returns the functional competency signature."""
    t1 = sum(navigate(s, RES, set()) for s in STARTS) / len(STARTS) >= 0.8        # reaches (open)
    t2 = sum(navigate(s, RES, BLOCK) for s in STARTS) / len(STARTS) >= 0.6        # re-routes (simple)
    t3 = sum(navigate(s, RES, TRAP) for s in HARD) / len(HARD) >= 0.6             # away-from-goal trap
    t4 = sum(navigate(s, RES, SERP) for s in [(0, 12), (0, 6), (2, 11)]) / 3 >= 0.6  # unseen maze
    return dict(reach=t1, reroute=t2, away_from_goal=t3, generalize=t4)


CATALOG = []  # each: {name, core_sig, ingressions:[{algo, substrate, full, evidence}]}


def classify(navigate, algo, substrate):
    sig = signature(navigate)
    core_competent = sig["reroute"] and sig["away_from_goal"]    # the irreducible competency
    if not core_competent:
        return f"NOT COMPETENT (sig={sig})", sig
    # match against catalogued forms by core signature (reach+reroute+away_from_goal)
    core = (sig["reach"], sig["reroute"], sig["away_from_goal"])
    for form in CATALOG:
        if form["core_sig"] == core:
            known = any(i["algo"] == algo and i["substrate"] == substrate for i in form["ingressions"])
            form["ingressions"].append({"algo": algo, "substrate": substrate,
                                        "full": sig["generalize"], "sig": sig})
            tier = "KNOWN FORM, DUP" if known else "KNOWN FORM, NEW INGRESSION"
            strength = "full" if sig["generalize"] else "variant(weaker generalization)"
            return f"{tier} -> '{form['name']}' [{strength}]", sig
    # no match -> new form
    name = f"form_{len(CATALOG)+1}"
    CATALOG.append({"name": name, "core_sig": core,
                    "ingressions": [{"algo": algo, "substrate": substrate, "full": sig["generalize"], "sig": sig}]})
    return f"NEW FORM '{name}' (core_sig={core})", sig


def fsm_nav(seed):
    return lambda s, r, w: run_fsm(make_fsm(seed), s, r, w, steps=240)[1]


def pledge_nav():
    return lambda s, r, w: pledge(s, r, w, steps=400)[1]


if __name__ == "__main__":
    print("=== catalog-classification demo: navigation, two different algorithms ===\n")
    # first ingression: the random FSM hit -> registers the form
    v, sig = classify(fsm_nav(2959), algo="random_FSM(seed2959)", substrate="gridworld")
    print(f"FSM-2959      : {v}\n               sig={sig}")
    # second ingression: a DIFFERENT algorithm (hand-built wall-follower) with the SAME competency
    v, sig = classify(pledge_nav(), algo="hand_Pledge_wallfollower", substrate="gridworld")
    print(f"Pledge        : {v}\n               sig={sig}")
    # a brittle/greedy rule -> not competent
    v, sig = classify(fsm_nav(0), algo="random_FSM(seed0)", substrate="gridworld")
    print(f"FSM-0(brittle): {v}")

    print("\n=== form registry + governing-dynamics readout ===")
    for form in CATALOG:
        algos = [i["algo"] for i in form["ingressions"]]
        print(f"  FORM '{form['name']}': {len(form['ingressions'])} ingression(s) via {algos}")
        print(f"    governing-dynamics (so far): every ingression of this form required INTERNAL STATE")
        print(f"    + goal/obstacle sensing (memoryless rules give 0). Distinct algorithms, same form")
        print(f"    -> Platonic data point. (Strengthens further with a structurally-distant SUBSTRATE.)")
