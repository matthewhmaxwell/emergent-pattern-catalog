"""Plain-English rendering of G1 programs (for class cards, the digest and literature checks)."""
from .sim import W, NA, NN

LAYER = {"C": f"a {W}x{W} lattice of cells (wrap-around edges, 8 neighbours each)",
         "A": f"{NA} self-propelled agents moving in a wrap-around square (20x20 when agents are alone, {W}x{W} when they share the lattice or fields)",
         "N": f"a random network of {NN} nodes (average 4 links each)",
         "F": f"continuous chemical field(s) diffusing on a {W}x{W} wrap-around grid"}


def _p(p): return "every step" if p >= 1 else f"with probability {p} per step"


def _t(v, L="cell"): return "any type" if v == "any" else f"type {v}"


def rule(r):
    d = dict(r.params); t = r.tmpl
    txt = {
        "C.COPY": lambda: f"each cell copies the type of a random neighbour ({_p(d['p'])})",
        "C.MAJ": lambda: f"each cell adopts its neighbourhood's majority type ({_p(d['p'])})",
        "C.CONTAG": lambda: f"a type-{d['a']} cell becomes type {d['b']} if at least {d['th']} neighbours are type {d['b']} ({_p(d['p'])})",
        "C.CYCLE": lambda: f"a cell of type t becomes type t+1 (cyclically) if at least {d['th']} neighbours already are t+1 ({_p(d['p'])})",
        "C.FLIP": lambda: f"a type-{d['a']} cell spontaneously becomes type {d['b']} ({_p(d['p'])})",
        "C.SWAP": lambda: f"neighbouring cells swap contents ({_p(d['p'])})",
        "C.SCHELL": lambda: f"an occupied cell (type 0 = empty) whose same-type share among occupied neighbours is below {d['th']} moves to a random empty cell",
        "C.SAND": lambda: f"grains drop on random cells ({_p(d['p'])} per cell); a cell reaching 4 grains topples one grain to each of its 4 neighbours (grains fall off the edge)",
        "C.GAME": lambda: f"cells play a two-strategy game with their 8 neighbours and themselves (type 0 cooperates, 1 defects; payoffs R=1, P=0, T={d['g'][0]}, S={d['g'][1]}) and copy the best-scoring neighbour",
        "C.KURA": lambda: f"each cell is a phase oscillator (natural frequencies ~ N(0, 0.05)) pulled toward its neighbours' phases with strength {d['K']}",
        "C.EMIT": lambda: f"type-{d['b']} cells release {d['p']} of field {d['f']} per step",
        "C.SENSE": lambda: f"a type-{d['a']} cell becomes type {d['b']} where field {d['f']} exceeds {d['th']}",
        "A.ALIGN": lambda: f"each agent turns to the average heading of agents within distance {d['r']}",
        "A.ATTRACT": lambda: f"each agent turns toward the centre of {_t(d['b'])} agents within distance {d['r']}",
        "A.REPEL": lambda: f"each agent turns away from the centre of {_t(d['b'])} agents within distance {d['r']}",
        "A.QUORUM": lambda: f"an agent with at least {d['th']} others within distance {d['r']} slows to 10% speed",
        "A.CONTAG": lambda: f"a type-{d['a']} agent becomes type {d['b']} if at least {d['th']} type-{d['b']} agents are within distance {d['r']}",
        "A.CYCLE": lambda: f"an agent of type t becomes t+1 (cyclically) if at least {d['th']} agents of type t+1 are within distance {d['r']}",
        "A.TURNCELL": lambda: "each agent turns left on an even-type cell and right on an odd-type cell",
        "A.WRITECELL": lambda: "each agent advances the type of the cell under it by one (cyclically)",
        "A.DEPOSIT": lambda: f"each agent deposits {d['p']} of field {d['f']} where it stands",
        "A.CLIMB": lambda: f"each agent turns up the gradient of field {d['f']}",
        "N.COPY": lambda: f"each node copies the type of a random neighbour ({_p(d['p'])})",
        "N.MAJ": lambda: f"each node adopts its neighbours' majority type ({_p(d['p'])})",
        "N.CONTAG": lambda: f"a type-{d['a']} node becomes type {d['b']} if at least {d['th']} neighbours are type {d['b']} ({_p(d['p'])})",
        "N.CYCLE": lambda: f"a node of type t becomes t+1 (cyclically) if at least {d['th']} neighbours are t+1 ({_p(d['p'])})",
        "N.FLIP": lambda: f"a type-{d['a']} node spontaneously becomes type {d['b']} ({_p(d['p'])})",
        "N.GAME": lambda: f"nodes play a two-strategy game with their neighbours (T={d['g'][0]}, S={d['g'][1]}) and copy the best-scoring neighbour",
        "N.KURA": lambda: f"each node is a phase oscillator pulled toward its neighbours' phases with strength {d['K']}",
        "N.REWIRE": lambda: f"a node linked to a different-type neighbour cuts that link and links to a random same-type node ({_p(d['p'])})",
        "F.DECAY": lambda: f"field {d['f']} decays by a fraction {d['d']} per step",
        "F.FEED": lambda: f"field {d['f']} is replenished toward 1 at rate {d['F']}",
        "F.AUTOCAT": lambda: f"field {d['g']} grows autocatalytically by consuming field {d['f']} (rate {d['p']} x f x g^2)",
    }[t]
    return txt()


def describe(p):
    parts = [LAYER[L] for L in p.layers]
    s = "World: " + "; plus ".join(parts) + "."
    if p.k: s += " " + " ".join(f"{ {'C': 'Cells', 'A': 'Agents', 'N': 'Nodes'}[L] } have {k} type{'s' if k > 1 else ''} (random start)." for L, k in p.k.items())
    if p.agent: s += f" Agents move at speed {p.agent[0]} with heading noise {p.agent[1]} rad per step."
    if p.nf: s += f" {p.nf} field(s) with diffusion {list(p.D)}."
    return s + " Rules, applied in order each step: " + "; ".join(f"({i + 1}) {rule(r)}" for i, r in enumerate(p.rules)) + "."
