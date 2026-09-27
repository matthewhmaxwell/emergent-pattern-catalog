"""G1 — the high-level census grammar (pilot draft v0; frozen at registration).

A program is a HEADER (which layers exist, how many entity types each has, layer-level constants) followed by one or
more RULES. Each rule is one mechanism-level template with a few quantized parameters. Rules are applied in order,
once per time step. Every choice is written with a prefix-free code (census/coding.py), so each program has an exact
length in bits; equivalent programs (type relabelings) are counted once, at their shortest encoding.

Layers   C = lattice cells (Moore-8 neighbourhood, torus)      A = self-propelled agents (continuous torus)
         N = network nodes (random graph, mean degree 4)        F = 1-2 continuous fields (diffusing, same grid as C)
Allowed layer combinations: C, A, N, F, CA, CF, AF, CAF (a network is its own world).

Parameter grids are coarse-to-fine (level 1 is cheapest). They were fixed before any census data existed; values are
round log-spaced numbers, not tuned to make particular known models short.
"""
from dataclasses import dataclass, field
from itertools import permutations
from .coding import tb_encode, tb_decode, tb_all, hier_codes, hier_decode

COMBOS = ["C", "A", "N", "F", "CA", "CF", "AF", "CAF"]
KMAX = 4                                                     # entity types per discrete layer: 1..4

GRIDS = {                                                    # kind -> levels (coarse -> fine)
    "rate":   [[1.0, 0.1], [0.3, 0.03, 0.01, 0.003], [0.5, 0.2, 0.05, 0.02, 0.005, 0.002, 0.001, 0.0003]],
    "count":  [[1, 3], [2, 4, 5, 6], [7, 8]],                # neighbour-count thresholds
    "frac":   [[0.5, 0.3], [0.2, 0.4, 0.6, 0.7]],            # fraction thresholds
    "radius": [[1.0, 3.0], [0.5, 2.0, 5.0, 8.0]],            # agent interaction radius (lattice units)
    "speed":  [[1.0, 0.3], [0.0, 0.1, 0.5, 2.0]],            # agent speed per step
    "angle":  [[0.3, 1.0], [0.0, 0.1, 2.0, 3.14159]],        # agent heading noise amplitude (rad)
    "diff":   [[0.2, 0.05], [0.1, 0.02, 0.15, 0.01]],        # field diffusion coefficient (explicit Euler, <0.25)
    "level":  [[0.5, 0.2], [0.1, 0.8, 0.05, 0.3]],           # field thresholds
    "game":   [[(1.9, 0.0), (1.5, 0.5)],                     # (T, S) with R=1, P=0: weak PD, snowdrift
               [(1.5, -0.5), (0.5, -0.5), (0.5, 0.5), (1.3, 0.0)]],   # PD, stag hunt, harmony, weak PD
}

# ------------------------------------------------------------------------------------------------ rule templates
# name: (target layer, required layers, [(param, kind)], validity(prog) -> bool, implicit-type semantics per layer)
# param kinds: GRIDS keys, or "type:X" (a type of layer X), "any:X" (a type of X or 'any'), "field".
T = {}


def tmpl(name, target, req, params, valid=lambda p: True, sem=None):
    T[name] = dict(target=target, req=req, params=params, valid=valid, sem=sem or {})


kc = lambda p: p.k.get("C", 0); ka = lambda p: p.k.get("A", 0); kn = lambda p: p.k.get("N", 0)
# lattice cells
tmpl("C.COPY",   "C", "C", [("p", "rate")])                                   # copy a random neighbour (voter)
tmpl("C.MAJ",    "C", "C", [("p", "rate")], lambda p: kc(p) >= 2)             # adopt the neighbourhood majority
tmpl("C.CONTAG", "C", "C", [("a", "type:C"), ("b", "type:C"), ("th", "count"), ("p", "rate")],
     lambda p: kc(p) >= 2)                                                    # a -> b if >= th neighbours are b
tmpl("C.CYCLE",  "C", "C", [("th", "count"), ("p", "rate")], lambda p: kc(p) >= 2,
     {"C": "cyclic"})                                                         # t -> t+1 if >= th neighbours are t+1
tmpl("C.FLIP",   "C", "C", [("a", "type:C"), ("b", "type:C"), ("p", "rate")], lambda p: kc(p) >= 2)   # a -> b
tmpl("C.SWAP",   "C", "C", [("p", "rate")], lambda p: kc(p) >= 2)             # exchange with a random neighbour
tmpl("C.SCHELL", "C", "C", [("th", "frac")], lambda p: kc(p) >= 3, {"C": "empty0"})   # unhappy -> random empty
tmpl("C.SAND",   "C", "C", [("p", "rate")], lambda p: kc(p) == 4, {"C": "fixed"})     # BTW: heights 0-3, topple at 4
tmpl("C.GAME",   "C", "C", [("g", "game")], lambda p: kc(p) == 2, {"C": "fixed"})     # imitate best payoff
tmpl("C.KURA",   "C", "C", [("K", "rate")])                                   # phase oscillators, neighbour coupling
tmpl("C.EMIT",   "C", "CF", [("b", "type:C"), ("f", "field"), ("p", "rate")])  # cells of type b feed field f
tmpl("C.SENSE",  "C", "CF", [("a", "type:C"), ("b", "type:C"), ("f", "field"), ("th", "level")],
     lambda p: kc(p) >= 2)                                                    # a -> b where field f > th
# agents (header carries speed + heading noise; with no rules agents are independent random walkers)
tmpl("A.ALIGN",  "A", "A", [("r", "radius")])                                 # heading -> mean neighbour heading
tmpl("A.ATTRACT", "A", "A", [("b", "any:A"), ("r", "radius")])                # steer to centroid of b-neighbours
tmpl("A.REPEL",  "A", "A", [("b", "any:A"), ("r", "radius")])                 # steer away from it
tmpl("A.QUORUM", "A", "A", [("th", "count"), ("r", "radius")])                # >= th neighbours -> stop this step
tmpl("A.CONTAG", "A", "A", [("a", "type:A"), ("b", "type:A"), ("th", "count"), ("r", "radius")],
     lambda p: ka(p) >= 2)
tmpl("A.CYCLE",  "A", "A", [("th", "count"), ("r", "radius")], lambda p: ka(p) >= 2, {"A": "cyclic"})
tmpl("A.TURNCELL", "A", "CA", [], sem={"C": "parity"})                        # turn +90 on even cell, -90 on odd
tmpl("A.WRITECELL", "A", "CA", [], sem={"C": "cyclic"})                       # cell under agent -> next state
tmpl("A.DEPOSIT", "A", "AF", [("f", "field"), ("p", "rate")])                 # add to field under agent
tmpl("A.CLIMB",  "A", "AF", [("f", "field")])                                 # steer up the field gradient
# network nodes
tmpl("N.COPY",   "N", "N", [("p", "rate")])
tmpl("N.MAJ",    "N", "N", [("p", "rate")], lambda p: kn(p) >= 2)
tmpl("N.CONTAG", "N", "N", [("a", "type:N"), ("b", "type:N"), ("th", "count"), ("p", "rate")], lambda p: kn(p) >= 2)
tmpl("N.CYCLE",  "N", "N", [("th", "count"), ("p", "rate")], lambda p: kn(p) >= 2, {"N": "cyclic"})
tmpl("N.FLIP",   "N", "N", [("a", "type:N"), ("b", "type:N"), ("p", "rate")], lambda p: kn(p) >= 2)
tmpl("N.GAME",   "N", "N", [("g", "game")], lambda p: kn(p) == 2, {"N": "fixed"})
tmpl("N.KURA",   "N", "N", [("K", "rate")])
tmpl("N.REWIRE", "N", "N", [("p", "rate")], lambda p: kn(p) >= 2)            # discordant edge -> random same-type
# fields
tmpl("F.DECAY",  "F", "F", [("f", "field"), ("d", "rate")])
tmpl("F.FEED",   "F", "F", [("f", "field"), ("F", "rate")])                   # f += F (1 - f)
tmpl("F.AUTOCAT", "F", "F", [("f", "field"), ("g", "field"), ("p", "rate")],
     lambda p: p.nf == 2)                                                     # r = p f g^2 : f -= r, g += r


def _ok_types(name, prm):
    """type arguments that must differ."""
    d = dict(prm)
    if "a" in d and "b" in d and isinstance(d["a"], int) and isinstance(d["b"], int) and d["a"] == d["b"]: return False
    if "f" in d and "g" in d and d["f"] == d["g"]: return False
    return True


@dataclass(frozen=True)
class Rule:
    tmpl: str
    params: tuple = ()                                       # ((name, value), ...)


@dataclass(frozen=True)
class Program:
    layers: str
    k: dict = field(default_factory=dict, hash=False, compare=False)
    ktuple: tuple = ()                                       # (kC, kA, kN) for hashing
    nf: int = 0
    D: tuple = ()                                            # per-field diffusion
    agent: tuple = ()                                        # (speed, noise)
    rules: tuple = ()

    def __post_init__(self):
        object.__setattr__(self, "k", {L: v for L, v in zip("CAN", self.ktuple) if L in self.layers and v})


def make(layers, k=None, nf=0, D=(), agent=(), rules=()):
    k = k or {}
    return Program(layers, ktuple=tuple(k.get(L, 0) for L in "CAN"), nf=nf, D=tuple(D), agent=tuple(agent),
                   rules=tuple(rules))


def valid_templates(p):
    return [n for n, t in T.items() if all(L in p.layers for L in t["req"]) and t["valid"](p)]


def _param_options(kind, p):
    """[(value, codeword)] for a parameter kind in the context of program p."""
    if kind.startswith("type:"):
        k = p.k[kind[5:]]; return [(i, tb_encode(i, k)) for i in range(k)]
    if kind.startswith("any:"):
        k = p.k[kind[4:]]; opts = ["any"] + list(range(k)) if k >= 2 else ["any"]
        return [(v, tb_encode(i, len(opts))) for i, v in enumerate(opts)]
    if kind == "field":
        return [(i, tb_encode(i, p.nf)) for i in range(p.nf)]
    return hier_codes(GRIDS[kind])


# ------------------------------------------------------------------------------------------------ encode / decode
def _header_codes():
    """[(partial program kwargs, codeword)] for every header."""
    out = []
    for ci, combo in enumerate(COMBOS):
        c0 = tb_encode(ci, len(COMBOS))
        ks = [[]]
        for L in "CAN":
            if L in combo: ks = [kk + [(L, kv, tb_encode(kv - 1, KMAX))] for kk in ks for kv in range(1, KMAX + 1)]
        for kk in ks:
            base = dict(layers=combo, k={L: v for L, v, _ in kk}); bits = c0 + "".join(b for _, _, b in kk)
            variants = [(base, bits)]
            if "A" in combo:
                variants = [({**b, "agent": (s, e)}, bb + cs + ce) for b, bb in variants
                            for s, cs in hier_codes(GRIDS["speed"]) for e, ce in hier_codes(GRIDS["angle"])]
            if "F" in combo:
                nv = []
                for b, bb in variants:
                    for nf in (1, 2):
                        for Ds in _prod([hier_codes(GRIDS["diff"])] * nf):
                            nv.append(({**b, "nf": nf, "D": tuple(v for v, _ in Ds)},
                                       bb + tb_encode(nf - 1, 2) + "".join(c for _, c in Ds)))
                variants = nv
            out.extend(variants)
    return out


def _prod(lists):
    res = [[]]
    for l in lists: res = [r + [x] for r in res for x in l]
    return res


def encode(p):
    bits = tb_encode(COMBOS.index(p.layers), len(COMBOS))
    for L in "CAN":
        if L in p.layers: bits += tb_encode(p.k[L] - 1, KMAX)
    if "A" in p.layers:
        bits += dict((v, c) for v, c in hier_codes(GRIDS["speed"]))[p.agent[0]]
        bits += dict((v, c) for v, c in hier_codes(GRIDS["angle"]))[p.agent[1]]
    if "F" in p.layers:
        bits += tb_encode(p.nf - 1, 2) + "".join(dict(hier_codes(GRIDS["diff"]))[d] for d in p.D)
    vt = valid_templates(p)
    for i, r in enumerate(p.rules):
        if i > 0: bits += "1"
        bits += tb_encode(vt.index(r.tmpl), len(vt))
        prm = dict(r.params)
        for name, kind in T[r.tmpl]["params"]:
            bits += dict(_param_options(kind, p))[prm[name]]
    return bits + "0"


def decode(bits):
    pos = 0; ci, pos = tb_decode(bits, pos, len(COMBOS)); combo = COMBOS[ci]; k = {}
    for L in "CAN":
        if L in combo: v, pos = tb_decode(bits, pos, KMAX); k[L] = v + 1
    agent, nf, D = (), 0, ()
    if "A" in combo:
        s, pos = hier_decode(bits, pos, GRIDS["speed"]); e, pos = hier_decode(bits, pos, GRIDS["angle"]); agent = (s, e)
    if "F" in combo:
        n, pos = tb_decode(bits, pos, 2); nf = n + 1; D = []
        for _ in range(nf): d, pos = hier_decode(bits, pos, GRIDS["diff"]); D.append(d)
    p = make(combo, k, nf, D, agent); vt = valid_templates(p); rules = []
    while True:
        ti, pos = tb_decode(bits, pos, len(vt)); name = vt[ti]; prm = []
        for pname, kind in T[name]["params"]:
            opts = _param_options(kind, p)
            for v, c in opts:
                if bits.startswith(c, pos): prm.append((pname, v)); pos += len(c); break
        rules.append(Rule(name, tuple(prm)))
        more = bits[pos]; pos += 1
        if more == "0": break
    assert pos == len(bits), "trailing bits"
    return make(combo, k, nf, D, agent, rules)


# ------------------------------------------------------------------------------------------------ canonical form
def _allowed_perms(p, L):
    k = p.k[L]; sems = {T[r.tmpl]["sem"].get(L) for r in p.rules} - {None}
    out = []
    for pi in permutations(range(k)):
        if "fixed" in sems and pi != tuple(range(k)): continue
        if "empty0" in sems and pi[0] != 0: continue
        if "cyclic" in sems and any(pi[(t + 1) % k] != (pi[t] + 1) % k for t in range(k)): continue
        if "parity" in sems and any(pi[t] % 2 != t % 2 for t in range(k)): continue
        out.append(pi)
    return out


def _relabel(p, perms):
    rules = []
    for r in p.rules:
        prm = []
        for (name, v), (_, kind) in zip(r.params, T[r.tmpl]["params"]):
            L = kind.split(":")[1] if ":" in kind else None
            prm.append((name, perms[L][v] if (L in perms and isinstance(v, int)) else v))
        rules.append(Rule(r.tmpl, tuple(prm)))
    return make(p.layers, p.k, p.nf, p.D, p.agent, rules)


def canonical_bits(p):
    """shortest (then lexicographically smallest) encoding over all semantics-preserving type relabelings."""
    Ls = [L for L in "CAN" if L in p.layers]
    groups = _prod([[(L, pi) for pi in _allowed_perms(p, L)] for L in Ls])
    return min(((len(b), b) for b in (encode(_relabel(p, dict(g))) for g in groups)))[1]


def is_valid(p):
    return all(_ok_types(r.tmpl, r.params) for r in p.rules)


# ------------------------------------------------------------------------------------------------ enumeration
def enumerate_programs(max_bits, canonical=True):
    """yield (bits, Program) for every valid (canonical) program of length <= max_bits, grouped by header."""
    for hdr, hbits in _header_codes():
        if len(hbits) + 1 > max_bits: continue
        p0 = make(hdr["layers"], hdr["k"], hdr.get("nf", 0), hdr.get("D", ()), hdr.get("agent", ()))
        vt = valid_templates(p0)
        rule_opts = []                                          # [(Rule, codeword-without-continuation)]
        for ti, name in enumerate(vt):
            tc = tb_encode(ti, len(vt))
            for combo in _prod([_param_options(kind, p0) for _, kind in T[name]["params"]]):
                prm = tuple((pn, v) for (pn, _), (v, _) in zip(T[name]["params"], combo))
                if not _ok_types(name, prm): continue
                rule_opts.append((Rule(name, prm), tc + "".join(c for _, c in combo)))
        rule_opts.sort(key=lambda x: len(x[1]))
        yield from _rules_dfs(p0, hbits, [], rule_opts, max_bits, canonical)


def _rules_dfs(p0, bits, rules, opts, max_bits, canonical):
    for r, c in opts:
        nb = bits + ("1" if rules else "") + c
        if len(nb) + 1 > max_bits: break                      # opts sorted by length
        rs = rules + [r]; p = make(p0.layers, p0.k, p0.nf, p0.D, p0.agent, rs); full = nb + "0"
        if not canonical or canonical_bits(p) == full: yield full, p
        yield from _rules_dfs(p0, nb, rs, opts, max_bits, canonical)


def describe(p):
    hdr = p.layers + "".join(f" k{L}={p.k[L]}" for L in p.k)
    if p.agent: hdr += f" v={p.agent[0]} eta={p.agent[1]}"
    if p.nf: hdr += f" fields={p.nf} D={list(p.D)}"
    return hdr + " | " + " ; ".join(r.tmpl + "(" + ",".join(f"{n}={v}" for n, v in r.params) + ")" for r in p.rules)
