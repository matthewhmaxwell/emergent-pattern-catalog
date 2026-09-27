"""G2 — the low-level census grammar (pilot draft v0; frozen at registration). Exists to test GRAMMAR DEPENDENCE.

Same header and substrates as G1 (grammar_g1), but each rule is composed from finer parts instead of being one
template:  rule := target layer, SUBJECT (any | type a), CONDITION, ACTION, RATE.  Every G1 template is expressible
as a G2 rule (usually at a different length), and G2 also admits combinations G1 has no template for (e.g.
"copy a random neighbour only if at least 3 neighbours are type b"). The emergence-complexity ordering is compared
across G1 and G2 by rank correlation on phenomena found in both.

G2 rules compile to the same interpreter (census/sim.py) through `compile_rule`, which maps each part onto generic
operations implemented in sim.State (_g2_* methods).
"""
from dataclasses import dataclass
from .coding import tb_encode, tb_decode, hier_codes, hier_decode
from . import grammar_g1 as G1

GRIDS = G1.GRIDS
COMBOS, KMAX = G1.COMBOS, G1.KMAX

# conditions / actions per target layer: name -> (params, required layers, validity(prog))
kc = G1.kc; ka = G1.ka; kn = G1.kn
COND = {
    "C": {"ALWAYS": ([], "C", lambda p: True),
          "CNT_GE": ([("b", "type:C"), ("th", "count")], "C", lambda p: kc(p) >= 2),
          "CNT_LE": ([("b", "type:C"), ("th", "count")], "C", lambda p: kc(p) >= 2),
          "NEXT_GE": ([("th", "count")], "C", lambda p: kc(p) >= 2),
          "SAME_LT": ([("th", "frac")], "C", lambda p: kc(p) >= 2),
          "FIELD_GT": ([("f", "field"), ("th", "level")], "CF", lambda p: True),
          "AGENT_ON": ([], "CA", lambda p: True)},
    "A": {"ALWAYS": ([], "A", lambda p: True),
          "CNT_GE": ([("b", "any:A"), ("th", "count"), ("r", "radius")], "A", lambda p: True),
          "CNT_LE": ([("b", "any:A"), ("th", "count"), ("r", "radius")], "A", lambda p: True),
          "NEXT_GE": ([("th", "count"), ("r", "radius")], "A", lambda p: ka(p) >= 2),
          "CELL_EVEN": ([], "CA", lambda p: True),
          "CELL_ODD": ([], "CA", lambda p: True),
          "FIELD_GT": ([("f", "field"), ("th", "level")], "AF", lambda p: True)},
    "N": {"ALWAYS": ([], "N", lambda p: True),
          "CNT_GE": ([("b", "type:N"), ("th", "count")], "N", lambda p: kn(p) >= 2),
          "CNT_LE": ([("b", "type:N"), ("th", "count")], "N", lambda p: kn(p) >= 2),
          "NEXT_GE": ([("th", "count")], "N", lambda p: kn(p) >= 2),
          "DISCORD": ([], "N", lambda p: kn(p) >= 2)},                     # a random neighbour has another type
    "F": {"ALWAYS": ([], "F", lambda p: True)},
}
ACT = {
    "C": {"SET": ([("b", "type:C")], "C", lambda p: kc(p) >= 2),
          "SET_NEXT": ([], "C", lambda p: kc(p) >= 2),
          "COPY_RAND": ([], "C", lambda p: True),
          "COPY_MAJ": ([], "C", lambda p: kc(p) >= 2),
          "SWAP_RAND": ([], "C", lambda p: kc(p) >= 2),
          "MOVE_EMPTY": ([], "C", lambda p: kc(p) >= 3),
          "IMITATE_BEST": ([("g", "game")], "C", lambda p: kc(p) == 2),
          "PHASE_COUPLE": ([("K", "rate")], "C", lambda p: True),
          "ADD_GRAIN": ([], "C", lambda p: kc(p) == 4),
          "EMIT": ([("f", "field"), ("amt", "rate")], "CF", lambda p: True)},
    "A": {"HEAD_MEAN": ([("b", "any:A"), ("r", "radius")], "A", lambda p: True),     # align
          "HEAD_TO": ([("b", "any:A"), ("r", "radius")], "A", lambda p: True),       # toward centroid
          "HEAD_AWAY": ([("b", "any:A"), ("r", "radius")], "A", lambda p: True),
          "HEAD_GRAD": ([("f", "field")], "AF", lambda p: True),
          "TURN_L": ([], "A", lambda p: True), "TURN_R": ([], "A", lambda p: True),
          "SLOW": ([], "A", lambda p: True),
          "SET": ([("b", "type:A")], "A", lambda p: ka(p) >= 2),
          "SET_NEXT": ([], "A", lambda p: ka(p) >= 2),
          "CELL_NEXT": ([], "CA", lambda p: True),
          "DEPOSIT": ([("f", "field"), ("amt", "rate")], "AF", lambda p: True)},
    "N": {"SET": ([("b", "type:N")], "N", lambda p: kn(p) >= 2),
          "SET_NEXT": ([], "N", lambda p: kn(p) >= 2),
          "COPY_RAND": ([], "N", lambda p: True),
          "COPY_MAJ": ([], "N", lambda p: kn(p) >= 2),
          "IMITATE_BEST": ([("g", "game")], "N", lambda p: kn(p) == 2),
          "PHASE_COUPLE": ([("K", "rate")], "N", lambda p: True),
          "REWIRE_SAME": ([], "N", lambda p: kn(p) >= 2)},
    "F": {"SCALE": ([("f", "field"), ("d", "rate")], "F", lambda p: True),          # f *= (1 - d)
          "RELAX1": ([("f", "field"), ("F", "rate")], "F", lambda p: True),         # f += F (1 - f)
          "AUTOCAT": ([("f", "field"), ("g", "field"), ("p", "rate")], "F", lambda p: p.nf == 2)},
}
# implicit type semantics (for canonical relabeling), per part
SEM = {("C", "NEXT_GE"): {"C": "cyclic"}, ("C", "SET_NEXT"): {"C": "cyclic"}, ("C", "MOVE_EMPTY"): {"C": "empty0"},
       ("C", "IMITATE_BEST"): {"C": "fixed"}, ("C", "ADD_GRAIN"): {"C": "fixed"},
       ("A", "NEXT_GE"): {"A": "cyclic"}, ("A", "SET_NEXT"): {"A": "cyclic"},
       ("A", "CELL_EVEN"): {"C": "parity"}, ("A", "CELL_ODD"): {"C": "parity"}, ("A", "CELL_NEXT"): {"C": "cyclic"},
       ("N", "NEXT_GE"): {"N": "cyclic"}, ("N", "SET_NEXT"): {"N": "cyclic"}, ("N", "IMITATE_BEST"): {"N": "fixed"}}
TARGETS = "CANF"


@dataclass(frozen=True)
class Rule2:
    target: str
    subject: object          # "any" or type index
    cond: str
    cparams: tuple
    act: str
    aparams: tuple
    rate: float
    tmpl: str = "G2"         # so G1 helpers that look at r.tmpl do not break


def _ok(p, parts_req):
    return all(L in p.layers for L in parts_req)


def targets(p):
    return [t for t in TARGETS if t in p.layers]


def subjects(p, t):
    if t == "F": return ["any"]
    k = p.k.get(t, 1); return ["any"] + list(range(k)) if k >= 2 else ["any"]


def conds(p, t): return [n for n, (_, req, v) in COND[t].items() if _ok(p, req) and v(p)]


def acts(p, t): return [n for n, (_, req, v) in ACT[t].items() if _ok(p, req) and v(p)]


def _param_combos(specs, p, budget):
    """yield (values, code) for the parameter list `specs`, pruning any prefix longer than budget."""
    if not specs:
        yield (), ""; return
    (n, k), rest = specs[0], specs[1:]
    for v, c in G1._param_options(k, p):
        if len(c) > budget: continue
        for vs, cs in _param_combos(rest, p, budget - len(c)):
            yield ((n, v),) + vs, c + cs


def rule_options(p, budget=64):
    """[(Rule2, codeword)] for every valid rule (codeword <= budget bits) in the context of program header p."""
    out = []; ts = targets(p); rates = hier_codes(GRIDS["rate"])
    for ti, t in enumerate(ts):
        c_t = tb_encode(ti, len(ts)); sb = subjects(p, t); cs = conds(p, t); as_ = acts(p, t)
        for si, s in enumerate(sb):
            c_s = tb_encode(si, len(sb))
            for ci, c in enumerate(cs):
                c_c = c_t + c_s + tb_encode(ci, len(cs))
                if len(c_c) > budget: continue
                for cp, cpc in _param_combos(COND[t][c][0], p, budget - len(c_c)):
                    for ai, a in enumerate(as_):
                        c_a = c_c + cpc + tb_encode(ai, len(as_))
                        if len(c_a) > budget: continue
                        for ap, apc in _param_combos(ACT[t][a][0], p, budget - len(c_a)):
                            if not _valid_combo(t, s, c, cp, a, ap): continue
                            base = c_a + apc
                            for rv, rc in ([(1.0, "")] if t == "F" else rates):
                                if len(base) + len(rc) <= budget: out.append((Rule2(t, s, c, cp, a, ap, rv), base + rc))
    return out


def _valid_combo(t, s, c, cp, a, ap):
    d = dict(ap)
    if a == "SET" and isinstance(s, int) and d.get("b") == s: return False                  # set to own type = no-op
    if a == "AUTOCAT" and d.get("f") == d.get("g"): return False
    return True


def _sem(p):
    out = {}
    for r in p.rules:
        for key in ((r.target, r.cond), (r.target, r.act)):
            for L, v in SEM.get(key, {}).items(): out.setdefault(L, set()).add(v)
    return out


def enumerate_programs(max_bits, canonical=True):
    for hdr, hbits in G1._header_codes():
        if len(hbits) + 1 > max_bits: continue
        p0 = G1.make(hdr["layers"], hdr["k"], hdr.get("nf", 0), hdr.get("D", ()), hdr.get("agent", ()))
        opts = sorted(rule_options(p0, budget=max_bits - len(hbits) - 1), key=lambda x: len(x[1]))
        yield from _dfs(p0, hbits, [], opts, max_bits, canonical)


def _dfs(p0, bits, rules, opts, max_bits, canonical):
    for r, c in opts:
        nb = bits + ("1" if rules else "") + c
        if len(nb) + 1 > max_bits: break
        rs = rules + [r]; p = G1.make(p0.layers, p0.k, p0.nf, p0.D, p0.agent, rs); full = nb + "0"
        if not canonical or canonical_bits(p) == full: yield full, p
        yield from _dfs(p0, nb, rs, opts, max_bits, canonical)


def rule_code(p, r):
    """codeword of one G2 rule in the context of program header p (computed directly, no option table)."""
    ts = targets(p); sb = subjects(p, r.target); cs = conds(p, r.target); as_ = acts(p, r.target)
    code = tb_encode(ts.index(r.target), len(ts)) + tb_encode(sb.index(r.subject), len(sb))
    code += tb_encode(cs.index(r.cond), len(cs))
    for (n, k), (_, v) in zip(COND[r.target][r.cond][0], r.cparams): code += dict(G1._param_options(k, p))[v]
    code += tb_encode(as_.index(r.act), len(as_))
    for (n, k), (_, v) in zip(ACT[r.target][r.act][0], r.aparams): code += dict(G1._param_options(k, p))[v]
    if r.target != "F": code += dict(hier_codes(GRIDS["rate"]))[r.rate]
    return code


def encode(p):
    return G1.header_bits(p) + "1".join(rule_code(p, r) for r in p.rules) + "0"


def _perm_ok(pi, sems):
    k = len(pi)
    if "fixed" in sems and pi != tuple(range(k)): return False
    if "empty0" in sems and pi[0] != 0: return False
    if "cyclic" in sems and any(pi[(t + 1) % k] != (pi[t] + 1) % k for t in range(k)): return False
    if "parity" in sems and any(pi[t] % 2 != t % 2 for t in range(k)): return False
    return True


def _relabel(p, perms):
    def rl(L, v): return perms[L][v] if (L in perms and isinstance(v, int)) else v
    rules = []
    for r in p.rules:
        cp = tuple((n, rl(k.split(":")[1] if ":" in k else None, v))
                   for (n, v), (_, k) in zip(r.cparams, COND[r.target][r.cond][0]))
        ap = tuple((n, rl(k.split(":")[1] if ":" in k else None, v))
                   for (n, v), (_, k) in zip(r.aparams, ACT[r.target][r.act][0]))
        rules.append(Rule2(r.target, rl(r.target, r.subject), r.cond, cp, r.act, ap, r.rate))
    return G1.make(p.layers, p.k, p.nf, p.D, p.agent, rules)


def canonical_bits(p):
    from itertools import permutations
    sems = _sem(p); Ls = [L for L in "CAN" if L in p.layers]
    groups = G1._prod([[(L, pi) for pi in permutations(range(p.k[L])) if _perm_ok(pi, sems.get(L, set()))] for L in Ls])
    return min(((len(b), b) for b in (encode(_relabel(p, dict(g))) for g in groups)))[1]


def describe(p):
    hdr = G1.describe(G1.make(p.layers, p.k, p.nf, p.D, p.agent, [])).split(" | ")[0]
    return hdr + " | " + " ; ".join(
        f"{r.target}[{r.subject}] {r.cond}({','.join(f'{n}={v}' for n, v in r.cparams)}) -> "
        f"{r.act}({','.join(f'{n}={v}' for n, v in r.aparams)}) p={r.rate}" for r in p.rules)
