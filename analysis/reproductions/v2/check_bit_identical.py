"""Bit-identity checks for the repro v2 group A fast loops, on small cases only.

1. Minority game: the frozen-commit model (git show 39688e3) and the current model give
   identical records under the default options (catalog runs unchanged).
2. Minority game: mg_fast.advance() vs simulate() in chunks, for all four
   (tie_break, payoff) combinations: attendance, winning side, final scores,
   history bits and rng state. Also: mu returned by advance() equals the index
   rebuilt from the concatenated winning-side series (incl. the m bits before).
3. Voter: voter_fast.sweep() vs VoterModel.step(): grid after every sweep, rng state.
4. Nowak-May: nm_fast.generation() vs NowakMayModel.step(): grid and f_C every generation.

Sizes and lengths are far too small to inform any pass rule.
Writes analysis/outputs/repro_v2/bit_identity_checks.json.
"""
from __future__ import annotations

import importlib.util
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common  # noqa: E402
import mg_fast  # noqa: E402
import nm_fast  # noqa: E402
import voter_fast  # noqa: E402

from epc.models.minority_game import MinorityGame, MinorityGameParams  # noqa: E402
from epc.models.nowak_may import NowakMayModel  # noqa: E402
from epc.models.voter import VoterModel  # noqa: E402


def _rng_state(rng):
    return repr(rng.bit_generator.state)


def check_mg_default_unchanged():
    src = subprocess.run(["git", "show", f"{common.FROZEN_PLAN_COMMIT}:epc/models/minority_game.py"],
                         cwd=common.REPO, capture_output=True, text=True, check=True).stdout
    tmp = Path(tempfile.mkdtemp()) / "mg_frozen.py"
    tmp.write_text(src)
    spec = importlib.util.spec_from_file_location("mg_frozen", tmp)
    old = importlib.util.module_from_spec(spec)
    sys.modules["mg_frozen"] = old
    spec.loader.exec_module(old)
    cases = []
    for (N, m, rounds, seed) in [(101, 6, 3000, s) for s in range(5)] + [(11, 3, 2000, 7), (31, 5, 2000, 3)]:
        h_old = old.MinorityGame(old.MinorityGameParams(n_agents=N, memory=m, n_rounds=rounds, seed=seed)).simulate()
        h_new = MinorityGame(MinorityGameParams(n_agents=N, memory=m, n_rounds=rounds, seed=seed)).simulate()
        ok = len(h_old) == len(h_new) and all(
            a['attendance'] == b['attendance'] and a['winning_side'] == b['winning_side']
            and np.array_equal(a['choices'], b['choices']) for a, b in zip(h_old, h_new))
        cases.append({"N": N, "m": m, "rounds": rounds, "seed": seed, "identical": bool(ok)})
    return {"name": "mg_default_options_identical_to_frozen_model", "cases": cases,
            "pass": all(c["identical"] for c in cases)}


def check_mg_fast():
    cases = []
    for tie in ("first", "random"):
        for payoff in ("sign", "linear"):
            for (N, m, n_chunks, seed, chunk_mode) in [
                    (11, 3, 300, 0, "P"), (31, 5, 80, 1, "P"), (101, 6, 40, 2, "P"),
                    (129, 7, 20, 3, "P"), (255, 8, 6, 4, "P"), (511, 9, 4, 5, "P"),
                    (33, 4, 120, 6, "mixed")]:
                P = 2 ** m
                ref = MinorityGame(MinorityGameParams(n_agents=N, memory=m, n_rounds=P, seed=seed,
                                                      tie_break=tie, payoff=payoff))
                fast = MinorityGame(MinorityGameParams(n_agents=N, memory=m, n_rounds=P, seed=seed,
                                                       tie_break=tie, payoff=payoff))
                init_bits = [int(b) for b in fast.history_bits]
                att_r, win_r, att_f, win_f, mu_f = [], [], [], [], []
                for k in range(n_chunks):
                    n = P if chunk_mode == "P" else (1, 7, P, 3 * P + 5)[k % 4]
                    ref.params.n_rounds = n
                    recs = ref.simulate()
                    att_r += [r['attendance'] for r in recs]
                    win_r += [r['winning_side'] for r in recs]
                    a, w, mu = mg_fast.advance(fast, n)
                    att_f += a.tolist(); win_f += w.tolist(); mu_f += mu.tolist()
                # rebuild mu from the concatenated winning-side series (incl. m bits before)
                series = init_bits + win_f
                mu_rebuilt = [mg_fast.history_index(series[t:t + m]) for t in range(len(win_f))]
                ok = {
                    "attendance": att_r == att_f,
                    "winning_side": win_r == win_f,
                    "scores": bool(np.array_equal(ref.scores, fast.scores)),
                    "history_bits": bool(np.array_equal(ref.history_bits, fast.history_bits)),
                    "rng_state": _rng_state(ref.rng) == _rng_state(fast.rng),
                    "mu_equals_rebuilt_from_winning_side": mu_f == mu_rebuilt,
                }
                cases.append({"tie_break": tie, "payoff": payoff, "N": N, "m": m, "P": P,
                              "rounds": len(att_f), "chunking": chunk_mode, "seed": seed, "checks": ok,
                              "identical": all(ok.values())})
    return {"name": "mg_fast_advance_vs_simulate", "cases": cases,
            "pass": all(c["identical"] for c in cases)}


def check_voter_fast():
    cases = []
    for nb in ("von_neumann", "moore"):
        for (R, C, sweeps, seed) in [(16, 16, 40, 0), (32, 32, 30, 1), (24, 40, 30, 2)]:
            ref = VoterModel(rows=R, cols=C, init_mode="random", neighborhood=nb, boundary="periodic", seed=seed)
            fast = VoterModel(rows=R, cols=C, init_mode="random", neighborhood=nb, boundary="periodic", seed=seed)
            ref.setup(); fast.setup()
            same = True
            for _ in range(sweeps):
                ref.step(); voter_fast.sweep(fast)
                if not np.array_equal(ref.grid, fast.grid):
                    same = False
                    break
            ok = same and _rng_state(ref.rng) == _rng_state(fast.rng) and ref._step_count == fast._step_count
            cases.append({"neighborhood": nb, "rows": R, "cols": C, "sweeps": sweeps, "seed": seed,
                          "identical": bool(ok)})
    return {"name": "voter_fast_sweep_vs_step", "cases": cases, "pass": all(c["identical"] for c in cases)}


def check_nm_fast():
    cases = []
    for b in (1.9, 1.8, 1.85, 1.6):
        for (R, C, gens, seed) in [(60, 60, 150, 0), (40, 72, 150, 1), (3, 5, 20, 2)]:
            ref = NowakMayModel(rows=R, cols=C, b=b, init_mode="random", init_coop_fraction=0.5,
                                boundary="periodic", seed=seed)
            fast = NowakMayModel(rows=R, cols=C, b=b, init_mode="random", init_coop_fraction=0.5,
                                 boundary="periodic", seed=seed)
            ref.setup(); fast.setup()
            same = True
            for _ in range(gens):
                st = ref.step(); n_c = nm_fast.generation(fast)
                if not (np.array_equal(ref.grid, fast.grid) and st["coop_fraction"] == n_c / (R * C)
                        and n_c == int(np.count_nonzero(fast.grid == 0))):
                    same = False
                    break
            cases.append({"b": b, "rows": R, "cols": C, "generations": gens, "seed": seed,
                          "identical": bool(same)})
    return {"name": "nm_fast_generation_vs_step", "cases": cases, "pass": all(c["identical"] for c in cases)}


def main():
    t0 = time.time()
    checks = [check_mg_default_unchanged(), check_mg_fast(), check_voter_fast(), check_nm_fast()]
    out = {"description": __doc__, "git": common.git_info(), "env": common.env_info(),
           "timestamp": common.now_iso(), "runtime_s": time.time() - t0,
           "all_pass": all(c["pass"] for c in checks), "checks": checks}
    common.write_json(common.OUT_DIR / "bit_identity_checks.json", out)
    for c in checks:
        print(c["name"], "PASS" if c["pass"] else "FAIL")
    print("all_pass", out["all_pass"], "runtime_s %.1f" % out["runtime_s"])
    return 0 if out["all_pass"] else 1


if __name__ == "__main__":
    sys.exit(main())
