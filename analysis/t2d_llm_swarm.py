"""T2d — external demonstration: a real multi-agent LLM swarm.

N independent LLM agents move on a 2D torus. Each agent, every round, is given
ONLY its local egocentric neighborhood (bearings, distances, and relative
headings of nearby agents within a sensing radius) and asked to choose a turn.
No agent ever sees the global state, and the *only* coordination is whatever the
agents do in response to the local rule stated in their prompt. Any global order
is therefore genuinely EMERGENT from local LLM decisions.

Three conditions, run blind through the catalog battery (T2d):
  align      "steer to match neighbors' heading"  -> expect flocking (P5)
  segregate  "move toward your own type, away from others" (Schelling-like,
             agents carry a binary type) -> expect aggregation/segregation (P1)
  random     "ignore neighbors, turn randomly" (control) -> expect NO-EMERGENCE

Output: a trace JSON {condition, params, history:[{positions, velocities[, types]}]}
consumed by analysis/t2d_profile.py (which runs profile_observation).

Pure stdlib + anthropic (no numpy) so it runs in the orchestrator venv. Reads
ANTHROPIC_API_KEY from the environment. Uses Haiku (cheap) by default.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import random
import re
from concurrent.futures import ThreadPoolExecutor

import anthropic

TWO_PI = 2.0 * math.pi


def wrap_angle(a: float) -> float:
    """Wrap to (-pi, pi]."""
    return (a + math.pi) % TWO_PI - math.pi


def torus_delta(a: float, b: float, L: float) -> float:
    """Minimum-image displacement b-a on a length-L torus."""
    d = (b - a) % L
    if d > L / 2:
        d -= L
    return d


SYSTEM = ("You are one agent in a swarm. You see only your local surroundings. Decide your "
          "turn this step in degrees, in [-60, 60] (positive = left, negative = right, "
          "0 = straight). Do any arithmetic SILENTLY and reply with ONLY the final integer "
          "(for example: -12). Do not show steps or write any words.")

RULES = {
    "align": ("RULE: match the group's direction. Look at the 'heading' value of each "
              "neighbor (their direction relative to yours, in degrees). Compute the "
              "AVERAGE of those heading values and turn by that amount so you end up "
              "moving the same way as the group. Output that average as a single integer, "
              "clipped to [-60, 60]. If you see no neighbors, output 0."),
    "segregate": ("RULE: you belong to a TYPE. Steer toward the middle of your OWN kind: "
                  "look at the 'bearing' values of the neighbors marked SAME and turn by "
                  "their AVERAGE so you head into your own group. Ignore DIFF neighbors. "
                  "Output that average of the SAME bearings as a single integer in [-60, 60]. "
                  "If you see no SAME neighbors, output 0."),
    "random": ("RULE: ignore your neighbors completely. Choose a RANDOM turn in [-60, 60]."),
}


def build_prompt(condition, bearings_deg, dists, relhead_deg, types_rel, max_list=8):
    """Compose the egocentric neighbor list + rule into a user message."""
    order = sorted(range(len(dists)), key=lambda k: dists[k])[:max_list]
    if not order:
        lines = "(no neighbors in range)"
    else:
        rows = []
        for k in order:
            row = f"- bearing {bearings_deg[k]:+d}deg, distance {dists[k]:.1f}, heading {relhead_deg[k]:+d}deg"
            if condition == "segregate":
                row += f", {types_rel[k]}"
            rows.append(row)
        lines = "\n".join(rows)
    return (f"Neighbors you can see:\n{lines}\n\n{RULES[condition]}\n"
            "Output ONLY your turn in degrees as a single integer in [-60, 60].")


_NUM = re.compile(r"-?\d+")


def parse_turn(text: str, max_turn: float = 60.0) -> float:
    """Extract the agent's chosen turn. If the model shows arithmetic, the result
    is after the last '='; otherwise take the last integer (results come last).
    Guards against grabbing an operand from a part-written expression."""
    if not text:
        return 0.0
    tail = text.rsplit("=", 1)[1] if "=" in text else text
    nums = _NUM.findall(tail) or _NUM.findall(text)
    if not nums:
        return 0.0
    try:
        v = float(nums[-1])
    except ValueError:
        return 0.0
    return max(-max_turn, min(max_turn, v))


class Swarm:
    def __init__(self, condition, n, box, radius, speed, model, workers, seed):
        self.condition = condition
        self.n = n
        self.L = box
        self.R = radius
        self.v = speed
        self.model = model
        self.workers = workers
        self.rng = random.Random(seed)
        self.client = anthropic.Anthropic(api_key=os.environ["ANTHROPIC_API_KEY"])
        self.pos = [[self.rng.uniform(0, box), self.rng.uniform(0, box)] for _ in range(n)]
        self.head = [self.rng.uniform(-math.pi, math.pi) for _ in range(n)]
        self.types = [self.rng.randint(0, 1) for _ in range(n)]

    def _neighbors(self, i):
        bearings, dists, relh, trel = [], [], [], []
        px, py = self.pos[i]
        hi = self.head[i]
        cos_i, sin_i = math.cos(-hi), math.sin(-hi)
        for j in range(self.n):
            if j == i:
                continue
            dx = torus_delta(px, self.pos[j][0], self.L)
            dy = torus_delta(py, self.pos[j][1], self.L)
            dist = math.hypot(dx, dy)
            if dist > self.R or dist < 1e-9:
                continue
            rx = dx * cos_i - dy * sin_i           # rotate into i's egocentric frame
            ry = dx * sin_i + dy * cos_i
            bearings.append(int(round(math.degrees(math.atan2(ry, rx)))))
            dists.append(dist)
            relh.append(int(round(math.degrees(wrap_angle(self.head[j] - hi)))))
            trel.append("SAME" if self.types[j] == self.types[i] else "DIFF")
        return bearings, dists, relh, trel

    def _same_frac(self, i):
        """Fraction of visible neighbors that share agent i's type (Schelling cue)."""
        _, _, _, trel = self._neighbors(i)
        return sum(1 for t in trel if t == "SAME") / len(trel) if trel else 0.0

    def _decide(self, i):
        if self.condition == "random":
            # still an LLM agent, but the rule yields no coordination
            return math.radians(self.rng.uniform(-60, 60))
        b, d, rh, tr = self._neighbors(i)
        if not d:
            return 0.0
        prompt = build_prompt(self.condition, b, d, rh, tr)
        for attempt in range(3):
            try:
                msg = self.client.messages.create(
                    model=self.model, max_tokens=40, system=SYSTEM,
                    messages=[{"role": "user", "content": prompt}])
                return math.radians(parse_turn(msg.content[0].text))
            except Exception:
                if attempt == 2:
                    return 0.0
        return 0.0

    def step(self):
        if self.condition == "random":
            turns = [self._decide(i) for i in range(self.n)]
        else:
            with ThreadPoolExecutor(max_workers=self.workers) as ex:
                turns = list(ex.map(self._decide, range(self.n)))
        # Schelling satisfaction: a segregate agent surrounded by its own kind
        # stays put (near-zero speed), so same-type clusters freeze instead of
        # drifting through each other. The LLM still chooses direction.
        fracs = [self._same_frac(i) for i in range(self.n)] if self.condition == "segregate" else None
        for i in range(self.n):
            self.head[i] = wrap_angle(self.head[i] + turns[i])
            spd = self.v
            if fracs is not None and fracs[i] >= 0.55:
                spd = self.v * 0.1
            vx, vy = spd * math.cos(self.head[i]), spd * math.sin(self.head[i])
            self.pos[i][0] = (self.pos[i][0] + vx) % self.L
            self.pos[i][1] = (self.pos[i][1] + vy) % self.L

    def frame(self):
        f = {"positions": [list(p) for p in self.pos],
             "velocities": [[self.v * math.cos(h), self.v * math.sin(h)] for h in self.head]}
        if self.condition == "segregate":
            f["types"] = list(self.types)
        return f

    def run(self, rounds):
        history = [self.frame()]
        for r in range(rounds):
            self.step()
            history.append(self.frame())
            print(f"  round {r+1}/{rounds} done", flush=True)
        return history


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--condition", required=True, choices=list(RULES))
    ap.add_argument("--n", type=int, default=30)
    ap.add_argument("--rounds", type=int, default=20)
    ap.add_argument("--box", type=float, default=12.0)
    ap.add_argument("--radius", type=float, default=4.0)
    ap.add_argument("--speed", type=float, default=0.7)
    ap.add_argument("--model", default="claude-haiku-4-5-20251001")
    ap.add_argument("--workers", type=int, default=30)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    print(f"T2d swarm: condition={a.condition} n={a.n} rounds={a.rounds} model={a.model}", flush=True)
    sw = Swarm(a.condition, a.n, a.box, a.radius, a.speed, a.model, a.workers, a.seed)
    history = sw.run(a.rounds)
    out = {"condition": a.condition,
           "params": {"n": a.n, "rounds": a.rounds, "box": a.box, "radius": a.radius,
                      "speed": a.speed, "model": a.model, "seed": a.seed},
           "history": history}
    json.dump(out, open(a.out, "w"))
    print(f"wrote {a.out}  ({len(history)} frames)", flush=True)


if __name__ == "__main__":
    main()
