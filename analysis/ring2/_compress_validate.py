"""Does compressibility earn its place? Long binary sequences: random bits (ratio~1 = null),
periodic (ratio<<1), Thue-Morse (deterministic aperiodic), rule-90 Sierpinski row (deterministic
fractal). Want: structured sequences ratio clearly <1; random ~1. (These are the cases the local
permutation-complexity lens reads as noise.)"""
import numpy as np
from epc.metrics.compressibility import compressibility

N = 2048


def F(seq):
    return [{"symbols": np.asarray(seq, dtype=np.int8)}]


def random_bits(seed=0):
    return np.random.default_rng(seed).integers(0, 2, N)


def periodic(seed=0):
    return np.tile([0, 1, 1, 0, 1], N // 5 + 1)[:N]


def thue_morse(seed=0):
    s = [0]
    while len(s) < N:
        s = s + [1 - x for x in s]
    return np.array(s[:N])


def rule90(seed=0):
    # elementary CA rule 90 (XOR) -> Sierpinski; flatten the space-time raster
    W = 64; T = N // W
    row = np.zeros(W, dtype=int); row[W // 2] = 1
    out = [row.copy()]
    for _ in range(T - 1):
        row = np.roll(row, 1) ^ np.roll(row, -1)
        out.append(row.copy())
    return np.array(out).ravel()[:N]


print(f"{'system':<14}{'lz_ratio':>10}{'lz_norm':>9}")
rows = {}
for nm, fn in [("random_bits", random_bits), ("periodic", periodic),
               ("thue_morse", thue_morse), ("rule90", rule90)]:
    r = compressibility(F(fn(1)))
    rows[nm] = r
    print(f"{nm:<14}{r['lz_ratio']:>10.3f}{r['lz_norm']:>9.3f}")

rnd = rows["random_bits"]["lz_ratio"]
struct = max(rows["periodic"]["lz_ratio"], rows["thue_morse"]["lz_ratio"], rows["rule90"]["lz_ratio"])
print(f"\nstructured(max) lz_ratio {struct:.3f} vs random {rnd:.3f}")
print(f"VERDICT: {'ADMIT' if struct < 0.8 * rnd and rnd > 0.85 else 'review'}")
