"""Smoke test: ring2_descriptor with the full 21-lens battery across substrates. Confirms no
import errors, lenses fire on the right substrates and self-guard off them, and tripped/classified
still compute. (Integrity: the tripwire + generic_emergence are unchanged, so tripped/classified
are determined exactly as before -- the new lenses only add fingerprint coordinates.)"""
import numpy as np
from epc.phase2a.ring2_descriptor import ring2_descriptor, SCHEMA

print(f"SCHEMA has {len(SCHEMA)} keys\n")


def positions_hist(n=24, N=300):
    rng = np.random.default_rng(0); pos = rng.uniform(0, 40, (N, 2)); out = []
    for _ in range(n):
        pos = pos + 0.3 * rng.standard_normal((N, 2)); out.append({"positions": pos.copy()})
    return out


def field_hist(n=14, G=64):
    y, x = np.indices((G, G)) * (2 * np.pi / G)
    return [{"field": np.cos(6 * x) + np.cos(6 * y) + 0.1 * np.random.default_rng(t).standard_normal((G, G))} for t in range(n)]


def phases_hist(n=14, G=48):
    base = np.zeros((G, G)); base[:, :G // 2] = 1.0
    return [{"phases": base * 2 + 0.2 * np.random.default_rng(t).standard_normal((G, G))} for t in range(n)]


def theta_field_hist(n=10, G=64):
    return [{"theta_field": (np.random.default_rng(t).uniform(0, np.pi, (G, G)))} for t in range(n)]


def twofield_hist(n=12, G=48):
    rng = np.random.default_rng(0)
    out = []
    for t in range(n):
        C = np.cos(5 * np.indices((G, G))[0] * 2 * np.pi / G)
        A = C + 2 * rng.standard_normal((G, G)); B = C + 2 * rng.standard_normal((G, G))
        out.append({"field_a": A, "field_b": B})
    return out


for name, h in [("positions", positions_hist()), ("field", field_hist()),
                ("phases", phases_hist()), ("theta_field", theta_field_hist()),
                ("two_field", twofield_hist())]:
    d = ring2_descriptor(h)
    print(f"{name:<12} fired({len(d['lenses_fired'])}): {sorted(d['lenses_fired'])}")
    print(f"{'':<12} tripped={d['tripped']} classified={d['classified']}")
print("\nOK: no crashes, descriptor runs across all substrates.")
