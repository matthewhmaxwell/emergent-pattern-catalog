"""Does the directed-info-flow lens earn its place? Coupled-logistic controls with
KNOWN flow structure: a unidirectional CASCADE (x1->x2->...) has high net_asymmetry;
a symmetric diffusive MESH has high mean_te but LOW asymmetry; INDEPENDENT maps have
~0 of both. ADMIT iff mean_te separates coupled-from-independent AND directionality
separates cascade-from-mesh."""
import numpy as np
from epc.metrics.directed_info_flow import directed_transfer_entropy

R, T, N, C = 3.99, 2000, 8, 0.35


def logistic(x):
    return R * x * (1 - x)


def independent(seed=0):
    rng = np.random.default_rng(seed)
    X = rng.random((T, N))
    for t in range(T - 1):
        X[t + 1] = logistic(X[t])
    return X


def cascade(seed=0):                                  # x_{k} drives x_{k+1}
    rng = np.random.default_rng(seed)
    X = rng.random((T, N))
    for t in range(T - 1):
        X[t + 1, 0] = logistic(X[t, 0])
        for k in range(1, N):
            X[t + 1, k] = logistic((1 - C) * X[t, k] + C * X[t, k - 1])
    return X


def mesh(seed=0):                                     # symmetric nearest-neighbour diffusion
    rng = np.random.default_rng(seed)
    X = rng.random((T, N))
    for t in range(T - 1):
        for k in range(N):
            left, right = X[t, (k - 1) % N], X[t, (k + 1) % N]
            X[t + 1, k] = logistic((1 - C) * X[t, k] + C * 0.5 * (left + right))
    return X


rows = []
for s in range(3):
    rows.append((f"independent[{s}]", "independent", directed_transfer_entropy(independent(s))))
    rows.append((f"cascade[{s}]", "cascade", directed_transfer_entropy(cascade(s))))
    rows.append((f"mesh[{s}]", "mesh", directed_transfer_entropy(mesh(s))))

print(f"{'system':<18}{'truth':<14}{'mean_te':>9}{'net_asym':>10}{'direction':>11}")
for nm, t, r in sorted(rows, key=lambda x: -((x[2] or {}).get('directionality', -1))):
    if r:
        print(f"{nm:<18}{t:<14}{r['mean_te']:>9.4f}{r['net_asymmetry']:>10.4f}{r['directionality']:>11.3f}")
    else:
        print(f"{nm:<18}{t:<14}   (n/a)")

g = lambda truth, key: [r[key] for nm, t, r in rows if r and t == truth]
print(f"\nmean_te    coupled(cascade,mesh): {[round(v,3) for v in g('cascade','mean_te')+g('mesh','mean_te')]}  "
      f"independent: {[round(v,3) for v in g('independent','mean_te')]}")
print(f"direction  cascade: {[round(v,2) for v in g('cascade','directionality')]}  "
      f"mesh: {[round(v,2) for v in g('mesh','directionality')]}")
mag_gap = min(g('cascade', 'mean_te') + g('mesh', 'mean_te')) - max(g('independent', 'mean_te'))
dir_gap = min(g('cascade', 'directionality')) - max(g('mesh', 'directionality'))
print(f"\nmean_te gap (coupled - independent) = {mag_gap:+.4f}")
print(f"directionality gap (cascade - mesh) = {dir_gap:+.3f}")
ok = mag_gap > 0.005 and dir_gap > 0.2
print(f"VERDICT: {'ADMIT — magnitude AND direction axes separate' if ok else 'REVIEW'}")
