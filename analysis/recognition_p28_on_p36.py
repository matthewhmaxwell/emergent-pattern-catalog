"""Wealth-condensation detector (P28) on the canonical P36 positives: tier and gate metrics per seed.
Follows up the recognition re-run v2, where P28 fired on all five P36 runs.
  python -m analysis.recognition_p28_on_p36
"""
import json
import numpy as np
import analysis.run_phase2a_panel as R
from analysis.battery_profile import build_detector_fns

def short(x):
    if isinstance(x, (float, np.floating)): return round(float(x), 4)
    if isinstance(x, (bool, np.bool_)): return bool(x)
    if isinstance(x, (int, np.integer)): return int(x)
    if isinstance(x, str) or x is None: return x
    if hasattr(x, "value") and hasattr(x, "name"): return str(x.name)
    if isinstance(x, dict): return {str(k): short(v) for k, v in x.items() if not isinstance(v, (list, np.ndarray))}
    return None

def main():
    fns = build_detector_fns(); runs, meta = R.build_p36_positives(n_seeds=5)
    md = meta if isinstance(meta, dict) else None; out = []
    for s, h in enumerate(runs):
        r28 = fns["P28"](h, md); r36 = fns["P36"](h, md)
        keep = lambda r: {k: short(v) for k, v in (r.items() if isinstance(r, dict) else vars(r).items()) if short(v) is not None}
        out.append({"seed": s, "P28": keep(r28), "P36": keep(r36)})
        k28 = keep(r28); print(s, k28.get("tier"), k28.get("primary_metric") or k28.get("primary"), k28.get("secondary_metrics") or k28.get("secondary"), flush=True)
    json.dump(out, open("analysis/outputs/recognition_p28_on_p36.json", "w"), indent=1)

if __name__ == "__main__":
    main()
