"""Self-recognition and cross-model recognition, re-run with five seeds each after the null-draw default fix
(P13, P22, P27 factories now 199 draws, as in their panels). Writes analysis/outputs/recognition_rerun.json.

  python -m analysis.recognition_rerun
"""
import json, time
import analysis.run_phase2a_panel as R
import analysis.t2c_systems as T
import analysis.transfer_matrix_audit as M
from analysis.battery_profile import profile_observation, build_detector_fns

PATTERNS = [f"p{i}" for i in range(1, 33) if i != 31]
CROSS = {"P16": ("recog_p16_boolean_grn", "BooleanGRN"), "P20": ("recog_p20_fraction_threshold", "FractionThreshold"),
         "P23": ("recog_p23_el_farol", "El Farol"), "P24": ("recog_p24_integral", "integral controller"),
         "P25": ("recog_p25_multibasin", "MultiBasinGRN"), "P26": ("recog_p26_threshold_unit", "ThresholdUnit")}


def row(res):
    v, top = res["verdict"], res.get("top") or {}
    fired = [r["pattern_id"] for r in res["profile"] if r.get("fired_detected")]
    return {"verdict": v.get("verdict"), "match": v.get("pattern_id"), "top": top.get("pattern_id"),
            "top_tier": top.get("detector_tier"), "fired": fired}


def main():
    t0 = time.time(); fns = build_detector_fns(); out = {"self": {}, "cross": {}}
    for p in PATTERNS:
        runs, meta = getattr(R, f"build_{p}_positives")(n_seeds=5)
        md = meta if isinstance(meta, dict) else None
        out["self"][p.upper()] = [row(profile_observation(h, md, detector_fns=fns, seed=s)) for s, h in enumerate(runs)]
        print(p.upper(), [(r["top"], r["top_tier"]) for r in out["self"][p.upper()]], f"{time.time()-t0:.0f}s", flush=True)
    for p, (fn, name) in CROSS.items():
        rs = []
        for s in range(5):
            h, meta = getattr(T, fn)(s); rs.append(row(profile_observation(h, meta if isinstance(meta, dict) else None, detector_fns=fns, seed=s)))
        out["cross"][p] = {"model": name, "runs": rs}
        print("cross", p, [(r["top"], r["top_tier"]) for r in rs], flush=True)
    rs = [row(profile_observation(h, md, detector_fns=fns, seed=s)) for s, (h, md) in enumerate(M._physarum())]
    out["cross"]["P29"] = {"model": "Physarum", "runs": rs}
    print("cross P29", [(r["top"], r["top_tier"]) for r in rs], flush=True)
    out["elapsed_s"] = round(time.time() - t0)
    json.dump(out, open("analysis/outputs/recognition_rerun.json", "w"), indent=1, default=str)
    print("RECOGNITION DONE", out["elapsed_s"], flush=True)


if __name__ == "__main__":
    main()
