"""Enrich the gallery manifest with the REAL validation evidence the catalogue
already has but the gallery wasn't surfacing: Phase-2a negative-control panel
counts (synthetic nulls + catalogue look-alikes + failed regimes) and cross-model
generalization. Pure data join over existing analysis/outputs + docs; no re-render.

Used both standalone (fast path: updates gallery/manifest.json in place) and from
generate_all.py (durable: called on the freshly built manifest before dump)."""
import json, os

PANEL_DIR = "analysis/outputs"
CROSS = "docs/validation_rebuild/cross_model_generalization_2026-06-16.json"


def _panel_counts(pid_low):
    f = os.path.join(PANEL_DIR, f"{pid_low}_phase2a_panel.json")
    if not os.path.exists(f):
        return None
    d = json.load(open(f))
    breakdown, rejected, total = {}, 0, 0
    for k, v in d.items():
        if not isinstance(v, list):
            continue
        cnt = 0
        for e in v:
            if isinstance(e, dict) and "rejection_stage" in e:
                cnt += 1
                total += 1
                if e.get("rejection_stage") != "FIRED" and not e.get("detected", False):
                    rejected += 1
        if cnt:
            breakdown[k] = cnt
    if not total:
        return None
    return {"neg_total": total, "neg_rejected": rejected, "neg_breakdown": breakdown}


def enrich(manifest, cross_path=CROSS):
    cross = json.load(open(cross_path)) if os.path.exists(cross_path) else {}
    for e in manifest:
        low = e["id"].lower()
        det = e.setdefault("detector", {})
        pc = _panel_counts(low)
        if pc:
            det.update(pc)
        if low in cross:
            c = cross[low]
            det["cross_model"] = {"alt": c.get("alt"), "tier": c.get("tier"),
                                  "verdict": c.get("verdict"), "generalizes": c.get("generalizes")}
    return manifest


if __name__ == "__main__":
    mp = "gallery/manifest.json"
    man = json.load(open(mp))
    enrich(man)
    json.dump(man, open(mp, "w"), indent=2)
    np = sum(1 for e in man if e["detector"].get("neg_total"))
    nc = sum(1 for e in man if e["detector"].get("cross_model"))
    print(f"enriched manifest: {np} patterns with panel counts, {nc} with cross-model")
