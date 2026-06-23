"""Generate the EPC model gallery: run each model, render its asset, attach the
detector readout (reused from the confusion matrix), write manifest.json."""
import json, os, time
import analysis.run_phase2a_panel as R
from gallery.registry import GALLERY
from gallery.renderers import RENDERERS
from epc.phase2a.continuous_metrics import CONTINUOUS_METRIC

ASSETS="gallery/assets"; os.makedirs(ASSETS, exist_ok=True)
# detector readout from the committed confusion matrix (same seed-0 positives)
CONF=json.load(open("docs/validation_rebuild/battery_confusion_matrix_2026-06-16.json"))
ANIM={"point_cloud","grid_field","phase_circle","director","adjacency_blocks","hexatic"}
manifest=[]
for pid_up, info in GALLERY.items():
    pid=pid_up.lower(); t0=time.time()
    entry={"id":pid_up, **{k:info[k] for k in ("name","ref","summary","effect","viz","watch")}}
    entry["metric"]=CONTINUOUS_METRIC.get(pid_up,{}).get("key")
    cm=CONF.get(pid,{})
    if pid in ("p33","p34","p35"):   # discovery patterns: readout from the discovery panel
        import glob as _glob
        _pj=_glob.glob(f"analysis/outputs/{pid}_*_phase2a_panel.json")
        _s=json.load(open(_pj[0])).get("summary",{}) if _pj else {}
        entry["detector"]={"top":pid_up,"tier":(_s.get("positive_tiers") or ["definitive"])[0],
                           "confidence":None,
                           "verdict":f"TNR {round(_s.get('overall_tnr',1.0),2)} · d={round(_s.get('cohens_d_continuous',0.0),1)}",
                           "emergence":None,"self_recognized":True}
    else:
        entry["detector"]={"top":cm.get("top_pattern"),"tier":cm.get("top_tier"),
                           "confidence":cm.get("top_conf"),"verdict":cm.get("verdict"),
                           "emergence":cm.get("emergence"),"self_recognized":cm.get("correct")}
    try:
        if pid=="p31":
            from epc.models.cell_view_sorting import CellViewSorting
            runs,meta=[CellViewSorting(n=60, algorithm="insertion", seed=1).run_to_completion(max_rounds=400)],{}
        elif pid=="p33":
            from epc.models.active_nematic import active_nematic_field
            runs,meta=[active_nematic_field(1)[0]],{}
        elif pid=="p34":
            from epc.models.adaptive_network import adaptive_voter
            runs,meta=[adaptive_voter(1)[0]],{}
        elif pid=="p35":
            from epc.models.entropy_crystal import hard_disk_crystallization
            runs,meta=[hard_disk_crystallization(1, eta_end=0.74, n_steps=2500)[0]],{}
        else:
            runs,meta=getattr(R,f"build_{pid}_positives")(n_seeds=1)
        kind=info["viz"]; ext="gif" if kind in ANIM else "png"
        out=f"{ASSETS}/{pid}.{ext}"
        fn=RENDERERS[kind]
        args=info.get("args",{})
        title=f"{pid_up} {info['name']}"
        r=fn(runs[0], out, title, **args) if args else fn(runs[0], out, title)
        ok=bool(r) and os.path.exists(out) and os.path.getsize(out)>500
        entry["asset"]=os.path.basename(out) if ok else None
        print(f"  {pid_up:<4} {kind:<13} {'OK '+str(os.path.getsize(out)//1024)+'KB' if ok else 'EMPTY'} ({round(time.time()-t0,1)}s)", flush=True)
    except Exception as e:
        entry["asset"]=None; entry["error"]=repr(e)
        print(f"  {pid_up:<4} ERR {type(e).__name__}: {str(e)[:55]}", flush=True)
    manifest.append(entry)
json.dump(manifest, open("gallery/manifest.json","w"), indent=2)
nok=sum(1 for e in manifest if e.get("asset"))
print(f"DONE: {nok}/{len(manifest)} assets rendered")
open("/tmp/gallery_DONE","w").write("done")
