"""Assemble the static EPC model gallery (index.html) from manifest.json + assets.
Self-contained: manifest embedded inline, assets referenced from ./assets/."""
import json, html, os

M=json.load(open("gallery/manifest.json"))
M.sort(key=lambda e:int(e["id"][1:]))

def esc(x): return html.escape(str(x if x is not None else ""))

VERDICT_COLOR={"MATCH":"#2f855a","EMERGENT-UNCLASSIFIED":"#b7791f","NO-EMERGENCE":"#718096"}

page=f"""<!doctype html><html lang=en><head><meta charset=utf-8>
<meta name=viewport content="width=device-width,initial-scale=1">
<title>EPC — Model Gallery</title>
<style>
:root{{--bg:#0f1419;--panel:#1a212b;--ink:#e6edf3;--mut:#8b98a5;--acc:#4a9eff;--line:#2b3440}}
*{{box-sizing:border-box}}body{{margin:0;font:15px/1.5 -apple-system,Segoe UI,Roboto,sans-serif;background:var(--bg);color:var(--ink)}}
header{{padding:14px 20px;border-bottom:1px solid var(--line);background:var(--panel)}}
header h1{{margin:0;font-size:18px}}header p{{margin:3px 0 0;color:var(--mut);font-size:13px}}
.wrap{{display:flex;height:calc(100vh - 64px)}}
.list{{width:300px;overflow-y:auto;border-right:1px solid var(--line);background:var(--panel)}}
.item{{padding:9px 16px;cursor:pointer;border-bottom:1px solid var(--line);font-size:14px}}
.item:hover{{background:#222c38}}.item.sel{{background:#243244;border-left:3px solid var(--acc)}}
.item b{{color:var(--acc);margin-right:6px}}.item small{{color:var(--mut);display:block;font-size:11px;margin-top:1px}}
.detail{{flex:1;overflow-y:auto;padding:24px 28px}}
.cols{{display:flex;gap:26px;flex-wrap:wrap}}
.viz{{flex:0 0 auto}}.viz img{{max-width:380px;border:1px solid var(--line);border-radius:8px;background:#fff}}
.info{{flex:1;min-width:280px}}
h2{{margin:0 0 2px;font-size:22px}}.ref{{color:var(--mut);font-style:italic;margin-bottom:16px}}
.k{{color:var(--mut);font-size:12px;text-transform:uppercase;letter-spacing:.04em;margin-top:16px}}
.v{{margin:3px 0 0}}
.card{{margin-top:18px;padding:14px 16px;background:var(--panel);border:1px solid var(--line);border-radius:8px}}
.card h3{{margin:0 0 10px;font-size:13px;color:var(--mut);text-transform:uppercase;letter-spacing:.04em}}
.badge{{display:inline-block;padding:3px 10px;border-radius:12px;color:#fff;font-weight:600;font-size:13px}}
.row{{display:flex;gap:24px;margin-top:10px;flex-wrap:wrap}}.row div b{{display:block;color:var(--mut);font-size:11px;font-weight:500;text-transform:uppercase}}
.row div span{{font-size:16px}}
code{{background:#243244;padding:1px 6px;border-radius:4px;font-size:13px}}
</style></head><body>
<header><h1>Emergent Pattern Catalog — Model Gallery</h1>
<p>32 minimal models of emergent behaviour, each with the calibrated detector readout from the validated battery. Select a pattern.</p></header>
<div class=wrap><div class=list id=list></div><div class=detail id=detail></div></div>
<script>
const M={json.dumps(M)};
const VC={json.dumps(VERDICT_COLOR)};
const list=document.getElementById('list'), detail=document.getElementById('detail');
function show(i){{
  document.querySelectorAll('.item').forEach((e,j)=>e.classList.toggle('sel',j===i));
  const m=M[i], d=m.detector||{{}};
  const vc=VC[d.verdict]||'#718096';
  const asset=m.asset?`<img src="assets/${{m.asset}}" alt="${{m.id}} visualization">`:'<p style=color:#8b98a5>(no visualization)</p>';
  detail.innerHTML=`<div class=cols>
    <div class=viz>${{asset}}</div>
    <div class=info>
      <h2>${{m.id}} · ${{m.name}}</h2><div class=ref>${{m.ref}}</div>
      <div class=k>What it is</div><div class=v>${{m.summary}}</div>
      <div class=k>Emergent effect</div><div class=v>${{m.effect}}</div>
      <div class=k>Canonical metric</div><div class=v><code>${{m.metric||'—'}}</code></div>
      <div class=card><h3>Detector readout (validated battery)</h3>
        ${{d.note&&!d.verdict?`<div style=color:#8b98a5>${{d.note}}</div>`:''}}
        <span class=badge style="background:${{vc}}">${{d.verdict||'—'}}</span>
        ${{d.self_recognized?'<span style="margin-left:10px;color:#2f855a">✓ self-recognized</span>':''}}
        <div class=row>
          <div><b>top pattern</b><span>${{d.top||'—'}}</span></div>
          <div><b>tier</b><span>${{d.tier||'—'}}</span></div>
          <div><b>calibrated confidence</b><span>${{d.confidence!=null?(+d.confidence).toFixed(2):'—'}}</span></div>
          <div><b>generic emergence</b><span>${{d.emergence!=null?(+d.emergence).toFixed(2):'—'}}</span></div>
        </div>
      </div>
    </div></div>`;
}}
M.forEach((m,i)=>{{const e=document.createElement('div');e.className='item';
  e.innerHTML=`<b>${{m.id}}</b>${{m.name}}<small>${{m.ref}}</small>`;e.onclick=()=>show(i);list.appendChild(e);}});
show(0);
</script></body></html>"""
open("gallery/index.html","w").write(page)
print("wrote gallery/index.html", len(page),"bytes;", sum(1 for e in M if e.get("asset")),"assets")
