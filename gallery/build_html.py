"""Assemble the EPC gallery: master list + detail (viz w/ playbar, info,
detector readout, core-algorithm code panel). Self-contained; assets in ./assets/."""
import json, html

M = json.load(open("gallery/manifest.json"))
M.sort(key=lambda e: int(e["id"][1:]))
VC = {"MATCH": "#2f855a", "EMERGENT-UNCLASSIFIED": "#b7791f", "NO-EMERGENCE": "#718096"}

CSS = """
:root{--bg:#0f1419;--panel:#1a212b;--ink:#e6edf3;--mut:#8b98a5;--acc:#4a9eff;--line:#2b3440}
*{box-sizing:border-box}body{margin:0;font:15px/1.5 -apple-system,Segoe UI,Roboto,sans-serif;background:var(--bg);color:var(--ink)}
header{padding:14px 20px;border-bottom:1px solid var(--line);background:var(--panel)}
header h1{margin:0;font-size:18px}header p{margin:3px 0 0;color:var(--mut);font-size:13px}
.wrap{display:flex;height:calc(100vh - 64px)}
.list{width:290px;overflow-y:auto;border-right:1px solid var(--line);background:var(--panel);flex:none}
.item{padding:9px 16px;cursor:pointer;border-bottom:1px solid var(--line);font-size:14px}
.item:hover{background:#222c38}.item.sel{background:#243244;border-left:3px solid var(--acc)}
.item b{color:var(--acc);margin-right:6px}.item small{color:var(--mut);display:block;font-size:11px;margin-top:1px}
.detail{flex:1;overflow-y:auto;padding:22px 26px}
.cols{display:flex;gap:26px;flex-wrap:wrap}
.viz{flex:0 0 auto}
.stage{border:1px solid var(--line);border-radius:8px;background-repeat:no-repeat;background-color:#fff}
.viz img{max-width:340px;border:1px solid var(--line);border-radius:8px;background:#fff}
.bar{display:flex;align-items:center;gap:10px;margin-top:10px;width:300px}
.bar button{background:#243244;color:var(--ink);border:1px solid var(--line);border-radius:6px;padding:4px 10px;cursor:pointer;font-size:14px}
.bar button:hover{background:#2d3e54}
.bar select{background:#243244;color:var(--ink);border:1px solid var(--line);border-radius:6px;padding:3px 6px}
a.dl{display:inline-block;margin-top:8px;background:#243244;color:var(--acc);border:1px solid var(--line);border-radius:6px;padding:5px 12px;text-decoration:none;font-size:13px}
a.dl:hover{background:#2d3e54}.bar input[type=range]{flex:1}
.watch{width:300px;margin-top:12px;padding:10px 12px;background:#10243a;border:1px solid #21466b;border-left:3px solid var(--acc);border-radius:6px;font-size:13px;color:#cfe3f7}
.watch b{color:var(--acc);display:block;font-size:11px;text-transform:uppercase;letter-spacing:.04em;margin-bottom:3px}
.info{flex:1;min-width:280px}
h2{margin:0 0 2px;font-size:22px}.ref{color:var(--mut);font-style:italic;margin-bottom:14px}
.k{color:var(--mut);font-size:12px;text-transform:uppercase;letter-spacing:.04em;margin-top:14px}.v{margin:3px 0 0}
.card{margin-top:16px;padding:14px 16px;background:var(--panel);border:1px solid var(--line);border-radius:8px}
.card h3{margin:0 0 10px;font-size:13px;color:var(--mut);text-transform:uppercase;letter-spacing:.04em}
.badge{display:inline-block;padding:3px 10px;border-radius:12px;color:#fff;font-weight:600;font-size:13px}
.row{display:flex;gap:24px;margin-top:10px;flex-wrap:wrap}.row div b{display:block;color:var(--mut);font-size:11px;font-weight:500;text-transform:uppercase}.row div span{font-size:16px}
code{background:#243244;padding:1px 6px;border-radius:4px;font-size:13px}
pre.codeblk{margin:8px 0 0;padding:14px 16px;background:#0b0f14;border:1px solid var(--line);border-radius:8px;overflow:auto;max-height:420px;font:12.5px/1.5 ui-monospace,SFMono-Regular,Menlo,monospace;color:#cdd9e5;white-space:pre}
.kw{color:#ff7b72}.cm{color:#6a9955}.st{color:#a5d6ff}.fn{color:#d2a8ff}
"""

JS = """
const M=__M__, VC=__VC__;
const list=document.getElementById('list'), detail=document.getElementById('detail');
let timer=null;
const KW=/\\b(def|return|for|while|if|elif|else|in|and|or|not|import|from|class|with|as|None|True|False|self|lambda|break|continue|np)\\b/g;
function hl(s){
  s=s.replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;');
  s=s.replace(/(#[^\\n]*)/g,'<span class=cm>$1</span>');
  s=s.replace(KW,'<span class=kw>$1</span>');
  return s;
}
function show(i){
  if(timer){clearInterval(timer);timer=null;}
  document.querySelectorAll('.item').forEach((e,j)=>e.classList.toggle('sel',j===i));
  const m=M[i], d=m.detector||{}, vc=VC[d.verdict]||'#718096';
  let vizHtml='';
  if(m.asset_type==='sprite' && m.sprite){
    const s=m.sprite;
    vizHtml=`<div class=stage id=stage style="width:${s.fw}px;height:${s.fh}px;background-image:url(assets/${m.asset})"></div>
      <div class=bar><button id=pp>⏸</button><button id=rs>⟲</button>
      <input type=range id=sl min=0 max=${s.frames-1} value=0><span id=fc style="color:var(--mut);font-size:12px;min-width:48px">0/${s.frames-1}</span>
      <select id=spd title="playback speed"><option value=0.25>0.25×</option><option value=0.5>0.5×</option><option value=1 selected>1×</option><option value=2>2×</option><option value=4>4×</option></select></div>
      ${m.gif?`<a class=dl href="assets/${m.gif}" download>\u2913 Download GIF</a>`:''}`;
  } else if(m.asset){
    vizHtml=`<img src="assets/${m.asset}" alt="${m.id}">`;
  } else { vizHtml='<p style=color:#8b98a5>(no visualization)</p>'; }
  const watchHtml=m.watch?`<div class=watch><b>&#128065; What to look for</b>${m.watch}</div>`:'';
  const det = d.verdict
    ? `<span class=badge style="background:${vc}">${d.verdict}</span>${d.self_recognized?'<span style="margin-left:10px;color:#2f855a">✓ self-recognized</span>':''}
       <div class=row><div><b>top</b><span>${d.top||'—'}</span></div><div><b>tier</b><span>${d.tier||'—'}</span></div>
       <div><b>calibrated confidence</b><span>${d.confidence!=null?(+d.confidence).toFixed(2):'—'}</span></div>
       <div><b>generic emergence</b><span>${d.emergence!=null?(+d.emergence).toFixed(2):'—'}</span></div></div>`
    : `<div style=color:#8b98a5>${(d.note)||'—'}</div>`;
  detail.innerHTML=`<div class=cols>
    <div class=viz>${vizHtml}${watchHtml}</div>
    <div class=info>
      <h2>${m.id} · ${m.name}</h2><div class=ref>${m.ref}</div>
      <div class=k>What it is</div><div class=v>${m.summary}</div>
      <div class=k>Emergent effect</div><div class=v>${m.effect}</div>
      <div class=k>Canonical metric</div><div class=v><code>${m.metric||'—'}</code></div>
      <div class=card><h3>Detector readout (validated battery)</h3>${det}</div>
    </div></div>
    <div class=card><h3>Core algorithm &nbsp;·&nbsp; <code>${m.code_module||''} ${m.code_where||''}</code></h3>
      <pre class=codeblk>${hl(m.code||'# (unavailable)')}</pre></div>`;
  if(m.asset_type==='sprite' && m.sprite){
    const s=m.sprite, stage=document.getElementById('stage'),
      pp=document.getElementById('pp'), rs=document.getElementById('rs'),
      sl=document.getElementById('sl'), fc=document.getElementById('fc');
    const spd=document.getElementById('spd'); const BASE=110; let cur=0;
    function setF(k){cur=((k%s.frames)+s.frames)%s.frames;
      const c=cur%s.cols, r=Math.floor(cur/s.cols);
      stage.style.backgroundPosition=`-${c*s.fw}px -${r*s.fh}px`;
      sl.value=cur; fc.textContent=cur+'/'+(s.frames-1);}
    function play(){if(timer)clearInterval(timer); timer=setInterval(()=>setF(cur+1), BASE/parseFloat(spd.value)); pp.textContent='⏸';}
    function pause(){if(timer){clearInterval(timer);timer=null;} pp.textContent='▶';}
    pp.onclick=()=>timer?pause():play();
    rs.onclick=()=>{pause();setF(0);};
    sl.oninput=()=>{pause();setF(+sl.value);};
    spd.onchange=()=>{if(timer)play();};
    setF(0); play();
  }
}
M.forEach((m,i)=>{const e=document.createElement('div');e.className='item';
  e.innerHTML=`<b>${m.id}</b>${m.name}<small>${m.ref}</small>`;e.onclick=()=>show(i);list.appendChild(e);});
show(0);
"""

page = ("<!doctype html><html lang=en><head><meta charset=utf-8>"
    "<meta name=viewport content='width=device-width,initial-scale=1'>"
    "<title>EPC — Model Gallery</title><style>" + CSS + "</style></head><body>"
    "<header><h1>Emergent Pattern Catalog — Model Gallery</h1>"
    "<p>32 minimal models of emergent behaviour — play the effect, read the simple rule, see the validated detector recognise it.</p></header>"
    "<div class=wrap><div class=list id=list></div><div class=detail id=detail></div></div>"
    "<script>" + JS.replace("__M__", json.dumps(M)).replace("__VC__", json.dumps(VC)) + "</script></body></html>")
open("gallery/index.html", "w").write(page)
print("wrote gallery/index.html", len(page), "bytes;",
      sum(1 for e in M if e.get("asset_type") == "sprite"), "playable,",
      sum(1 for e in M if e.get("code")), "with code")
