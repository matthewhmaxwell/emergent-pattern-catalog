"""Assemble the EPC gallery: master list + detail (viz w/ playbar, info,
detector readout, core-algorithm code panel). Self-contained; assets in ./assets/.
Accessible & responsive: keyboard-operable listbox, ARIA on viz/controls,
deep-linkable via URL hash, single mobile breakpoint, reduced-motion aware."""
import json, html
from gallery.education import OVERVIEW_HTML

M = json.load(open("gallery/manifest.json"))
M.sort(key=lambda e: int(e["id"][1:]))
VC = {"MATCH": "#2f855a", "EMERGENT-UNCLASSIFIED": "#b7791f", "NO-EMERGENCE": "#718096"}

CSS = """
:root{--bg:#0f1419;--panel:#1a212b;--ink:#e6edf3;--mut:#8b98a5;--acc:#4a9eff;--line:#2b3440}
*{box-sizing:border-box}body{margin:0;font:15px/1.5 -apple-system,Segoe UI,Roboto,sans-serif;background:var(--bg);color:var(--ink);overflow-x:hidden}
header{padding:14px 20px;border-bottom:1px solid var(--line);background:var(--panel)}
header h1{margin:0;font-size:18px}header p{margin:3px 0 0;color:var(--mut);font-size:13px}
.wrap{display:flex;height:calc(100vh - 64px)}
.list{width:290px;overflow-y:auto;border-right:1px solid var(--line);background:var(--panel);flex:none}
.item{padding:9px 16px;cursor:pointer;border-bottom:1px solid var(--line);font-size:14px}
.item.about-item b{color:#b7791f}
.about{max-width:780px;line-height:1.65}
.about h2{font-size:23px;margin:0 0 14px}
.about p{margin:0 0 13px;color:#cdd9e5}
.about b{color:#fff}
.item:hover{background:#222c38}.item.sel{background:#243244;border-left:3px solid var(--acc)}
.item b{color:var(--acc);margin-right:6px}.item small{color:var(--mut);display:block;font-size:11px;margin-top:1px}
.detail{flex:1;overflow-y:auto;padding:22px 26px}
.cols{display:flex;gap:26px;flex-wrap:wrap}
.viz{flex:0 0 auto;max-width:100%;overflow-x:auto}
.stage{border:1px solid var(--line);border-radius:8px;background-repeat:no-repeat;background-color:#fff}
.viz img{max-width:340px;border:1px solid var(--line);border-radius:8px;background:#fff}
.bar{display:flex;align-items:center;gap:10px;margin-top:10px;width:300px}
.bar button{background:#243244;color:var(--ink);border:1px solid var(--line);border-radius:6px;padding:4px 10px;cursor:pointer;font-size:14px}
.bar button:hover{background:#2d3e54}
.bar select{background:#243244;color:var(--ink);border:1px solid var(--line);border-radius:6px;padding:3px 6px}
a.dl{display:inline-block;margin-top:8px;margin-right:8px;background:#243244;color:var(--acc);border:1px solid var(--line);border-radius:6px;padding:5px 12px;text-decoration:none;font-size:13px}
a.dl:hover{background:#2d3e54}.bar input[type=range]{flex:1}
.watch{width:300px;margin-top:12px;padding:10px 12px;background:#10243a;border:1px solid #21466b;border-left:3px solid var(--acc);border-radius:6px;font-size:13px;color:#cfe3f7}
.watch b{color:var(--acc);display:block;font-size:11px;text-transform:uppercase;letter-spacing:.04em;margin-bottom:3px}
.info{flex:1;min-width:280px}
h2{margin:0 0 2px;font-size:22px}.ref{color:var(--mut);font-style:italic;margin-bottom:14px}
.k{color:var(--mut);font-size:12px;text-transform:uppercase;letter-spacing:.04em;margin:14px 0 0;font-weight:600}.v{margin:3px 0 0}
.card{margin-top:16px;padding:14px 16px;background:var(--panel);border:1px solid var(--line);border-radius:8px}
.card h3{margin:0 0 10px;font-size:13px;color:var(--mut);text-transform:uppercase;letter-spacing:.04em}
.badge{display:inline-block;padding:3px 10px;border-radius:12px;color:#fff;font-weight:600;font-size:13px}
.row{display:flex;gap:24px;margin-top:10px;flex-wrap:wrap}.row div b{display:block;color:var(--mut);font-size:11px;font-weight:500;text-transform:uppercase}.row div span{font-size:16px}
code{background:#243244;padding:1px 6px;border-radius:4px;font-size:13px}
pre.codeblk{margin:8px 0 0;padding:14px 16px;background:#0b0f14;border:1px solid var(--line);border-radius:8px;overflow:auto;max-height:420px;font:12.5px/1.5 ui-monospace,SFMono-Regular,Menlo,monospace;color:#cdd9e5;white-space:pre}
.kw{color:#ff7b72}.cm{color:#6a9955}.st{color:#a5d6ff}.fn{color:#d2a8ff}
details.card>summary.codesum{cursor:pointer;font-size:13px;color:var(--mut);text-transform:uppercase;letter-spacing:.04em}
details.card>summary.codesum code{text-transform:none}
details.card[open]>summary.codesum{margin-bottom:10px}
.skip{position:absolute;left:-9999px;top:0;z-index:20;background:var(--acc);color:#04121f;font-weight:600;padding:8px 14px;border-radius:0 0 6px 0;text-decoration:none}
.skip:focus{left:0}
.item:focus-visible,.bar button:focus-visible,.bar select:focus-visible,.bar input[type=range]:focus-visible,a.dl:focus-visible,summary.codesum:focus-visible{outline:2px solid var(--acc);outline-offset:2px}
.detail:focus{outline:none}h2:focus{outline:none}
@media(max-width:760px){
  .wrap{flex-direction:column;height:auto}
  .list{width:100%;flex:none;max-height:42vh;border-right:none;border-bottom:1px solid var(--line)}
  .detail{padding:16px 14px}
  .cols{gap:16px}
  .viz,.bar,.watch{width:100%;max-width:100%}
  header h1{font-size:16px}header p{font-size:12px}
}
"""

JS = """
const M=__M__, VC=__VC__, OVERVIEW=__OVERVIEW__;
const list=document.getElementById('list'), detail=document.getElementById('detail');
let timer=null, suppressHash=false;
const RM = matchMedia('(prefers-reduced-motion: reduce)').matches;
const KW=/\\b(def|return|for|while|if|elif|else|in|and|or|not|import|from|class|with|as|None|True|False|self|lambda|break|continue|np)\\b/g;
function hl(s){
  s=s.replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;');
  s=s.replace(/(#[^\\n]*)/g,'<span class=cm>$1</span>');
  s=s.replace(KW,'<span class=kw>$1</span>');
  return s;
}
function slug(s){return (s||'').replace(/[^A-Za-z0-9]+/g,'-').replace(/^-+|-+$/g,'');}
function esc(s){return (s||'').replace(/[<>]/g,'').replace(/"/g,'&quot;');}
function clearTimer(){if(timer){clearInterval(timer);timer=null;}}
function setSelected(key){
  document.querySelectorAll('.item').forEach(el=>{
    const on = el.dataset.key===String(key);
    el.classList.toggle('sel', on);
    el.setAttribute('aria-selected', on?'true':'false');
    el.tabIndex = on?0:-1;
  });
}
function focusDetail(){const h=detail.querySelector('h2')||detail; h.setAttribute('tabindex','-1'); h.focus();}

function show(i){
  clearTimer();
  try{
    const m=M[i], d=m.detector||{}, vc=VC[d.verdict]||'#718096';
    setSelected(i);
    const desc=esc(m.watch||m.effect);
    let vizHtml='';
    if(m.asset_type==='sprite' && m.sprite){
      const s=m.sprite;
      const mp4name=m.id+'_'+slug(m.name)+'.mp4', pngname=m.id+'_'+slug(m.name)+'_frames.png';
      vizHtml=`<div class=stage id=stage role=img aria-label="${m.id} ${esc(m.name)} animation. ${desc}" style="width:${s.fw}px;height:${s.fh}px;max-width:100%;background-image:url(assets/${m.asset})"></div>
        <div class=bar>
        <button id=pp aria-label="Pause animation">⏸</button>
        <button id=rs aria-label="Reset to first frame">⟲</button>
        <input type=range id=sl min=0 max=${s.frames-1} value=0 aria-label="Animation frame">
        <span id=fc aria-hidden=true style="color:var(--mut);font-size:12px;min-width:48px">0/${s.frames-1}</span>
        <select id=spd aria-label="Playback speed"><option value=0.25>0.25×</option><option value=0.5>0.5×</option><option value=1 selected>1×</option><option value=2>2×</option><option value=4>4×</option></select></div>
        ${m.mp4?`<a class=dl href="assets/${m.mp4}" download="${mp4name}">⤓ MP4</a>`:''}${m.asset?`<a class=dl href="assets/${m.asset}" download="${pngname}" title="One PNG with all ${s.frames} frames tiled in a grid (a sprite sheet, not a single still)">⤓ Frame sheet (PNG)</a>`:''}`;
    } else if(m.asset){
      vizHtml=`<img src="assets/${m.asset}" alt="${m.id} ${esc(m.name)}: ${desc}">`;
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
        <h2 tabindex=-1>${m.id} · ${m.name}</h2><div class=ref>${m.ref}</div>
        <h3 class=k>What it is</h3><div class=v>${m.summary}</div>
        <h3 class=k>Emergent effect</h3><div class=v>${m.effect}</div>
        <h3 class=k>How it emerges</h3><div class=v>${m.mechanism||''}</div>
        <h3 class=k>Canonical metric</h3><div class=v><code>${m.metric||'—'}</code></div>
        <div class=card><h3>Detector readout (validated battery)</h3>${det}</div>
      </div></div>
      <details class=card open><summary class=codesum>Core algorithm · <code>${m.code_module||''} ${m.code_where||''}</code></summary>
        <pre class=codeblk>${hl(m.code||'# (unavailable)')}</pre></details>`;
    if(m.asset_type==='sprite' && m.sprite) initSprite(m);
  }catch(err){
    detail.innerHTML='<div class=about><h2>Could not render this pattern</h2><p style=color:#8b98a5>'+(err&&err.message||err)+'</p></div>';
  }
}

function initSprite(m){
  const s=m.sprite, stage=document.getElementById('stage'),
    pp=document.getElementById('pp'), rs=document.getElementById('rs'),
    sl=document.getElementById('sl'), fc=document.getElementById('fc'), spd=document.getElementById('spd');
  const BASE=110; let cur=0;
  function setF(k){cur=((k%s.frames)+s.frames)%s.frames;
    const c=cur%s.cols, r=Math.floor(cur/s.cols);
    stage.style.backgroundPosition=`-${c*s.fw}px -${r*s.fh}px`;
    sl.value=cur; fc.textContent=cur+'/'+(s.frames-1);
    sl.setAttribute('aria-valuetext',(cur+1)+' of '+s.frames);}
  function play(){clearTimer(); timer=setInterval(()=>setF(cur+1), BASE/parseFloat(spd.value)); pp.textContent='⏸'; pp.setAttribute('aria-label','Pause animation');}
  function pause(){clearTimer(); pp.textContent='▶'; pp.setAttribute('aria-label','Play animation');}
  pp.onclick=()=>timer?pause():play();
  rs.onclick=()=>{pause();setF(0);};
  sl.oninput=()=>{pause();setF(+sl.value);};
  spd.onchange=()=>{if(timer)play();};
  setF(0);
  if(RM) pause(); else play();
}

function showAbout(){clearTimer(); setSelected('about'); detail.innerHTML='<div class=about>'+OVERVIEW+'</div>';}

function renderKey(k){
  k=(k||'').trim();
  if(k.toLowerCase()==='about'){showAbout();return true;}
  const i=M.findIndex(m=>m.id.toLowerCase()===k.toLowerCase());
  if(i>=0){show(i);return true;}
  return false;
}
function navigate(k){
  const target='#'+k;
  if(location.hash!==target){suppressHash=true; location.hash=k;}
  renderKey(k); focusDetail();
}
window.addEventListener('hashchange',()=>{ if(suppressHash){suppressHash=false;return;} renderKey((location.hash||'').replace(/^#/,'')); });

function makeItem(key, nav, htmlStr){
  const e=document.createElement('div');
  e.className='item'+(key==='about'?' about-item':'');
  e.dataset.key=String(key); e.dataset.nav=nav;
  if(key==='about') e.id='about-item';
  e.setAttribute('role','option'); e.setAttribute('aria-selected','false'); e.tabIndex=-1;
  e.innerHTML=htmlStr;
  e.onclick=()=>navigate(nav);
  list.appendChild(e);
}
makeItem('about','about','<b>✦</b>About emergence<small>start here — what all 32 share</small>');
M.forEach((m,i)=>makeItem(i, m.id, `<b>${m.id}</b>${m.name}<small>${m.ref}</small>`));

list.addEventListener('keydown', e=>{
  const its=[...list.querySelectorAll('.item')], idx=its.indexOf(document.activeElement);
  if(idx<0) return;
  if(e.key==='Enter'||e.key===' '){e.preventDefault(); navigate(its[idx].dataset.nav); return;}
  let n=-1;
  if(e.key==='ArrowDown') n=Math.min(its.length-1,idx+1);
  else if(e.key==='ArrowUp') n=Math.max(0,idx-1);
  else if(e.key==='Home') n=0;
  else if(e.key==='End') n=its.length-1;
  else return;
  e.preventDefault();
  its.forEach(el=>el.tabIndex=-1); its[n].tabIndex=0; its[n].focus();
});

if(!renderKey((location.hash||'').replace(/^#/,''))) showAbout();
"""

page = ("<!doctype html><html lang=en><head><meta charset=utf-8>"
    "<meta name=viewport content='width=device-width,initial-scale=1'>"
    "<title>EPC — Model Gallery</title><style>" + CSS + "</style></head><body>"
    "<a class=skip href='#detail'>Skip to content</a>"
    "<header><h1>Emergent Pattern Catalog — Model Gallery</h1>"
    "<p>32 minimal models of emergent behaviour — start with <b>About emergence</b>, then explore each: play the effect, read how &amp; why it emerges, see the simple rule, and watch the validated detector recognise it.</p></header>"
    "<div class=wrap><div class=list id=list role=listbox aria-label='Emergent pattern models — use arrow keys to browse, Enter to open'></div>"
    "<main class=detail id=detail tabindex=-1></main></div>"
    "<noscript><p style='padding:20px;color:#8b98a5'>This interactive gallery requires JavaScript to be enabled.</p></noscript>"
    "<script>" + JS.replace("__M__", json.dumps(M)).replace("__VC__", json.dumps(VC)).replace("__OVERVIEW__", json.dumps(OVERVIEW_HTML)) + "</script></body></html>")
open("gallery/index.html", "w").write(page)
print("wrote gallery/index.html", len(page), "bytes;",
      sum(1 for e in M if e.get("asset_type") == "sprite"), "playable,",
      sum(1 for e in M if e.get("code")), "with code")
