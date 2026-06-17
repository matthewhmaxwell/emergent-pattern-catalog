"""Assemble the EPC gallery: master list + detail (viz w/ playbar, info,
detector readout, core-algorithm code panel). Self-contained; assets in ./assets/.
Accessible & responsive: keyboard-operable listbox, ARIA on viz/controls,
deep-linkable via URL hash, single mobile breakpoint, reduced-motion aware."""
import json, html
from gallery.education import OVERVIEW_HTML, METHODS_HTML

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
.dl{display:inline-block;margin-top:8px;margin-right:8px;background:#243244;color:var(--acc);border:1px solid var(--line);border-radius:6px;padding:5px 12px;text-decoration:none;font:inherit;font-size:13px;cursor:pointer}
.dl:hover{background:#2d3e54}.bar input[type=range]{flex:1}
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
.about h3{font-size:16px;margin:20px 0 7px;color:#fff;text-transform:none;letter-spacing:0}
.about ul{margin:0 0 13px;padding-left:20px;color:#cdd9e5}.about li{margin:0 0 7px}
.about a{color:var(--acc);text-decoration:underline}
.item.methods-item b{color:#3fae6b}
.tier{display:inline-block;margin-left:8px;padding:2px 8px;border:1px solid var(--line);border-radius:10px;font-size:11px;color:var(--mut);text-transform:uppercase;letter-spacing:.04em;vertical-align:middle;cursor:help}
.disc{margin-top:11px;font-size:13px}.disc .lab{display:block;color:var(--mut);font-size:11px;text-transform:uppercase;letter-spacing:.04em;margin-bottom:2px;cursor:help}
.disc strong{color:#7ee2a8}
.cmline{margin-top:8px;font-size:13px;color:#cfe3f7;cursor:help}
.foot{margin-top:9px;font-size:12px;color:var(--mut)}.foot .hint{cursor:help;border-bottom:1px dotted var(--line)}
.mut2{color:var(--mut)}
.methlink{color:var(--acc);text-decoration:none}.methlink:hover{text-decoration:underline}
.ovl{position:fixed;inset:0;z-index:50;display:none;flex-direction:column;align-items:center;justify-content:center;gap:14px;background:rgba(6,9,13,.93);padding:18px}
.ovl.open{display:flex}
.ovl-top{width:100%;max-width:1100px;display:flex;align-items:center;justify-content:space-between;gap:12px}
.ovl-cap{font-size:15px;color:var(--ink);font-weight:600}
.ovl-actions{display:flex;gap:8px;flex:none}
.ovl-btn{background:#243244;color:var(--ink);border:1px solid var(--line);border-radius:6px;padding:6px 12px;cursor:pointer;font:inherit;font-size:13px}
.ovl-btn:hover{background:#2d3e54}
.ovl-stage{background-repeat:no-repeat;background-color:#fff;border:1px solid var(--line);border-radius:8px;max-width:100%}
.ovl-bar{display:flex;align-items:center;gap:10px;width:100%;max-width:560px}
.ovl-cbtn{background:#243244;color:var(--ink);border:1px solid var(--line);border-radius:6px;padding:4px 10px;cursor:pointer;font:inherit;font-size:14px}
.ovl-cbtn:hover{background:#2d3e54}
.ovl-range{flex:1}
.ovl-fc{color:var(--mut);font-size:12px;min-width:48px}
.ovl-sel{background:#243244;color:var(--ink);border:1px solid var(--line);border-radius:6px;padding:3px 6px;font:inherit}
.ovl-btn:focus-visible,.ovl-cbtn:focus-visible,.ovl-range:focus-visible,.ovl-sel:focus-visible{outline:2px solid var(--acc);outline-offset:2px}
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
const M=__M__, VC=__VC__, OVERVIEW=__OVERVIEW__, METHODS=__METHODS__;
const list=document.getElementById('list'), detail=document.getElementById('detail');
let timer=null, suppressHash=false, paintBigRef=null, lastEnl=null;
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

// ---- in-page enlarge / fullscreen overlay (no popup window — cannot be blocked) ----
const ovl=document.createElement('div');
ovl.id='ovl'; ovl.className='ovl';
ovl.setAttribute('role','dialog'); ovl.setAttribute('aria-modal','true');
ovl.setAttribute('aria-label','Enlarged visualization'); ovl.setAttribute('aria-hidden','true');
ovl.innerHTML=`<div class=ovl-top><span id=ovlcap class=ovl-cap></span><span class=ovl-actions>
  <button id=ovlfs class=ovl-btn type=button>⛶ Fullscreen</button>
  <button id=ovlx class=ovl-btn type=button aria-label="Close enlarged view">✕ Close</button></span></div>
  <div id=ovlstage class=ovl-stage role=img></div>
  <div class=ovl-bar><button id=ovlpp class=ovl-cbtn type=button aria-label="Pause animation">⏸</button>
  <button id=ovlrs class=ovl-cbtn type=button aria-label="Reset to first frame">⟲</button>
  <input id=ovlsl class=ovl-range type=range min=0 value=0 aria-label="Animation frame">
  <span id=ovlfc class=ovl-fc aria-hidden=true></span>
  <select id=ovlspd class=ovl-sel aria-label="Playback speed"><option value=0.25>0.25×</option><option value=0.5>0.5×</option><option value=1 selected>1×</option><option value=2>2×</option><option value=4>4×</option></select></div>`;
document.body.appendChild(ovl);
const OVL={el:ovl, cap:ovl.querySelector('#ovlcap'), stage:ovl.querySelector('#ovlstage'),
  pp:ovl.querySelector('#ovlpp'), rs:ovl.querySelector('#ovlrs'), sl:ovl.querySelector('#ovlsl'),
  fc:ovl.querySelector('#ovlfc'), spd:ovl.querySelector('#ovlspd'), fs:ovl.querySelector('#ovlfs'), x:ovl.querySelector('#ovlx')};
function fsActive(){return document.fullscreenElement||document.webkitFullscreenElement;}
function closeOverlay(){
  if(fsActive() && document.exitFullscreen) document.exitFullscreen().catch(()=>{});
  ovl.classList.remove('open'); ovl.setAttribute('aria-hidden','true'); paintBigRef=null;
  if(lastEnl && document.contains(lastEnl)) lastEnl.focus();
  lastEnl=null;
}
OVL.x.onclick=closeOverlay;
ovl.addEventListener('click', e=>{ if(e.target===ovl) closeOverlay(); });
document.addEventListener('keydown', e=>{ if(e.key==='Escape' && ovl.classList.contains('open')) closeOverlay(); });
const reqFS=ovl.requestFullscreen||ovl.webkitRequestFullscreen;
if(!reqFS) OVL.fs.style.display='none';
OVL.fs.onclick=()=>{ if(!fsActive()){ reqFS&&reqFS.call(ovl); } else if(document.exitFullscreen){ document.exitFullscreen(); } };
document.addEventListener('fullscreenchange', ()=>{ OVL.fs.textContent=fsActive()?'⛶ Exit fullscreen':'⛶ Fullscreen'; if(ovl.classList.contains('open')&&paintBigRef) paintBigRef(true); });
window.addEventListener('resize', ()=>{ if(ovl.classList.contains('open')&&paintBigRef) paintBigRef(true); });

function methLink(){return `<div class=foot><a href="#methods" class=methlink onclick="navigate('methods');return false;">How we validate ↗</a></div>`;}
function detectorHtml(d){
  const vc=VC[d.verdict]||'#718096';
  const TG={definitive:'strongest tier — the detector fired its most specific gate',
            confirmation:'mid tier — a corroborating gate fired',
            screening:'weakest tier — only a coarse screening gate fired',
            none:'did not reach a firing tier on this run'};
  if(!d.verdict) return `<div style=color:#8b98a5>${d.note||'—'}</div>`+methLink();
  let disc='';
  if(d.neg_total){
    const b=d.neg_breakdown||{}, parts=[];
    if(b.synthetic) parts.push(b.synthetic+' synthetic nulls');
    if(b.catalog) parts.push(b.catalog+' look-alikes');
    if(b.failed_regime) parts.push(b.failed_regime+' failed regimes');
    disc=`<div class=disc><span class=lab title="Negative controls the detector had to reject: structureless nulls, the other catalogue patterns, and the model run where the pattern never forms.">Discrimination · negative controls</span>rejected <strong>${d.neg_rejected}/${d.neg_total}</strong>${parts.length?` <span class=mut2>(${parts.join(' · ')})</span>`:''}</div>`;
  } else if(d.note){
    disc=`<div class=disc><span class=lab>Discrimination</span><span class=mut2>${d.note}</span></div>`;
  }
  let cm='';
  if(d.cross_model){
    const c=d.cross_model;
    cm = c.verdict==='MATCH'
      ? `<div class=cmline title="Re-tested on a second, independent implementation it was never tuned on, and still recognised the phenomenon.">✓ Generalises to an independent model <span class=mut2>(${c.alt})</span></div>`
      : `<div class=cmline title="On an independent model it ranked its own pattern first but stayed below the firing threshold — a disclosed generalisation limit.">~ Ranks top on an independent model <span class=mut2>(${c.alt}, below threshold)</span></div>`;
  }
  const ge=d.emergence!=null?(+d.emergence).toFixed(2):'—';
  const foot=`<div class=foot><span class=hint title="A coarse, cross-pattern emergence index — NOT this detector's verdict. A definitive match can sit beside a low value here; trust the verdict and the panel.">generic emergence index ${ge}</span>`
    +(d.self_recognized?` · <span class=hint title="Among all 32 detectors run on this pattern, the correct one ranked first. In-sample evidence — the weakest of the three tests.">self-identifies in the 32-detector battery</span>`:'')
    +`</div>`;
  return `<span class=badge style="background:${vc}">${d.verdict}</span><span class=tier title="${TG[d.tier]||''}">${d.tier||'—'}</span>`
    + disc + cm + foot + methLink();
}

function show(i){
  clearTimer(); if(ovl.classList.contains('open')) closeOverlay();
  try{
    const m=M[i], d=m.detector||{};
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
        <button class=dl id=enl type=button title="Open a large, in-page view (supports fullscreen) — no popup window">⛶ Enlarge</button>${m.mp4?`<a class=dl href="assets/${m.mp4}" download="${mp4name}">⤓ MP4</a>`:''}${m.asset?`<a class=dl href="assets/${m.asset}" download="${pngname}" title="One PNG with all ${s.frames} frames tiled in a grid (a sprite sheet, not a single still)">⤓ Frame sheet (PNG)</a>`:''}`;
    } else if(m.asset){
      vizHtml=`<img src="assets/${m.asset}" alt="${m.id} ${esc(m.name)}: ${desc}">`;
    } else { vizHtml='<p style=color:#8b98a5>(no visualization)</p>'; }
    const watchHtml=m.watch?`<div class=watch><b>&#128065; What to look for</b>${m.watch}</div>`:'';
    const det = detectorHtml(d);
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
    sl=document.getElementById('sl'), fc=document.getElementById('fc'), spd=document.getElementById('spd'),
    enl=document.getElementById('enl');
  const BASE=110, rows=Math.ceil(s.frames/s.cols); let cur=0, bigSize=0;
  function paintBig(resize){
    if(resize){
      const avail=Math.min(window.innerWidth-32, window.innerHeight-150);
      bigSize=Math.max(240, Math.min(avail, 1100));
      OVL.stage.style.width=bigSize+'px'; OVL.stage.style.height=bigSize+'px';
      OVL.stage.style.backgroundImage='url(assets/'+m.asset+')';
      OVL.stage.style.backgroundSize=(s.cols*bigSize)+'px '+(rows*bigSize)+'px';
    }
    const c=cur%s.cols, r=Math.floor(cur/s.cols);
    OVL.stage.style.backgroundPosition='-'+(c*bigSize)+'px -'+(r*bigSize)+'px';
  }
  function setF(k){
    cur=((k%s.frames)+s.frames)%s.frames;
    const c=cur%s.cols, r=Math.floor(cur/s.cols);
    stage.style.backgroundPosition=`-${c*s.fw}px -${r*s.fh}px`;
    sl.value=cur; fc.textContent=cur+'/'+(s.frames-1); sl.setAttribute('aria-valuetext',(cur+1)+' of '+s.frames);
    if(ovl.classList.contains('open')){ paintBig(false); OVL.sl.value=cur; OVL.fc.textContent=cur+'/'+(s.frames-1); }
  }
  function setBars(sym,lab){pp.textContent=sym; pp.setAttribute('aria-label',lab); OVL.pp.textContent=sym; OVL.pp.setAttribute('aria-label',lab);}
  function play(){clearTimer(); timer=setInterval(()=>setF(cur+1), BASE/parseFloat(spd.value)); setBars('⏸','Pause animation');}
  function pause(){clearTimer(); setBars('▶','Play animation');}
  const toggle=()=>timer?pause():play();
  pp.onclick=OVL.pp.onclick=toggle;
  rs.onclick=OVL.rs.onclick=()=>{pause();setF(0);};
  sl.oninput=()=>{pause();setF(+sl.value);};
  OVL.sl.oninput=()=>{pause();setF(+OVL.sl.value);};
  spd.onchange=()=>{OVL.spd.value=spd.value; if(timer)play();};
  OVL.spd.onchange=()=>{spd.value=OVL.spd.value; if(timer)play();};
  if(enl) enl.onclick=()=>{
    lastEnl=enl; paintBigRef=paintBig;
    OVL.cap.textContent=m.id+' · '+m.name;
    OVL.stage.setAttribute('aria-label', m.id+' '+m.name+' — enlarged animation');
    OVL.sl.max=s.frames-1; OVL.sl.value=cur; OVL.spd.value=spd.value;
    OVL.fc.textContent=cur+'/'+(s.frames-1);
    OVL.pp.textContent=timer?'⏸':'▶'; OVL.pp.setAttribute('aria-label', timer?'Pause animation':'Play animation');
    ovl.classList.add('open'); ovl.setAttribute('aria-hidden','false');
    paintBig(true); OVL.x.focus();
  };
  setF(0);
  if(RM) pause(); else play();
}

function showAbout(){clearTimer(); setSelected('about'); detail.innerHTML='<div class=about>'+OVERVIEW+'</div>';}
function showMethods(){clearTimer(); setSelected('methods'); detail.innerHTML='<div class=about>'+METHODS+'</div>';}

function renderKey(k){
  k=(k||'').trim(); const kl=k.toLowerCase();
  if(kl==='about'){showAbout();return true;}
  if(kl==='methods'){showMethods();return true;}
  const i=M.findIndex(m=>m.id.toLowerCase()===kl);
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
  e.className='item'+(key==='about'?' about-item':'')+(key==='methods'?' methods-item':'');
  e.dataset.key=String(key); e.dataset.nav=nav;
  if(key==='about') e.id='about-item';
  e.setAttribute('role','option'); e.setAttribute('aria-selected','false'); e.tabIndex=-1;
  e.innerHTML=htmlStr;
  e.onclick=()=>navigate(nav);
  list.appendChild(e);
}
makeItem('about','about','<b>✦</b>About emergence<small>start here — what all 32 share</small>');
makeItem('methods','methods','<b>✓</b>How models are validated<small>what “validated” means here</small>');
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
    "<script>" + JS.replace("__M__", json.dumps(M)).replace("__VC__", json.dumps(VC)).replace("__OVERVIEW__", json.dumps(OVERVIEW_HTML)).replace("__METHODS__", json.dumps(METHODS_HTML)) + "</script></body></html>")
open("gallery/index.html", "w").write(page)
print("wrote gallery/index.html", len(page), "bytes;",
      sum(1 for e in M if e.get("asset_type") == "sprite"), "playable,",
      sum(1 for e in M if e.get("code")), "with code")
