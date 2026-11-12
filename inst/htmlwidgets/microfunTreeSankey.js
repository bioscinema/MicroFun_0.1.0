HTMLWidgets.widget({name:"microfunTreeSankey",type:"output",factory:function(el,width,height){return {renderValue:function(x){el.innerHTML=x.svg;el.style.overflow="auto";var svg=el.querySelector("svg");svg.style.maxWidth="100%";svg.style.height="auto";(function(svg){
  if(!svg)return;
  const byTaxon=new Map(),byFunction=new Map();
  const add=(m,k,p)=>{if(!m.has(k))m.set(k,[]);m.get(k).push(p);};
  svg.querySelectorAll('.mf-link').forEach(p=>{add(byTaxon,p.dataset.taxon,p);add(byFunction,p.dataset.function,p);});
  let active=[],current=null;
  const clear=()=>{active.forEach(p=>p.classList.remove('mf-active'));active=[];svg.classList.remove('mf-hover');current=null;};
  svg.addEventListener('mouseover',ev=>{
    const el=ev.target.closest('.mf-node,.mf-link');
    if(!el||el===current)return;
    clear();current=el;
    const t=el.dataset.taxon,f=el.dataset.function;
    active=(t?byTaxon.get(t):byFunction.get(f))||[];
    if(t&&f)active=active.filter(p=>p.dataset.function===f);
    svg.classList.add('mf-hover');active.forEach(p=>p.classList.add('mf-active'));
  });
  svg.addEventListener('mouseout',ev=>{if(current&&!current.contains(ev.relatedTarget))clear();});
  svg.addEventListener('mouseleave',clear);
})(svg);},resize:function(width,height){}};}});
