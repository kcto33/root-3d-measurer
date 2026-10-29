// Additive original-input analysis: orientation changes never rewrite measurement.
export async function installMorphology({apiURL, after, initialSettings={}}) {
 const card=document.createElement('section');card.className='card';card.id='morphology_card';
 const heading=document.createElement('h2');heading.textContent='植物形态指标';card.append(heading);
 const controls=document.createElement('div');controls.style.cssText='display:flex;gap:8px;align-items:center;flex-wrap:wrap';
 const axis=document.createElement('select');axis.id='morph_vertical_axis';axis.setAttribute('aria-label','形态指标竖直轴');
 for(const v of ['x','y','z']){const o=document.createElement('option');o.value=v;o.textContent=v.toUpperCase()+' 轴竖直';axis.append(o);}
 axis.value=initialSettings.vertical_axis??'z';
 const sign=document.createElement('select');sign.id='morph_down_sign';sign.setAttribute('aria-label','形态指标向下方向');
 for(const [v,t] of [[-1,'负方向向下'],[1,'正方向向下']]){const o=document.createElement('option');o.value=v;o.textContent=t;sign.append(o);}
 sign.value=initialSettings.down_sign??-1;
 const apply=document.createElement('button');apply.textContent='更新形态指标';
 controls.append(axis,sign,apply);card.append(controls);
 const body=document.createElement('div');card.append(body);after.before(card);
 const format=v=>v==null?'不可用':Number(v).toLocaleString('zh-CN',{maximumSignificantDigits:5});
 const text=(parent,tag,content)=>{const e=document.createElement(tag);e.textContent=content;parent.append(e);return e;};
 const unit=u=>u.replace(/normalized_model_coordinates|normalized_model_units|model_units/g,'模型单位').replace('dimensionless','无量纲').replace('degree','°');
 const status=s=>s.includes('UNAVAILABLE')?'缺少合格依据':s.includes('PARTIAL')?'仅已提取轴':s.includes('LOW_LINEARITY')?'有限尺度／线性较弱':s.includes('INSUFFICIENT')?'尺度范围不足':s.includes('BUDGET')?'计算预算不足':s.includes('FINITE_SCALE')?'有限尺度估计':s.includes('AXIS_REQUIRED')?'缺少中心线':s.includes('INPUT_RESOLUTION')?'缺少分辨率':s.includes('UNQUALIFIED')?'未认证':'几何描述／精度另验';
 function table(parent,headers,rows){const wrap=document.createElement('div');wrap.style.overflowX='auto';const t=document.createElement('table');t.style.cssText='width:100%;border-collapse:collapse;font-size:12px';const head=document.createElement('tr');for(const h of headers)text(head,'th',h).style.cssText='text-align:left;padding:7px;border-bottom:1px solid #425967';t.append(head);for(const row of rows){const tr=document.createElement('tr');for(const v of row)text(tr,'td',v).style.cssText='padding:7px;border-bottom:1px solid #2b414d';t.append(tr);}wrap.append(t);parent.append(wrap);}
 function details(title){const d=document.createElement('details');text(d,'summary',title);body.append(d);return d;}
 async function load(){
  apply.disabled=true;body.replaceChildren();text(body,'p','正在核对原输入并计算形态指标…');
  const query='?'+new URLSearchParams({vertical_axis:axis.value,down_sign:sign.value});
  try{
   const response=await fetch(apiURL('/morphology.json')+query);const value=await response.json();if(!response.ok)throw Error(value.error||response.status);
   body.replaceChildren();const metrics=Object.fromEntries(value.metrics.map(m=>[m.parameter,m]));
   const summary=document.createElement('div');summary.style.cssText='display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:8px;margin-top:12px';body.append(summary);
   for(const key of ['fractal_dimension_3d','root_depth','horizontal_spread','volume_solidity']){const m=metrics[key],cell=document.createElement('div');cell.style.cssText='padding:10px;background:#101e28;border-radius:6px';text(cell,'small',m.label);text(cell,'p',format(m.value)+' '+unit(m.unit)).style.cssText='font-size:20px;margin:6px 0';text(cell,'small',status(m.status));summary.append(cell);}
   const note=text(body,'p',value.settings.vertical_axis.toUpperCase()+' 轴'+(value.settings.down_sign<0?'负':'正')+'方向向下。深度零点取实体最上边界；未指定土壤表面。'+(value.partial_axis?'中心线尚不完整：轴相关指标仅描述当前提取部分，凸包内根长密度不予输出。':''));note.className='muted';
   const fd=value.fractal, f=details('分形维数的尺度与拟合依据');
   text(f,'p','三维中心线盒计数；不是实体体积的维数。使用全部合格二倍尺度、8 种网格偏移，不自动挑选最好看的拟合区间。');
   text(f,'p','有效尺度 '+fd.scales.length+' 级；跨度 '+format(fd.scale_ratio)+' 倍；R² = '+format(fd.r_squared)+'。'+(fd.dimension==null?'当前分辨率或尺度范围不足，未给出维数。':'这是有限尺度估计，拟合优度不能证明具有真实分形性质。'));
   if(fd.phase_dimensions)text(f,'p','各网格偏移的斜率范围：'+format(Math.min(...fd.phase_dimensions))+'–'+format(Math.max(...fd.phase_dimensions))+'。');
   table(f,['盒边长 ('+unit(value.unit)+')','平均盒数','8 种偏移盒数'],fd.scales.map(s=>[format(s.box_size),format(s.mean_count),s.phase_counts.join(' / ')]));
   const all=details('全部形态指标与资格状态（'+value.metrics.length+' 项）');
   table(all,['指标','数值','单位','状态'],value.metrics.map(m=>[m.label,format(m.value),unit(m.unit),status(m.status)]));
   const profile=details('根长与实体体积的深度分布');
   text(profile,'p','10 个等深区间；根长按完整线段在区间内的弧长积分，体积按原始占据单元与区间重叠体积积分。');
   table(profile,['深度区间 ('+unit(value.unit)+')','根长','根长占比','实体体积 ('+unit(value.unit)+'³)'],value.depth_profile.map(r=>[format(r.depth_start)+'–'+format(r.depth_end),format(r.axis_length),r.axis_length_fraction==null?'不可用':format(r.axis_length_fraction*100)+'%',format(r.solid_volume)]));
   const downloads=document.createElement('p');for(const [route,label] of [['morphology','形态指标 CSV'],['morphology_json','形态完整 JSON'],['depth_profile','深度分布 CSV']]){const a=text(downloads,'a',label);a.className='action';a.href=apiURL('/export/'+route)+query;a.download='';downloads.append(' ');}body.append(downloads);
   text(body,'small','以上导出使用当前方向设置。原始参数 CSV 保留测量时设置。分叉为几何图节点，未推断主根／侧根级次、干重或真实土芯根长密度。');
   card.dataset.ready='true';card.dataset.verticalAxis=value.settings.vertical_axis;card.dataset.downSign=value.settings.down_sign;
  }catch(error){body.replaceChildren();text(body,'p','形态指标不可用：'+error.message);card.dataset.ready='false';}
  finally{apply.disabled=false;}
 }
 apply.onclick=load;await load();
}
