"""Delivery annotations and read-only visualization; no geometry is remeasured."""
import hashlib,json,os
from pathlib import Path
from root_graph_source_length import q

def delivery_annotations(report):
    graph=report.get('candidate_centerline');profile=report.get('radius_curve');centering=[]
    if profile:
        for index,row in enumerate(profile['rows']):
            section=row['original_mesh_section'];point=section.get('exact_query_point');centroid=section.get('exact_centroid')
            equal=(tuple(map(q,point))==tuple(map(q,centroid))) if point and centroid else None
            centering.append(dict(radius_row_index=index,generated_edge_id=row['generated_edge_id'],status='MEASURED_SECTION_CENTROID_MATCH' if equal else 'OFF_MEASURED_SECTION_CENTROID' if equal is False else 'CENTRING_UNKNOWN',exact_query_equals_area_centroid=equal))
    counts={s:sum(x['status']==s for x in centering) for s in ('MEASURED_SECTION_CENTROID_MATCH','OFF_MEASURED_SECTION_CENTROID','CENTRING_UNKNOWN')}
    return dict(status=report['status'],automatic_nodes=len(graph['points']) if graph else 0,path_count=len(profile['paths']) if profile else 0,
        candidate_length=report.get('length'),Source_radius_stations=profile['station_denominator'] if profile else 0,available_Source_radii=profile['available_radii'] if profile else 0,
        failure_indices=list(range(len(report.get('failures',[])))),failure_location_count=len(report.get('failures',[])),all_original_failure_locations_referenced=True,
        centering=centering,centering_counts=counts,centering_definition='exact query equals the exact area centroid of this measured section only; no medial/anatomical/globalcurve assertion',
        no_model_specific_node_ids_required=report.get('new_model_without_manual_node_ids',False),original_model_modified=False,geometry_or_section_kernels_reexecuted=0,
        first_version_functional_outputs_available=bool(graph and report.get('length') and profile),remaining_blockers=['centering unmet/unknown portions explicitly marked','endpoint and continuous Source containment remain unqualified','values between radius sample stations remain unmeasured'],
        physical_or_biological_admission=False)

def write_delivery(report_path,output_prefix,experiment_root,*,model_path=None,case_label=None):
    report_path=Path(report_path);prefix=Path(output_prefix);root=Path(experiment_root)
    report=json.loads(report_path.read_text(encoding='utf-8-sig'));raw=report_path.read_bytes();annotation=delivery_annotations(report)
    data=dict(schema_version='root-model-integrated-delivery-v001',measurement=dict(experiment_relative_path=report_path.resolve().relative_to(root.resolve()).as_posix(),sha256=hashlib.sha256(raw).hexdigest(),bytes=len(raw)),summary=annotation,
        viewer=dict(type='read-only HTTP canvas orbit/radius/failure viewer',large_source_mesh_or_graph_copied=False),original_source_report_preserved=True,case_label=case_label)
    if model_path is not None:
        model=Path(model_path).resolve();model_raw=model.read_bytes();binding=report['input_file_binding']
        if hashlib.sha256(model_raw).hexdigest()!=binding['sha256'] or len(model_raw)!=binding['bytes']:raise ValueError('ORIGINAL_MODEL_OVERLAY_SHA_BYTES_MISMATCH')
        data['original_model']=dict(relative_to_delivery=os.path.relpath(model,prefix.parent),sha256=binding['sha256'],bytes=binding['bytes'],preview_byte_guard=2000000,overlay_uses_original_coordinates=True)
    path=prefix.with_name(prefix.name+'_delivery.json');html=prefix.with_name(prefix.name+'_viewer.html')
    with path.open('x',encoding='utf8') as f:json.dump(data,f,indent=2,allow_nan=False)
    with html.open('x',encoding='utf8') as f:f.write(VIEWER)
    return data,path,html

VIEWER='''<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>模型测量交付</title>
<style>*{box-sizing:border-box}body{margin:0;background:#101820;color:#e3edf3;font:14px system-ui}header{padding:18px 24px;border-bottom:1px solid #324452}h1{font-size:21px;margin:0 0 7px}small,.muted{color:#a8bdc9}main{display:grid;grid-template-columns:minmax(400px,1fr) 360px;gap:12px;padding:12px}.card{background:#17242e;border:1px solid #304453;border-radius:8px;padding:14px}canvas{display:block;width:100%;height:460px;background:#101b23;border-radius:5px;touch-action:none}#chart{height:210px}#stats{display:flex;gap:25px;flex-wrap:wrap}.stat{font-size:25px}.stat small{font-size:12px;display:block}button,select{background:#233c4b;border:1px solid #577180;color:#fff;border-radius:5px;padding:6px;cursor:pointer}#issues{max-height:565px;overflow:auto}.issue{display:block;width:100%;text-align:left;margin:6px 0;font-size:11px;word-break:break-word}.orange{color:#ffb65c}.green{color:#65d7b0}.red{color:#ff7e85}#detail{white-space:pre-wrap;font:12px ui-monospace;max-height:220px;overflow:auto}label{margin-right:14px}footer{padding:0 24px 18px}@media(max-width:900px){main{display:block}.card{margin-bottom:12px}canvas{height:360px}}</style>
<header><h1>新模型 → 候选中心线 / 路径长度 / 采样半径 / 失败位置</h1><small id="binding">正在读取结果；请通过项目只读查看服务打开本页。</small></header>
<main><section><div class="card" id="stats"></div><div class="card"><p><span class="green">● 采样截面质心重合</span>　<span class="orange">● 未满足截面居中</span>　<span class="red">● 缺失/失败</span>　<span class="muted">● 未核验</span></p>
<p><label><input type="checkbox" id="surface" checked>叠加原模型</label><label><input type="checkbox" id="rings">面积等效半径环</label><button id="reset">复位视角</button></p><canvas id="scene"></canvas><small id="surface-note">拖动旋转，滚轮缩放，点击采样点查看。原模型按原始坐标叠加；半径环仅为面积等效示意。</small></div>
<div class="card"><h3>沿线原三角面采样半径</h3><select id="path"></select><canvas id="chart"></canvas><small>横轴：候选路径弧长；纵轴：采样面积等效半径。只画实测点，缺失不补值、不插值。</small></div></section>
<aside><div class="card"><h3>选中位置</h3><div id="detail">点击采样点或失败项。</div></div><div class="card"><h3>全部缺失与未满足项</h3><div id="issues"></div></div></aside></main>
<footer class="muted">中心线与长度是有声明分辨率的几何候选结果。居中标记只比较已测截面的精确面积质心；不是完整轴、物理尺度或真实准确率资格。</footer>
<script>
let report,delivery,graph,modelPreview=null,rows=[],yaw=.7,pitch=-.45,zoom=1,selected=null,screen=[];const $=s=>document.querySelector(s),scene=$('#scene'),chart=$('#chart');
const big=x=>x.startsWith('-')?-BigInt(x.slice(1)):BigInt(x);const q=x=>Number(big(x[0]))/Number(big(x[1]));const vec=x=>x.map(q);function detail(x){$('#detail').textContent=JSON.stringify(x,null,2)}
function setup(c){const d=window.devicePixelRatio||1,r=c.getBoundingClientRect();c.width=r.width*d;c.height=r.height*d;const ctx=c.getContext('2d');ctx.setTransform(d,0,0,d,0,0);return [ctx,r.width,r.height]}
let center=[0,0,0],scale=1,bounds;function project(p,w,h){let [x,y,z]=p.map((v,i)=>(v-center[i])*scale*zoom),a=x*Math.cos(yaw)+z*Math.sin(yaw),b=-x*Math.sin(yaw)+z*Math.cos(yaw),v=y*Math.cos(pitch)-b*Math.sin(pitch);return [w/2+a,h/2-v]}
function line(ctx,a,b,color,width=1){ctx.strokeStyle=color;ctx.lineWidth=width;ctx.beginPath();ctx.moveTo(...a);ctx.lineTo(...b);ctx.stroke()}
function draw(){if(!report)return;let [ctx,w,h]=setup(scene);ctx.clearRect(0,0,w,h);screen=[];let pts=graph?graph.points:[],adj=graph?graph.adjacency:[],colors=new Map();delivery.summary.centering.forEach(x=>colors.set(x.generated_edge_id,x.status));let eid=0;
if(modelPreview?.available&&$('#surface').checked){const vs=modelPreview.vertices,ps=vs.map(p=>project(p,w,h));ctx.fillStyle='rgba(133,166,184,0.085)';ctx.strokeStyle='rgba(146,178,197,0.24)';ctx.lineWidth=.6;for(const face of modelPreview.faces){ctx.beginPath();ctx.moveTo(...ps[face[0]]);ctx.lineTo(...ps[face[1]]);ctx.lineTo(...ps[face[2]]);ctx.closePath();ctx.fill();ctx.stroke()}}
if(bounds){let corners=[];for(let i=0;i<8;i++)corners.push([0,1,2].map(k=>i&(1<<k)?bounds.maximum[k]:bounds.minimum[k]));for(let i=0;i<8;i++)for(let k=0;k<3;k++){let j=i^(1<<k);if(i<j)line(ctx,project(corners[i],w,h),project(corners[j],w,h),'#34505e')}}
for(let a=0;a<adj.length;a++)for(const b of [...adj[a]].sort((x,y)=>x-y))if(a<b){const s=colors.get(eid++),c=s==='MEASURED_SECTION_CENTROID_MATCH'?'#65d7b0':s==='OFF_MEASURED_SECTION_CENTROID'?'#ffb65c':s==='CENTRING_UNKNOWN'?'#ff7e85':'#78919e';line(ctx,project(pts[a],w,h),project(pts[b],w,h),c,2)}
rows.forEach((r,i)=>{let p=project(r.position,w,h),s=delivery.summary.centering[i]?.status,c=r.available?(s==='MEASURED_SECTION_CENTROID_MATCH'?'#65d7b0':'#ffb65c'):'#ff7e85';ctx.fillStyle=c;ctx.beginPath();ctx.arc(...p,selected===i?6:3.5,0,Math.PI*2);ctx.fill();screen.push({x:p[0],y:p[1],i});if($('#rings').checked&&r.available){let n=vec(r.exact_normal),norm=Math.hypot(...n);n=n.map(x=>x/norm);let a=Math.abs(n[0])<.8?[1,0,0]:[0,1,0],u=[n[1]*a[2]-n[2]*a[1],n[2]*a[0]-n[0]*a[2],n[0]*a[1]-n[1]*a[0]],l=Math.hypot(...u);u=u.map(x=>x/l);let v=[n[1]*u[2]-n[2]*u[1],n[2]*u[0]-n[0]*u[2],n[0]*u[1]-n[1]*u[0]],prev=null;for(let j=0;j<=32;j++){let t=j*Math.PI/16,x=r.position.map((x,k)=>x+r.area_equivalent_radius*(u[k]*Math.cos(t)+v[k]*Math.sin(t))),cur=project(x,w,h);if(prev)line(ctx,prev,cur,c,1);prev=cur}}});
report.failures.filter(x=>x.position).forEach(x=>{let p=project(x.position,w,h);line(ctx,[p[0]-3,p[1]-3],[p[0]+3,p[1]+3],'#ed6572');line(ctx,[p[0]-3,p[1]+3],[p[0]+3,p[1]-3],'#ed6572')});}
function plot(){let [ctx,w,h]=setup(chart);ctx.clearRect(0,0,w,h);const paths=report.radius_curve?.paths||[],p=paths[Number($('#path').value)||0];if(!p)return;let samples=p.samples,xmax=Math.max(1e-99,...samples.map(s=>s.arc_position)),ymax=Math.max(1e-99,...samples.map(s=>s.radius||0))*1.12;line(ctx,[45,15],[45,h-30],'#75929f');line(ctx,[45,h-30],[w-15,h-30],'#75929f');ctx.fillStyle='#bcd0da';ctx.font='11px system-ui';ctx.fillText('0',25,h-15);ctx.fillText(xmax.toPrecision(5),w-85,h-10);ctx.fillText(ymax.toPrecision(4),5,12);samples.forEach(s=>{let x=45+s.arc_position/xmax*(w-65);ctx.fillStyle=s.radius===null?'#ff7e85':'#64d8ce';ctx.beginPath();ctx.arc(x,s.radius===null?h-30:h-30-s.radius/ymax*(h-50),4,0,Math.PI*2);ctx.fill()})}
Promise.all([fetch('measurement.json').then(r=>{if(!r.ok)throw Error(r.status);return r.json()}),fetch('delivery.json').then(r=>r.json()),fetch('model_preview.json').then(r=>r.ok?r.json():({available:false,status:'ORIGINAL_MODEL_PREVIEW_UNAVAILABLE'})).catch(()=>({available:false}))]).then(([a,b,m])=>{report=a;delivery=b;modelPreview=m;graph=a.candidate_centerline;rows=a.radius_curve?.rows||[];bounds=a.automatic_generation_evidence?.model_bounds||a.model_bounds;
if(bounds){center=bounds.minimum.map((x,i)=>(x+bounds.maximum[i])/2);scale=300/Math.max(...bounds.maximum.map((x,i)=>x-bounds.minimum[i]),1e-99)}
$('#binding').textContent=(b.case_label?b.case_label+' | ':'')+'输入 SHA '+a.input_file_binding.sha256.slice(0,16)+'… | '+a.coordinate_units+' | '+a.status;const s=b.summary;if(!m.available)$('#surface-note').textContent='原模型预览未提供或受guard限制；灰框仅为包围盒。中心线和数值保持可查看。';
for(const [value,label] of [[s.automatic_nodes,'自动节点'],[s.path_count,'路径'],[s.candidate_length?.value?.toPrecision(7)||'缺失','候选长度'],[s.available_Source_radii+'/'+s.Source_radius_stations,'原三角面半径采样'],[s.centering_counts.OFF_MEASURED_SECTION_CENTROID,'未满足截面居中']]){let e=document.createElement('div');e.className='stat';e.textContent=value;let t=document.createElement('small');t.textContent=label;e.append(t);$('#stats').append(e)}
(a.radius_curve?.paths||[]).forEach((p,i)=>{let o=document.createElement('option');o.value=i;o.textContent=p.path_id+' | 长度 '+p.sampled_length.toPrecision(6);$('#path').append(o)});
let issues=$('#issues');a.failures.forEach((x,i)=>{let e=document.createElement('button');e.className='issue red';e.textContent=(i+1)+'. '+x.code+(x.position?' @ '+x.position.map(v=>v.toPrecision(5)).join(', '):'（模型级）');e.onclick=()=>detail(x);issues.append(e)});
b.summary.centering.filter(x=>x.status!=='MEASURED_SECTION_CENTROID_MATCH').forEach(x=>{let e=document.createElement('button');e.className='issue orange';e.textContent='采样边 '+x.generated_edge_id+': '+x.status;e.onclick=()=>{selected=x.radius_row_index;detail({centering:x,station:rows[selected]});draw()};issues.append(e)});draw();plot();
}).catch(e=>{$('#binding').textContent='读取失败：'+e.message+'。请运行 serve_model_delivery.py --prefix <交付前缀>，通过 http://127.0.0.1 打开；不要直接用 file://。'});
let drag=null;scene.onpointerdown=e=>{drag=[e.clientX,e.clientY];scene.setPointerCapture(e.pointerId)};scene.onpointermove=e=>{if(!drag)return;yaw+=(e.clientX-drag[0])*.01;pitch=Math.max(-1.5,Math.min(1.5,pitch+(e.clientY-drag[1])*.01));drag=[e.clientX,e.clientY];draw()};scene.onpointerup=()=>drag=null;scene.onwheel=e=>{e.preventDefault();zoom=Math.max(.1,Math.min(20,zoom*Math.exp(-e.deltaY*.001)));draw()};scene.onclick=e=>{let r=scene.getBoundingClientRect(),x=e.clientX-r.left,y=e.clientY-r.top,best=screen.map(p=>({...p,d:Math.hypot(p.x-x,p.y-y)})).sort((a,b)=>a.d-b.d)[0];if(best&&best.d<15){selected=best.i;detail({centering:delivery.summary.centering[selected],station:rows[selected]});draw()}};$('#reset').onclick=()=>{yaw=.7;pitch=-.45;zoom=1;draw()};$('#rings').onchange=draw;$('#surface').onchange=draw;$('#path').onchange=plot;window.onresize=()=>{draw();plot()};
</script></html>'''
