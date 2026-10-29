// Additive phenotype analysis and measurement-bound human annotations.
// All model/user strings enter the DOM through textContent, never HTML.
export async function installPhenotyping({apiURL, report, graph, paths, scene, project, highlight, after}) {
  const existing = document.getElementById('phenotyping_card');
  if (existing) return;
  const card = document.createElement('section');
  card.className = 'card';
  card.id = 'phenotyping_card';
  card.innerHTML = `
    <style>
      #phenotyping_card{font-size:13px;line-height:1.55}
      #phenotyping_card h3{font-size:14px;margin:15px 0 8px}
      #phenotyping_card p{margin:8px 0}
      #phenotyping_card details{border-top:1px solid #344b58;padding-top:10px;margin-top:12px}
      #phenotyping_card summary{cursor:pointer;font-weight:600;color:#dcebf1}
      #phenotyping_card label{display:block;margin:7px 0 3px;font-size:12px}
      #phenotyping_card input,#phenotyping_card select{display:block;width:100%;max-width:100%;background:#10202a;color:#e6f1f5;border:1px solid #507185;border-radius:5px;padding:7px;margin:3px 0 7px;font:inherit}
      #phenotyping_card input:focus,#phenotyping_card select:focus,#phenotyping_card button:focus-visible{outline:2px solid #74ddc8;outline-offset:2px}
      #phenotyping_card button{font:inherit;line-height:1.35;margin:2px 0}
      #phenotyping_card button:disabled{cursor:default;opacity:.45}
      #phenotyping_card .ph-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(125px,1fr));gap:5px 10px}
      #phenotyping_card .ph-actions{display:flex;flex-wrap:wrap;gap:6px;margin:9px 0}
      #phenotyping_card .ph-primary{background:#216153;border-color:#68c3ac}
      #phenotyping_card .ph-note{color:#a9bdc8;font-size:12px}
      #phenotyping_card .ph-status{padding:9px 10px;border-radius:5px;background:#10212b;border-left:3px solid #76c5b4;overflow-wrap:anywhere}
      #phenotyping_card .ph-error{border-left-color:#ef927d;color:#ffd0bd}
      #phenotyping_card .ph-warning{border-left-color:#e2ae67;color:#f4cf9b}
      #phenotyping_card .ph-chip{display:flex;gap:8px;align-items:center;justify-content:space-between;padding:7px;background:#12212c;border-radius:5px;margin:5px 0}
      #phenotyping_card .ph-chip span{overflow-wrap:anywhere}
      #phenotyping_card .ph-stat{padding:9px;background:#11212b;border-radius:5px;min-width:0}
      #phenotyping_card .ph-stat strong{display:block;font-size:19px;font-weight:500;color:#8cdec9;overflow-wrap:anywhere}
      #phenotyping_card .ph-stat small{display:block;font-size:11px}
      #phenotyping_card .ph-table{overflow:auto;max-height:330px;margin-top:8px}
      #phenotyping_card table{width:100%;border-collapse:collapse;font-size:11px;white-space:nowrap}
      #phenotyping_card th,#phenotyping_card td{text-align:left;padding:6px 7px;border-bottom:1px solid #2e4451}
      #phenotyping_card th{color:#bdd0dc;position:sticky;top:0;background:#192c37;z-index:1}
      #phenotyping_card .ph-bars{margin:10px 0}
      #phenotyping_card .ph-bar{display:grid;grid-template-columns:minmax(55px,90px) 1fr minmax(45px,68px);gap:7px;align-items:center;font-size:11px;margin:6px 0}
      #phenotyping_card .ph-track{height:9px;background:#0e1c25;border-radius:3px;overflow:hidden}
      #phenotyping_card .ph-fill{height:100%;background:#70cdb6;border-radius:3px}
      #phenotyping_card .ph-junction{padding:8px;background:#12232d;border-radius:6px;margin:7px 0}
      #phenotyping_card .ph-selected{color:#8de3c8}
      #phenotyping_card canvas{width:100%;height:auto;border:1px solid #344b58;border-radius:5px;margin-top:8px;background:#0e1d26}
      #phenotyping_card .ph-nowrap{white-space:nowrap}
    </style>
    <h2>根系表型与根轴标注</h2>
    <p class="ph-note">先查看直径与方向分布，再指定根基并核对根轴。候选结构、人工确认和几何精度验收分别记录。</p>
    <p id="ph_status" class="ph-status" role="status" aria-live="polite">正在载入当前测量的表型分析…</p>
    <div id="ph_stats" class="ph-grid"></div>
    <div class="ph-actions"><button id="ph_analyze" class="ph-primary">应用设置并分析</button><button id="ph_save">保存分析与标注</button><button id="ph_reload">恢复已保存版本</button></div>
    <details id="ph_settings"><summary>分箱、坐标方向与测量窗口</summary>
      <p id="ph_units" class="ph-note"></p>
      <div class="ph-grid">
        <div><label for="ph_vertical_axis">竖直轴</label><select id="ph_vertical_axis"><option value="x">X 轴</option><option value="y">Y 轴</option><option value="z">Z 轴</option></select></div>
        <div><label for="ph_down_sign">向下方向</label><select id="ph_down_sign"><option value="-1">负方向向下</option><option value="1">正方向向下</option></select></div>
      </div>
      <label for="ph_diameter_edges">直径内部分界值</label><input id="ph_diameter_edges" type="text" inputmode="decimal" placeholder="0.5, 1, 2">
      <p class="ph-note">自动包含 0 和无上界的末箱；分界值仅表示数值类别。</p>
      <div class="ph-grid"><div><label for="ph_depth_bins">等深层数</label><input id="ph_depth_bins" type="number" min="1" max="512" step="1"></div><div><label for="ph_depth_edges">或指定完整深度边界</label><input id="ph_depth_edges" type="text" placeholder="留空使用等深层"></div></div>
      <label for="ph_inclination_edges_deg">竖直夹角完整边界（0–90°）</label><input id="ph_inclination_edges_deg" type="text">
      <label for="ph_azimuth_edges_deg">无向根段朝向完整边界（0–180°）</label><input id="ph_azimuth_edges_deg" type="text">
      <label for="ph_bearing_edges_deg">位置方位完整边界（0–360°）</label><input id="ph_bearing_edges_deg" type="text">
      <label for="ph_origin">位置方位与径向分布参考点 X, Y, Z（可选）</label><input id="ph_origin" type="text" placeholder="留空时不计算相对位置分布">
      <div class="ph-actions"><button id="ph_origin_from_base">用当前根基节点设参考点</button></div>
      <label for="ph_radial_edges">水平径向距离内部分界值</label><input id="ph_radial_edges" type="text">
      <div class="ph-grid">
        <div><label for="ph_tangent_window">局部角度窗口</label><input id="ph_tangent_window" type="number" min="0" step="any"></div>
        <div><label for="ph_diameter_offset">距分叉粗度避让距离</label><input id="ph_diameter_offset" type="number" min="0" step="any"></div>
        <div><label for="ph_diameter_window">粗度统计窗口</label><input id="ph_diameter_window" type="number" min="0" step="any"></div>
        <div><label for="ph_minimum_window_fraction">最小窗口保留比例</label><input id="ph_minimum_window_fraction" type="number" min="0" max="1" step="0.05"></div>
        <div><label for="ph_minimum_radius_coverage">最小有效半径覆盖率</label><input id="ph_minimum_radius_coverage" type="number" min="0" max="1" step="0.05"></div>
      </div>
      <p class="ph-note">深度零点来自原实体上边界；只有轴数据时使用已提取轴的上边界。局部父子夹角和整根方向分别计算，缺失半径保留为空。</p>
    </details>
    <details id="ph_annotations" open><summary>根基、续接与人工确认</summary>
      <div class="ph-grid"><div><label for="ph_base_node">根基节点编号</label><input id="ph_base_node" type="number" min="0" step="1"></div><div><label for="ph_base_type">根基类型标签（可选）</label><input id="ph_base_type" type="text" maxlength="200" placeholder="例如 primary / nodal"></div></div>
      <div class="ph-actions"><button id="ph_pick" aria-pressed="false">在模型上选节点</button><button id="ph_add_base">添加根基并分析</button></div>
      <p id="ph_pick_note" class="ph-note">输入原图节点编号，或开启选点后在模型中左键点击；选中后仍需点击“添加根基”。</p>
      <div id="ph_bases"></div>
      <h3>连通分量</h3><div id="ph_components"></div>
      <p class="ph-note">“确认此分量”表示你已核对根基、分叉续接与侧根归属；确认不等于精度验收。修改根基或续接会清除旧根轴的确认与类型标签；修改端点类型会清除确认。</p>
      <details id="ph_junction_section"><summary>核对分叉续接</summary><div id="ph_junctions"></div></details>
      <label for="ph_axis_select">选择根轴并在模型上高亮</label><select id="ph_axis_select"></select>
      <div class="ph-grid"><div><label for="ph_axis_type">所选根轴类型</label><input id="ph_axis_type" type="text" maxlength="200" placeholder="unknown"></div><div class="ph-actions"><button id="ph_set_type">更新类型标签</button></div></div>
      <div class="ph-grid"><div><label for="ph_endpoint_node">端点节点编号</label><input id="ph_endpoint_node" type="number" min="0" step="1"></div><div><label for="ph_endpoint_type">端点类型</label><select id="ph_endpoint_type"><option value="unknown">未知端点</option><option value="root_tip">已确认生物根尖</option><option value="cut_end">切断／截断端</option><option value="root_base">已指定根基</option></select></div></div>
      <div class="ph-actions"><button id="ph_use_apex">填入所选轴末端</button><button id="ph_set_endpoint">更新端点类型</button><button id="ph_clear_highlight">清除高亮</button></div>
    </details>
    <details id="ph_distribution_section" open><summary>直径、方向与空间分布</summary>
      <div id="ph_distribution_note" class="ph-note"></div><div id="ph_diameter_chart"></div>
      <details><summary>竖直夹角与无向根段朝向</summary><div id="ph_direction_charts"></div></details>
      <details><summary>直径 × 深度、朝向与位置 × 深度</summary><div id="ph_joint_tables"></div></details>
      <details><summary>相对参考点的位置方位与径向分布</summary><div id="ph_bearing_chart"></div><div id="ph_radial_chart"></div></details>
    </details>
    <details id="ph_axis_section" open><summary>所选根轴的分枝、粗度与渐细</summary><div id="ph_axis_details"></div></details>
    <details id="ph_relation_section"><summary>父子根夹角与直径关系</summary><div id="ph_relations"></div></details>
    <details id="ph_summary_section"><summary>按根序、类型与确认状态汇总</summary><div id="ph_summaries"></div></details>
    <details id="ph_export_section"><summary>导入标注与导出结果</summary>
      <p class="ph-note">标注文件绑定当前测量的 SHA-256。导入会检查来源；导出使用当前已分析的设置。</p>
      <div class="ph-actions"><button id="ph_download_json">完整分析 JSON</button><button id="ph_download_annotations">标注 JSON</button><button id="ph_import">导入标注 JSON</button></div>
      <input id="ph_import_file" type="file" accept=".json,application/json" hidden style="display:none">
      <div id="ph_csv_buttons" class="ph-actions"></div><p id="ph_binding" class="ph-note"></p>
    </details>`;
  after.before(card);
  const $ = selector => card.querySelector(selector);
  const defaults = {vertical_axis:'z', down_sign:-1, diameter_edges:[.5,1,2], depth_bins:10,
    depth_edges:null, inclination_edges_deg:[0,15,30,45,60,75,90], azimuth_edges_deg:[0,15,30,45,60,75,90,105,120,135,150,165,180], bearing_edges_deg:[0,45,90,135,180,225,270,315,360],
    radial_edges:[1,2,5,10], origin:null, tangent_window:1, diameter_offset:1, diameter_window:1,
    minimum_window_fraction:.5, minimum_radius_coverage:.5};
  const clone = value => JSON.parse(JSON.stringify(value));
  const emptyAnnotations = () => ({schema_version:'root-axis-annotations-v001', root_bases:[],
    junction_continuations:[], endpoint_types:[], root_types:[], confirmed_axes:[], confirmed_components:[]});
  let token = null, result = null, annotations = emptyAnnotations(), measurementHash = null;
  let dirty = false, editVersion = 0, operationId = 0, busy = false, activeAction = null, controller = null;
  let selectedAxis = null, pickMode = false, pickDown = null;
  const finite = value => typeof value === 'number' && Number.isFinite(value);
  const format = value => finite(value) ? value!==0&&(Math.abs(value)<1e-4||Math.abs(value)>=1e8)?value.toExponential(4):value.toLocaleString('zh-CN', {maximumSignificantDigits:5}) : '不可用';
  const percent = value => finite(value) ? (100 * value).toFixed(1) + '%' : '不可用';
  const unitText = value => String(value || '').replace(/normalized_model_coordinates|normalized_model_units|model_units|model_unit/g, '模型单位');
  const unit = () => unitText(result?.distributions?.unit || (report.scale_mm_per_model_unit == null ? report.coordinate_units : 'mm'));
  const add = (parent, tag, content, className) => {const node=document.createElement(tag); if(content!=null)node.textContent=String(content);if(className)node.className=className;parent.append(node);return node;};
  const button = (parent, label, action, blocked=false, mutates=false) => {const node=add(parent,'button',label);node.type='button';node.dataset.phBlocked=String(blocked);node.dataset.phNeedsToken=String(mutates);node.disabled=busy||blocked||mutates&&!token;node.onclick=()=>guard(action);return node;};
  function statusLabel(status) {
    const exact = {GEOMETRY_ONLY:'仅几何路径／未定根轴', ROOT_BASE_REQUIRED:'尚未指定根基',
      ROOTED_CANDIDATE:'候选结构，待核对', ROOTED_CONFIRMED:'人工确认', ROOTED_PARTIALLY_CONFIRMED:'部分根轴已确认',
      ROOTED_AUTO_INFERRED:'多证据自动推断（待生物学验证）', ROOTED_REVIEW_RECOMMENDED:'部分续接证据不足，建议核对',
      ROOTED_UNRESOLVED_CONTINUATION:'分叉续接有歧义', BLOCKED_CYCLIC_COMPONENT:'保留闭环，未推断层级',
      PARTIAL_ROOTED_CYCLIC_COMPONENT:'保留环内区域，分析环外有向根段', PARTIAL_ROOTED_ARCHITECTURE:'含环局部根段分析',
      CYCLIC_CORE_INTERRUPTS_COMPLETE_MOTHER_AXIS:'环内区域中断完整母根，长度未定',
      LOCAL_FRAGMENT_CANDIDATE_RELATIONSHIPS:'局部根段候选关系，非整根统计',
      BLOCKED_MULTIPLE_BASES_IN_COMPONENT:'同一分量有多个竞争根基', NO_EDGES:'没有可测边',
      CONFIRMED_ARCHITECTURE:'人工确认的根轴解释', CANDIDATE_ARCHITECTURE:'候选根轴解释',
      AUTO_INFERRED_ARCHITECTURE:'多证据自动推断的根轴解释',
      UNRESOLVED_ARCHITECTURE:'根轴解释仍有歧义', PARTIAL_ARCHITECTURE:'部分结构可解释',
      CONFIRMED:'人工确认', AUTO_INFERRED:'自动推断／非人工确认', CANDIDATE:'候选', UNRESOLVED:'待消除歧义',
      ANNOTATION_CONFIRMED:'人工确认／精度另验', AUTO_INFERRED_AXIS_RELATIONSHIPS:'自动推断父子关系', CANDIDATE_AXIS_RELATIONSHIPS:'候选父子关系',
      EXPLICIT_CONTINUATION:'已指定续接', EXPLICIT_TERMINATION:'已指定本根终止',
      INFERRED_STRAIGHTEST_CANDIDATE:'自动直顺候选', AUTO_INFERRED_CONTINUATION:'多证据自动推断',
      INFERRED_MULTI_CUE_CANDIDATE:'多证据候选，建议复核', UNRESOLVED_CONTINUATION:'未确定续接',
      HIGH_GEOMETRIC:'高几何证据／未验证', REVIEW_RECOMMENDED:'建议人工核对',
      UNRESOLVED_AXIS_CONTINUATION:'根轴续接未决',
      SAMPLE_SUPPORTED_ESTIMATE:'采样支持估计', OBSERVED_AXIS_DESCRIPTOR:'已提取轴的几何分布',
      TRUE_ROOT_TIP_NOT_IDENTIFIED:'未确认真实根尖', NATURAL_AXIS_BASE_NOT_IDENTIFIED:'未确认自然根基',
      AT_LEAST_TWO_LATERALS_REQUIRED:'不足两条侧根', NO_LATERAL_ROOTS:'没有已识别侧根',
      NO_POSITIVE_BRANCHED_ZONE:'没有正长度的分枝区', ORIENTED_OPEN_AXIS_REQUIRED:'需要有方向的开放根轴',
      AT_LEAST_THREE_DISTINCT_RADIUS_STATIONS_REQUIRED:'不足三个不同位置的有效半径样本',
      UNAVAILABLE_EXPLICIT_RADIAL_ORIGIN_REQUIRED:'请先指定径向参考点', UNAVAILABLE_ZERO_DEPTH_SPAN:'当前方向没有深度跨度',
      UNAVAILABLE_EXPLICIT_BEARING_ORIGIN_REQUIRED:'请先指定位置方位参考点',
      UNAVAILABLE_NO_QUALIFIED_RADIUS:'没有合格单分支半径', UNAVAILABLE_NO_MEASURABLE_AXIS:'没有可测中心线'};
    if (exact[status]) return exact[status];
    if (!status) return '不可用';
    if (status.includes('PARTIAL')) return '部分弧段／按覆盖率解释';
    if (status.includes('CANDIDATE')) return '候选解释';
    if (status.includes('CONFIRMED')) return '人工确认／精度另验';
    if (status==='AVAILABLE') return '可用／精度另验';
    if (status.includes('UNAVAILABLE')||status.includes('REQUIRED')||status.includes('INSUFFICIENT')||status.includes('AMBIGUOUS')||status.startsWith('NO_')) return '缺少合格依据';
    return '请查看完整记录';
  }
  const errors = {UNRESOLVED_CONTINUATION_REQUIRES_EXPLICIT_OVERRIDE:'该分量仍有歧义，请先逐个指定分叉续接。',
    COMPONENT_CONFIRMATION_REQUIRES_UNIQUE_ROOTED_TREE:'仅能确认已有唯一根基、无闭环且续接明确的分量。',
    CONTINUATION_CONFLICTS_WITH_ROOTED_DIRECTION:'续接与当前根基方向冲突，请重新核对该分叉。',
    UNKNOWN_AXIS_ID_OR_CHANGED_AXIS_PARTITION:'根轴划分已改变，请清除旧确认或重新导入匹配的标注。',
    ROOT_TIP_OR_CUT_REQUIRES_NONBASE_DEGREE_ONE_NODE:'生物根尖或截断端必须选择非根基的图端点。',
    ROOT_BASE_LABEL_REQUIRES_ROOT_BASE_ANNOTATION:'请先把这个节点添加为根基。'};
  function message(text, kind='') {const p=$('#ph_status');p.textContent=text;p.className='ph-status'+(kind?' ph-'+kind:'');}
  function guard(action) {try {const value=action();if(value?.catch)value.catch(error=>message(errors[error.message]||error.message,'error'));}catch(error){message(errors[error.message]||error.message,'error');}}
  function updateButtons() {
    for (const node of card.querySelectorAll('button')) node.disabled = busy || node.dataset.phBlocked==='true' || !token&&node.dataset.phNeedsToken==='true';
    for (const id of ['ph_download_json','ph_download_annotations']) $('#'+id).disabled=busy||!result||dirty;
    for (const node of $('#ph_csv_buttons').querySelectorAll('button')) node.disabled=busy||!result||dirty||!token;
    for (const id of ['ph_analyze','ph_save','ph_add_base','ph_set_type','ph_set_endpoint']) $('#'+id).disabled=busy||!token;
    $('#ph_use_apex').disabled=busy||!selectedAxis;
    $('#ph_set_type').disabled=busy||!token||!selectedAxis||dirty;
    $('#ph_import').disabled=busy||!measurementHash||!token;
    $('#ph_pick').disabled=busy||!graph?.points?.length||typeof project!=='function'||!scene;
  }
  function changed(note='设置已更改，请点击“应用设置并分析”；下方仍显示上次结果。') {
    dirty=true;editVersion++;
    if(activeAction==='analyze'||activeAction==='load')controller?.abort();
    card.dataset.dirty='true';message(note,'warning');updateButtons();
  }
  function setOptions(options) {
    const value={...defaults,...options};
    for(const key of Object.keys(defaults)) {
      const node=$('#ph_'+key);
      node.value=Array.isArray(value[key])?value[key].join(', '):value[key]??'';
    }
  }
  function parseList(key, label, nullable=false) {
    const text=$('#ph_'+key).value.trim();
    if(!text)return nullable?null:[];
    const values=text.split(/[,，;；\s]+/).map(Number);
    if(values.some(v=>!Number.isFinite(v)))throw Error(label+'必须填写有限数字，以逗号分隔。');
    if(values.length>512)throw Error(label+'最多允许 512 个分界值。');
    return values;
  }
  function readOptions() {
    const value={vertical_axis:$('#ph_vertical_axis').value,down_sign:Number($('#ph_down_sign').value)};
    if(!['x','y','z'].includes(value.vertical_axis)||![-1,1].includes(value.down_sign))throw Error('请选择有效的竖直轴和向下方向。');
    for(const key of ['diameter_edges','radial_edges','inclination_edges_deg','azimuth_edges_deg','bearing_edges_deg']) value[key]=parseList(key,'分箱边界');
    value.depth_edges=parseList('depth_edges','深度边界',true);
    value.origin=parseList('origin','径向参考点',true);
    if(value.origin&&value.origin.length!==3)throw Error('径向参考点需要 X、Y、Z 三个坐标，或全部留空。');
    for(const key of ['depth_bins','tangent_window','diameter_offset','diameter_window','minimum_window_fraction','minimum_radius_coverage']) {
      const raw=$('#ph_'+key).value.trim();value[key]=Number(raw);
      if(!raw||!Number.isFinite(value[key]))throw Error('层数、窗口与比例必须填写有限数字。');
    }
    if(!Number.isInteger(value.depth_bins)||value.depth_bins<1||value.depth_bins>512)throw Error('等深层数需要是 1–512 的整数。');
    for(const key of ['diameter_edges','radial_edges'])if(value[key].some(v=>v<=0))throw Error('直径和径向内部分界值必须为正数，0 与末端无上界会自动添加。');
    for(const key of ['diameter_edges','radial_edges','inclination_edges_deg','azimuth_edges_deg','bearing_edges_deg','depth_edges']) {
      const a=value[key];if(a&&a.some((v,i)=>i>0&&v<=a[i-1]))throw Error('各组分箱边界必须严格递增。');
    }
    for(const [key,max] of [['inclination_edges_deg',90],['azimuth_edges_deg',180],['bearing_edges_deg',360]]) {
      const a=value[key];if(a.length<2||a[0]!==0||a.at(-1)!==max)throw Error('角度边界必须完整覆盖 0–'+max+'°。');
    }
    if(value.depth_edges&&value.depth_edges.length<2)throw Error('指定深度边界时至少填写两个值。');
    if(value.tangent_window<=0||value.diameter_window<=0||value.diameter_offset<0)throw Error('角度和粗度窗口必须大于 0，粗度避让距离不能小于 0。');
    if(['minimum_window_fraction','minimum_radius_coverage'].some(k=>value[k]<=0||value[k]>1))throw Error('窗口与覆盖率的最小比例必须在 0–1 之间且大于 0。');
    return value;
  }
  function profile(options=readOptions()) {
    if(!measurementHash)throw Error('尚未读取当前测量的来源绑定，暂不能保存或导入标注。');
    return {schema_version:'root-phenotype-annotations-v001',measurement_sha256:measurementHash,settings:options,annotations:clone(annotations)};
  }
  async function responseError(response) {
    let value;try{value=await response.json();}catch{throw Error('请求失败（HTTP '+response.status+'）。');}
    throw Error(errors[value.error]||value.message||value.error||('请求失败（HTTP '+response.status+'）。'));
  }
  async function operate(action='analyze', kind=null) {
    if(!token)throw Error('本地交互会话不可用，请从测量工作台打开这个结果。');
    if(busy&&activeAction!=='analyze'&&activeAction!=='load')throw Error('正在保存或导出，请等待完成。');
    const options=readOptions(), annotation_profile=profile(options), version=editVersion;
    const body=JSON.stringify({action,options,annotation_profile,...(kind?{kind}:{})});
    if(new Blob([body]).size>8*1024*1024)throw Error('设置与标注请求超过 8 MiB，未发送。请减少标注文件中的冗余记录。');
    controller?.abort();controller=new AbortController();const aborter=controller;
    const id=++operationId;busy=true;activeAction=action;updateButtons();
    message(action==='save'?'正在分析并保存当前设置与标注…':action==='export'?'正在导出当前已分析结果…':'正在计算表型并核对根轴解释…');
    try {
      const response=await fetch(apiURL('/phenotyping'),{method:'POST',headers:{'Content-Type':'application/json','X-Root-Local-Token':token},
        body,signal:aborter.signal});
      if(!response.ok)await responseError(response);
      if(action==='export') {
        const blob=await response.blob();
        if(id!==operationId||version!==editVersion){message('设置已变化，本次导出未下载；请重新分析后导出。','warning');return false;}
        download(blob,'root_'+kind+'_'+measurementHash.slice(0,10)+'.csv');message('已导出当前分析的 CSV。');return true;
      }
      const value=await response.json();
      if(id!==operationId||version!==editVersion) {
        if(action==='save')message('请求开始时的版本已保存；你随后修改的设置仍未保存，请再次分析保存。','warning');
        return false;
      }
      accept(value);
      if(action==='save')message(value.saved_files?'分析与标注已保存，重新打开该结果可恢复。':'分析完成；服务器没有返回保存文件记录。',value.saved_files?'':'warning');
      else message('分析已更新。请核对候选根轴；确认状态不代表几何精度验收。');
      return true;
    } catch(error) {
      if(error.name!=='AbortError'&&id===operationId)message(errors[error.message]||error.message,'error');
      return false;
    } finally {if(id===operationId){busy=false;activeAction=null;controller=null;updateButtons();}}
  }
  async function load() {
    if(busy&&activeAction==='save')throw Error('正在保存，请完成后再恢复。');
    controller?.abort();controller=new AbortController();const aborter=controller;
    const id=++operationId,version=editVersion;busy=true;activeAction='load';updateButtons();
    message('正在读取该测量的分析与保存标注…');
    try {
      const response=await fetch(apiURL('/phenotyping.json'),{signal:aborter.signal});
      if(!response.ok)await responseError(response);const value=await response.json();
      if(id!==operationId||version!==editVersion)return;
      accept(value);message('已载入当前测量的表型分析。先指定根基，再核对候选根轴。');
    } catch(error) {if(error.name!=='AbortError'&&id===operationId)message('表型分析未载入：'+error.message,'error');}
    finally {if(id===operationId){busy=false;activeAction=null;controller=null;updateButtons();}}
  }
  function accept(value) {
    if(!value||typeof value!=='object'||!value.distributions||!value.architecture)throw Error('服务器返回的表型分析结构不完整。');
    const hash=value.measurement?.sha256||value.annotation_profile?.measurement_sha256;
    if(typeof hash!=='string'||!/^[a-f0-9]{64}$/i.test(hash))throw Error('返回结果没有有效的测量来源绑定。');
    if(measurementHash&&hash!==measurementHash)throw Error('测量来源已改变，请重新打开结果后再标注。');
    result=value;measurementHash=hash;annotations={...emptyAnnotations(),...clone(value.annotation_profile?.annotations||value.architecture.annotations||{})};
    dirty=false;card.dataset.dirty='false';card.dataset.ready='true';setOptions(value.settings||value.annotation_profile?.settings||{});
    render();updateButtons();
  }
  function table(parent, columns, records, limit=180) {
    const wrap=add(parent,'div',null,'ph-table'),t=add(wrap,'table'),head=add(t,'thead'),hr=add(head,'tr');
    for(const [label] of columns)add(hr,'th',label);
    const body=add(t,'tbody');
    for(const row of records.slice(0,limit)) {
      const tr=add(body,'tr');for(const [,get] of columns){const td=add(tr,'td');const value=get(row);if(value instanceof Node)td.append(value);else td.textContent=String(value??'不可用');}
    }
    if(!records.length){const tr=add(body,'tr');const td=add(tr,'td','暂无可显示记录');td.colSpan=columns.length;}
    if(records.length>limit)add(parent,'p','界面显示前 '+limit+' / '+records.length+' 条；完整记录可导出 CSV。','ph-note');
  }
  function interval(row,prefix,angle=false) {
    const a=row[prefix+'_start'],b=row[prefix+'_end'];
    return format(a)+'–'+(b==null?'∞':format(b))+(angle?'°':'');
  }
  function bars(parent,hist,title,prefix,angle=false) {
    parent.replaceChildren();add(parent,'h3',title);
    if(!hist){add(parent,'p','当前没有该分布。','ph-note');return;}
    const note=add(parent,'p',statusLabel(hist.status)+' · 已分配弧长 '+format(hist.assigned_axis_length)+' '+unit()+' · 占全部已提取轴 '+percent(hist.covered_arc_fraction),'ph-note');note.title=hist.status;
    if(finite(hist.unassigned_axis_length)&&hist.unassigned_axis_length>Math.abs(hist.total_axis_length||0)*1e-12)add(parent,'p','未分配弧长 '+format(hist.unassigned_axis_length)+' '+unit()+'，保留在总分母中。','ph-note');
    const rows=hist.bins||[],container=add(parent,'div',null,'ph-bars');
    let maximum=0;for(const row of rows)if(finite(row.axis_length))maximum=Math.max(maximum,row.axis_length);
    for(const row of rows.slice(0,64)) {
      const line=add(container,'div',null,'ph-bar');add(line,'span',interval(row,prefix,angle));
      const track=add(line,'div',null,'ph-track'),fill=add(track,'div',null,'ph-fill');
      fill.style.width=(maximum&&finite(row.axis_length)?100*row.axis_length/maximum:0)+'%';
      add(line,'span',format(row.axis_length));line.title='占总轴长 '+percent(row.fraction_of_total_axis)+'；占已分配轴长 '+percent(row.fraction_of_assigned_axis);
    }
    if(rows.length>64)add(parent,'p','图中显示前 64 箱；完整分箱见 CSV。','ph-note');
  }
  function render() {
    const d=result.distributions,a=result.architecture,t=result.axis_traits;
    $('#ph_units').textContent='所有长度、窗口、直径分界和参考点坐标均使用 '+unit()+'。'+(report.scale_mm_per_model_unit==null?'当前模型尚未给定毫米尺度。':'毫米尺度来自原测量报告。');
    $('#ph_binding').textContent='测量 SHA-256：'+measurementHash;
    const stats=$('#ph_stats');stats.replaceChildren();
    const interpretable=(a.axes||[]).filter(axis=>['CANDIDATE','CONFIRMED','AUTO_INFERRED'].includes(axis.confirmation_status)).length;
    for(const [label,value,note] of [['已提取轴长',format(d.total_axis_length)+' '+unit(),d.partial_axis?'仅当前提取部分':'完整性与精度另验'],
      ['合格半径覆盖',percent(d.sample_coverage?.qualified_radius_arc_fraction),'以已提取轴长为分母'],
      [a.analysis_regions?'可分析有向根段':'可解释根轴',interpretable?String(interpretable)+' 条':'待指定',statusLabel(a.status)+'；'+(a.axes?.length??0)+' 个几何分段']]) {
      const cell=add(stats,'div',null,'ph-stat');add(cell,'small',label);add(cell,'strong',value);add(cell,'small',note);
    }
    renderBases();renderComponents();renderJunctions();renderAxisSelect();
    const note=$('#ph_distribution_note');note.replaceChildren();
    add(note,'p','直径统计按合格单分支半径样本的弧长权重计算。深度零点：'+(d.depth_reference?.source==='SOURCE_MATERIAL_BOUNDS'?'原实体上边界':'已提取轴上边界')+'。');
    const h=d.distributions||{};bars($('#ph_diameter_chart'),h.diameter,'直径分级根长','diameter');
    const directions=$('#ph_direction_charts');directions.replaceChildren();
    bars(add(directions,'div'),h.inclination,'与竖直方向的无向夹角','inclination_deg',true);
    bars(add(directions,'div'),h.azimuth,'水平无向根段朝向（0–180°）','azimuth_deg',true);
    const horizontal=(d.horizontal_axes||['x','y']).map(x=>x.toUpperCase()).join(' → ');
    add(directions,'p','朝向在 '+horizontal+' 平面计算，反向视为同一轴向。这描述根段切线方向，不代表根位于哪个空间象限。竖直根段的水平朝向未定义，不放入 0° 箱。','ph-note');
    bars($('#ph_bearing_chart'),h.bearing,'相对参考点的位置方位（0–360°）','bearing_deg',true);
    add($('#ph_bearing_chart'),'p','位置方位从参考点指向根段所在位置；正角按 '+horizontal+' 增加。必须显式指定参考点，参考竖直轴上的根段保留为未分配。','ph-note');
    bars($('#ph_radial_chart'),h.radial,'距参考竖直轴的径向根长','radial_distance');
    if(d.settings?.origin)add($('#ph_radial_chart'),'p','参考点：'+d.settings.origin.map(format).join(', ')+' '+unit()+'。','ph-note');
    const joint=$('#ph_joint_tables');joint.replaceChildren();
    add(joint,'h3','直径 × 深度');add(joint,'p','样本位置所在深度箱承接其弧长权重；这是采样支持估计，未插补缺失半径。','ph-note');
    if(h.diameter_depth){add(joint,'p',statusLabel(h.diameter_depth.status)+' · 覆盖 '+percent(h.diameter_depth.covered_arc_fraction),'ph-note');table(joint,[['直径 / '+unit(),r=>interval(r,'diameter')],['深度 / '+unit(),r=>interval(r,'depth')],['根长 / '+unit(),r=>format(r.axis_length)],['占总轴长',r=>percent(r.fraction_of_total_axis)]],h.diameter_depth.bins||[]);}
    add(joint,'h3','深度 × 无向根段朝向（0–180°）');add(joint,'p','按深度边界切分线段后计长；表示各层根段的切线方向，不是各层的空间象限。','ph-note');
    if(h.depth_azimuth){add(joint,'p',statusLabel(h.depth_azimuth.status)+' · 覆盖 '+percent(h.depth_azimuth.covered_arc_fraction),'ph-note');table(joint,[['深度 / '+unit(),r=>interval(r,'depth')],['朝向',r=>interval(r,'azimuth_deg',true)],['根长 / '+unit(),r=>format(r.axis_length)],['占总轴长',r=>percent(r.fraction_of_total_axis)]],h.depth_azimuth.bins||[]);}
    add(joint,'h3','深度 × 相对位置方位（0–360°）');add(joint,'p','相对显式参考点按深度与位置方位共同分箱，可比较各层不同空间扇区的根长。','ph-note');
    if(h.depth_bearing){add(joint,'p',statusLabel(h.depth_bearing.status)+' · 覆盖 '+percent(h.depth_bearing.covered_arc_fraction),'ph-note');table(joint,[['深度 / '+unit(),r=>interval(r,'depth')],['位置方位',r=>interval(r,'bearing_deg',true)],['根长 / '+unit(),r=>format(r.axis_length)],['占总轴长',r=>percent(r.fraction_of_total_axis)]],h.depth_bearing.bins||[]);}
    renderAxisDetails();renderRelations();renderSummaries(t?.summaries||{});
  }
  function topologyChanged(clearContinuations=false) {
    annotations.confirmed_axes=[];annotations.confirmed_components=[];annotations.root_types=[];
    if(clearContinuations)annotations.junction_continuations=[];
    selectedAxis=null;
  }
  async function annotationEdit(edit,{topology=false,base=false,confirmation=false}={}) {
    if(busy&&activeAction==='save')throw Error('正在保存，请完成后再修改标注。');
    edit();if(topology)topologyChanged(base);
    else if(confirmation){annotations.confirmed_axes=[];annotations.confirmed_components=[];}
    changed(topology?'根轴划分已更改，旧确认和根类型标签已清除，正在重新分析…':confirmation?'端点类型已更改，旧确认已清除，正在重新分析…':'标注已更改，正在重新分析…');
    renderBases();await operate('analyze');
  }
  function nodeInput(id) {
    const text=$(id).value.trim(),node=Number(text);
    if(!text||!Number.isInteger(node)||node<0||node>=(graph?.points?.length||0))throw Error('请填写当前中心线中有效的整数节点编号。');
    return node;
  }
  function renderBases() {
    const box=$('#ph_bases');box.replaceChildren();
    if(!annotations.root_bases.length)add(box,'p','尚未指定根基。直径与无向方向分布仍可使用；根序和父子关系暂不可用。','ph-note');
    for(const base of annotations.root_bases) {
      const row=add(box,'div',null,'ph-chip');add(row,'span','节点 '+base.node_id+' · '+(base.root_type||'unknown'));
      button(row,'设原点',()=>setOriginFromBase(base.node_id));
      button(row,'移除',()=>annotationEdit(()=>{
        annotations.root_bases=annotations.root_bases.filter(b=>b.node_id!==base.node_id);
        annotations.endpoint_types=annotations.endpoint_types.filter(e=>e.node_id!==base.node_id||e.type!=='root_base');
      },{topology:true,base:true}),false,true);
    }
  }
  function renderComponents() {
    const box=$('#ph_components');box.replaceChildren();
    for(const component of result.architecture.components||[]) {
      const row=add(box,'div',null,'ph-junction');add(row,'strong','分量 '+component.component_id+' · '+statusLabel(component.status));
      add(row,'p',component.node_count+' 个节点，'+component.edge_count+' 条边；'+(component.axis_ids?.length||0)+' 个当前几何分段。','ph-note');
      if(component.status==='PARTIAL_ROOTED_CYCLIC_COMPONENT')add(row,'p','环外有向根段覆盖 '+percent(component.rooted_fragment_arc_fraction)+'；环内 '+component.cyclic_core_edge_count+' 条原边完整保留。未选择穿环路径；跨环全局根序及完整母根长度不可用。','ph-note');
      const confirmed=component.status==='ROOTED_CONFIRMED';
      const eligible=['ROOTED_CANDIDATE','ROOTED_PARTIALLY_CONFIRMED','ROOTED_CONFIRMED','ROOTED_AUTO_INFERRED','ROOTED_REVIEW_RECOMMENDED'].includes(component.status);
      button(row,confirmed?'撤回此分量确认':'已核对：确认此分量',()=>annotationEdit(()=>{
        const id=component.component_id;
        if(confirmed){annotations.confirmed_components=annotations.confirmed_components.filter(x=>x!==id);annotations.confirmed_axes=annotations.confirmed_axes.filter(x=>!component.axis_ids.includes(x));}
        else annotations.confirmed_components=[...new Set([...annotations.confirmed_components,id])];
      }),!eligible||dirty,true);
      if(!eligible)add(row,'p',component.status==='ROOTED_UNRESOLVED_CONTINUATION'?'请在“核对分叉续接”中处理歧义后再确认。':'该分量尚不具备确认根轴层级的条件。','ph-note');
    }
  }
  function renderJunctions() {
    const box=$('#ph_junctions');box.replaceChildren();const decisions=result.architecture.junction_decisions||[];
    if(!decisions.length)add(box,'p','指定根基后，有多个后续分支的节点会列在这里。','ph-note');
    for(const decision of decisions.slice(0,150)) {
      const row=add(box,'div',null,'ph-junction');add(row,'strong','分叉节点 '+decision.node_id+' · '+statusLabel(decision.status));
      add(row,'p','入边来自节点 '+decision.incoming_node_id+'；'+statusLabel(decision.confidence_tier)+'。自动结果由局部方向、合格半径连续性、下游长度和层级综合推断，不代表生物学验证。','ph-note');
      const select=add(row,'select');select.setAttribute('aria-label','节点 '+decision.node_id+' 的续接');
      for(const [value,label] of [['auto','采用自动推断'],['stop','本根在此终止；子分支另起根轴'],...(decision.alternatives||[]).map(x=>[String(x.outgoing_node_id),'继续到节点 '+x.outgoing_node_id+'（综合分数 '+format(x.ranking_score)+'，方向 '+format(x.alignment_cosine)+'，径比 '+format(x.radius_ratio_to_mother)+'，下游长 '+format(x.downstream_max_path_length)+'）'])]){const option=add(select,'option',label);option.value=value;}
      const override=annotations.junction_continuations.find(x=>x.node_id===decision.node_id);
      select.value=override?(override.outgoing_node_id===null?'stop':String(override.outgoing_node_id)):'auto';
      const actions=add(row,'div',null,'ph-actions');
      button(actions,'高亮分叉',()=>{const p=graph.points[decision.node_id];highlight?.((decision.alternatives||[]).flatMap(x=>[graph.points[decision.incoming_node_id],p,graph.points[x.outgoing_node_id],null]));});
      button(actions,'应用续接',()=>annotationEdit(()=>{
        annotations.junction_continuations=annotations.junction_continuations.filter(x=>x.node_id!==decision.node_id);
        if(select.value!=='auto')annotations.junction_continuations.push({node_id:decision.node_id,incoming_node_id:decision.incoming_node_id,outgoing_node_id:select.value==='stop'?null:Number(select.value)});
      },{topology:true}),dirty,true);
    }
    if(decisions.length>150)add(box,'p','界面显示前 150 个分叉；其余可通过导出、编辑并重新导入标注 JSON 核对。','ph-note');
  }
  function renderAxisSelect() {
    const select=$('#ph_axis_select');select.replaceChildren();const axes=result.architecture.axes||[];
    if(!axes.some(a=>a.axis_id===selectedAxis))selectedAxis=axes[0]?.axis_id||null;
    for(const [index,axis] of axes.entries()) {
      const option=add(select,'option',(index+1)+' · 分量 '+axis.component_id+' · '+statusLabel(axis.confirmation_status)+' · '+format(axis.length_model_units*(report.scale_mm_per_model_unit??1))+' '+unit());
      option.value=axis.axis_id;
    }
    if(!axes.length){const option=add(select,'option','暂无根轴');option.value='';}
    select.value=selectedAxis||'';select.disabled=!axes.length;
    $('#ph_axis_type').value=axes.find(a=>a.axis_id===selectedAxis)?.root_type||'unknown';
  }
  function axisRecord() {return result?.architecture?.axes?.find(a=>a.axis_id===selectedAxis);}
  function showAxis() {
    const axis=axisRecord();if(!axis)return;
    highlight?.(axis.points||axis.node_ids.map(i=>graph.points[i]));
    $('#ph_axis_type').value=axis.root_type||'unknown';renderAxisDetails();renderRelations();updateButtons();
  }
  function propertyRows(parent,rows) {table(parent,[['指标',r=>r[0]],['数值',r=>r[1]],['依据',r=>r[2]||'']],rows);}
  function renderAxisDetails() {
    const box=$('#ph_axis_details');box.replaceChildren();const geometry=axisRecord();
    if(!geometry){add(box,'p','请选择一个根轴。','ph-note');return;}
    const axis=result.axis_traits?.axes?.find(x=>x.axis_id===selectedAxis);
    const short=selectedAxis.slice(-8);add(box,'h3','根轴 …'+short+' · '+statusLabel(geometry.confirmation_status));
    add(box,'p','根类型：'+geometry.root_type+'；发育根序：'+(geometry.developmental_order!=null?geometry.developmental_order:geometry.candidate_developmental_order!=null?geometry.candidate_developmental_order+'（候选）':'不可用')+'。基轴从 0 计，侧根为母根级次加 1。','ph-note');
    if(!axis){add(box,'p','该根轴的详细性状尚不可用。','ph-note');return;}
    if(axis.axis_scope==='ROOTED_BRIDGE_FRAGMENT')add(box,'p','当前显示环外局部有向根段。夹角、粗度和渐细只在本段内计算；穿环续接与全局根序未推断，分段数不能当作生物根数。','ph-note');
    if(axis.axis_scope==='ROOTED_BRIDGE_FRAGMENT'&&axis.local_candidate_developmental_order!=null)add(box,'p','本片段内的候选级次：'+axis.local_candidate_developmental_order+'（从片段入口计0，不是整株发育根序）。','ph-note');
    if(axis.axis_scope==='CYCLIC_CORE')add(box,'p','本段属于保留的环内几何。原边和长度完整计入总量；穿环方向、父子身份和完整根序尚未确定。','ph-note');
    const branching=axis.branching||{},diameter=axis.diameter_distribution||{},taper=axis.taper||{};
    const metric=(name,value,u='')=>[name,format(value?.value)+(value?.value==null?'':' '+u),statusLabel(value?.status)];
    propertyRows(box,[['根轴长度',format(axis.length)+' '+unit(),statusLabel(axis.confirmation_status)],
      ['平均直径',format(diameter.mean)+' '+unit(),'有效弧长 '+percent(diameter.valid_arc_fraction)],
      ['直径中位数 / P95',format(diameter.median)+' / '+format(diameter.p95)+' '+unit(),'仅合格单分支样本'],
      metric('侧根数量',branching.lateral_count,'条'),metric('全母根侧根密度',branching.lrbd_whole_axis,unit()+'⁻¹'),
      metric('分枝区侧根密度',branching.lrbd_branched_zone,unit()+'⁻¹'),
      ...(branching.lrbd_observed_fragment?[metric('局部观测根段侧分枝密度',branching.lrbd_observed_fragment,unit()+'⁻¹')]:[]),
      ['侧根间距中位数',format(branching.interbranch_distance_statistics?.median)+' '+unit(),statusLabel(branching.interbranch_distance_statistics?.status)],
      ['侧根间距 CV',format(branching.interbranch_distance_statistics?.cv),'沿母根弧长，不是空间直线距离'],
      metric('基部未分枝区长度',branching.basal_unbranched_length,unit()),
      metric('根尖未分枝区长度',branching.apical_unbranched_length,unit()),
      metric('末侧根至观测末端',branching.observed_apical_unbranched_length,unit()),
      ['渐细率（直径 / 弧长）',format(taper.taper_rate),statusLabel(taper.status)+'；正值表示向末端变细'],['渐细回归 R²',format(taper.r_squared),'描述拟合，不代表测量精度'],
      [axis.axis_scope==='ROOTED_BRIDGE_FRAGMENT'?'当前根段对向下方向夹角':'整根对向下方向夹角',format(axis.whole_axis_direction?.growth_angle_from_down_deg)+'°',statusLabel(axis.whole_axis_direction?.status)]]);
    add(box,'p','全母根密度 = 侧根数 / 母根长度；分枝区密度 = 侧根数 / 首末侧根之间弧长。真实根尖未确认时，根尖未分枝区保持不可用。','ph-note');
    drawProfile(box,axis);
    const details=add(box,'details');add(details,'summary','沿根轴直径样本与缺失');
    table(details,[['轴向距离 / '+unit(),r=>format(r.s)],['直径 / '+unit(),r=>format(r.diameter)],['支持弧长 / '+unit(),r=>format(r.arc_weight)],['状态',r=>r.available?'有效单分支样本':'缺失／不合格']],axis.diameter_profile||[]);
  }
  function drawProfile(parent,axis) {
    const samples=axis.diameter_profile||[],canvas=add(parent,'canvas');canvas.width=640;canvas.height=360;
    canvas.setAttribute('role','img');canvas.setAttribute('aria-label','所选根轴沿弧长的合格直径样本；红线为缺失弧段');
    const ctx=canvas.getContext('2d');if(!ctx)return;ctx.scale(2,2);const width=320,height=180,left=44,right=12,top=21,bottom=30;
    let ymax=0;for(const sample of samples)if(finite(sample.diameter))ymax=Math.max(ymax,sample.diameter);
    const xmax=axis.length||1;ymax=ymax?ymax*1.1:1;
    const x=s=>left+(width-left-right)*s/xmax,y=d=>height-bottom-(height-top-bottom)*d/ymax;
    ctx.strokeStyle='#657f8f';ctx.beginPath();ctx.moveTo(left,top);ctx.lineTo(left,height-bottom);ctx.lineTo(width-right,height-bottom);ctx.stroke();
    ctx.fillStyle='#b7ccd6';ctx.font='10px system-ui';ctx.fillText('直径 / '+unit(),6,12);ctx.fillText(format(ymax),3,top+5);ctx.fillText('0',left-10,height-bottom+13);ctx.fillText(format(xmax)+' '+unit(),width-85,height-7);
    const step=Math.max(1,Math.ceil(samples.length/1200));
    for(let i=0;i<samples.length;i++) {
      const row=samples[i];if(!row.available){ctx.strokeStyle='#e68c82';ctx.beginPath();ctx.moveTo(x(row.support_start),height-bottom+3);ctx.lineTo(x(row.support_end),height-bottom+3);ctx.stroke();}
      else if(i%step===0&&finite(row.diameter)){ctx.fillStyle='#83ddc3';ctx.beginPath();ctx.arc(x(row.s),y(row.diameter),2,0,Math.PI*2);ctx.fill();}
    }
    add(parent,'p','横轴为沿线弧长，'+(axis.arc_origin==='axis_base'?'从根轴基部计起；':axis.arc_origin==='root_relative_bridge_fragment_start'?'从当前有向片段入口计起；':'从几何路径的任意端点计起，未定义生长方向；')+'绿点为合格直径样本，红线为缺失支持区间。'+(step>1?'显示已抽稀，导出保留全部样本。':''),'ph-note');
  }
  function renderRelations() {
    const box=$('#ph_relations');box.replaceChildren();const all=result?.axis_traits?.branch_relations||[];
    const rows=selectedAxis?all.filter(r=>r.parent_axis_id===selectedAxis||r.child_axis_id===selectedAxis):all;
    if(result.architecture.analysis_regions)add(box,'p','含环结果只列同一桥片段内的候选关系与当前片段PCA角；不推断穿环母侧身份。','ph-note');
    add(box,'p','显示与所选根轴有关的 '+rows.length+' 条关系。局部切线角与整根 PCA 角分别列出。','ph-note');
    table(box,[['母根 → 子根',r=>r.parent_axis_id.slice(-6)+' → '+r.child_axis_id.slice(-6)],
      ['附着弧长 / '+unit(),r=>format(r.attachment_arc)],['局部夹角',r=>format(r.local_insertion_angle_deg?.value)+'°'],
      ['整根夹角',r=>format(r.whole_axis_pca_angle_deg?.value)+'°'],['母根直径 / '+unit(),r=>format(r.parent_diameter?.value)],
      ['子根直径 / '+unit(),r=>format(r.child_diameter?.value)],['子/母直径',r=>format(r.child_parent_diameter_ratio?.value)],
      ['状态',r=>statusLabel(r.confirmation_status)]],rows);
    if(rows.length) {
      const details=add(box,'details');add(details,'summary','实际角窗、粗度窗与有效覆盖');
      add(details,'p','下表窗口单位为 '+unit()+'，格式为“实际 / 设定”。窗口靠近端点时可缩短；请结合覆盖率比较父子直径。','ph-note');
      const windowText=value=>format(value?.effective_window)+' / '+format(value?.requested_window)+(value?.window_shortened?'（缩短）':'');
      table(details,[['子根',r=>r.child_axis_id.slice(-6)],['母根角窗',r=>windowText(r.parent_tangent)],['子根角窗',r=>windowText(r.child_tangent)],
        ['局部角状态',r=>statusLabel(r.local_insertion_angle_deg?.status)],['母根粗度窗',r=>windowText(r.parent_diameter)],['子根粗度窗',r=>windowText(r.child_diameter)],
        ['母根半径覆盖',r=>percent(r.parent_diameter?.valid_arc_fraction)],['子根半径覆盖',r=>percent(r.child_diameter?.valid_arc_fraction)],
        ['直径比状态',r=>statusLabel(r.child_parent_diameter_ratio?.status)]],rows);
    }
    if(rows.some(r=>r.local_insertion_angle_deg?.value==null||r.child_parent_diameter_ratio?.value==null))add(box,'p','部分角度或粗度缺失：可能没有足够窗口、根轴方向未确定或有效半径覆盖不足；请查看分支 CSV 中的逐项状态。','ph-note');
  }
  function renderSummaries(summaries) {
    const box=$('#ph_summaries');box.replaceChildren();
    for(const [key,title,label] of [['by_developmental_order','发育根序','order'],['by_root_type','人工根类型','root_type'],['by_confirmation_status','确认状态','confirmation_status'],['by_strahler_edge_order','Strahler 边级次','order']]) {
      add(box,'h3',title);table(box,[['分组',r=>label==='confirmation_status'?statusLabel(r[label]):r[label]],
        ['状态',r=>statusLabel(r.order_status||r.confirmation_status)],['根轴数',r=>r.axis_count],
        ['几何分段 / 边段数',r=>r.geometric_partition_count??r.geometric_segment_count],
        ['总长度 / '+unit(),r=>format(r.total_length)],['平均直径 / '+unit(),r=>format(r.diameter_distribution?.mean)],
        ['半径覆盖',r=>percent(r.diameter_distribution?.valid_arc_fraction)]],summaries[key]||[]);
    }
    add(box,'p','候选与确认后的统计分开列出。Strahler 是图边级次，同一条根轴可包含不同级次；未知根类型保留 unknown。','ph-note');
  }
  function download(blob,name) {
    const url=URL.createObjectURL(blob),link=document.createElement('a');link.href=url;link.download=name;document.body.append(link);link.click();link.remove();setTimeout(()=>URL.revokeObjectURL(url),1000);
  }
  function jsonDownload(value,name) {download(new Blob([JSON.stringify(value,null,2)+'\n'],{type:'application/json;charset=utf-8'}),name);}
  function setOriginFromBase(node) {
    if(!annotations.root_bases.some(base=>base.node_id===node))throw Error('请先将这个节点添加为根基，或在参考点框手动填写坐标。');
    const factor=report.scale_mm_per_model_unit??1;
    $('#ph_origin').value=graph.points[node].map(v=>v*factor).join(', ');
    $('#ph_settings').open=true;
    changed('已将根基节点 '+node+' 的坐标明确设为位置分布参考点（'+unit()+'）。请应用设置重新分析。');
  }
  function validateImportedProfile(value) {
    if(!value||typeof value!=='object'||Array.isArray(value)||Object.keys(value).some(k=>!['schema_version','measurement_sha256','settings','annotations'].includes(k)))throw Error('标注文件顶层字段不受支持。');
    const settings=value.settings===undefined?{}:value.settings;
    if(!settings||typeof settings!=='object'||Array.isArray(settings)||Object.keys(settings).some(k=>!Object.prototype.hasOwnProperty.call(defaults,k)))throw Error('标注中存在当前版本不支持的设置，请使用匹配版本。');
    const lists=['diameter_edges','radial_edges','inclination_edges_deg','azimuth_edges_deg','bearing_edges_deg','depth_edges','origin'];
    for(const [key,item] of Object.entries(settings)) {
      if(lists.includes(key)) {
        if(item===null&&['depth_edges','origin'].includes(key))continue;
        if(!Array.isArray(item)||item.some(x=>!finite(x)))throw Error('标注中的 '+key+' 必须是有限数字数组，不能含布尔值或空值。');
        if(key==='origin'&&item.length!==3||key==='depth_edges'&&item.length<2)throw Error('标注中的参考点或深度边界数量不正确。');
      } else if(key==='vertical_axis') {
        if(!['x','y','z'].includes(item))throw Error('标注中的竖直轴必须是 x、y 或 z。');
      } else if(!finite(item))throw Error('标注中的 '+key+' 必须是有限数字，不能用字符串或布尔值代替。');
    }
    const a=value.annotations,keys=Object.keys(emptyAnnotations());
    if(!a||typeof a!=='object'||Array.isArray(a)||Object.keys(a).some(k=>!keys.includes(k)))throw Error('标注结构中有不受支持的字段。');
    if(a.schema_version!=='root-axis-annotations-v001')throw Error('根轴标注版本不匹配。');
    for(const key of keys.filter(k=>k!=='schema_version'))if(key in a&&!Array.isArray(a[key]))throw Error('标注中的 '+key+' 必须是列表。');
    if(a.confirmed_axes?.some(x=>typeof x!=='string'||!x.trim()))throw Error('根轴确认列表必须包含非空标识符。');
    if(a.confirmed_components?.some(x=>!Number.isInteger(x)||x<0))throw Error('分量确认列表必须包含非负整数编号。');
    for(const key of ['root_bases','junction_continuations','endpoint_types','root_types'])if(a[key]?.some(x=>!x||typeof x!=='object'||Array.isArray(x)))throw Error('标注列表中的每项必须是记录对象。');
  }
  async function importFile(file) {
    if(!file)return;if(file.size>8*1024*1024)throw Error('标注文件超过 8 MB，请检查是否选择了正确的标注 JSON。');
    const raw=JSON.parse(await file.text()),value=raw?.annotation_profile||raw;
    validateImportedProfile(value);
    if(value.schema_version!=='root-phenotype-annotations-v001')throw Error('请选择本测量器导出的表型标注 JSON。');
    if(!measurementHash||value.measurement_sha256!==measurementHash)throw Error('标注绑定的测量 SHA-256 与当前文件不同，未导入。');
    const previousOptions=readOptions(),previousAnnotations=annotations,previousDirty=dirty;
    try{setOptions(value.settings||previousOptions);readOptions();annotations={...emptyAnnotations(),...clone(value.annotations)};}
    catch(error){setOptions(previousOptions);annotations=previousAnnotations;throw error;}
    changed('来源绑定匹配，正在检查导入标注…');const version=editVersion;let accepted=false;
    try{accepted=await operate('analyze');}catch(error){message(errors[error.message]||error.message,'error');}
    if(!accepted&&version===editVersion){setOptions(previousOptions);annotations=previousAnnotations;dirty=previousDirty;card.dataset.dirty=String(dirty);renderBases();updateButtons();$('#ph_status').textContent+=' 导入未通过，已恢复导入前的设置与标注。';}
  }
  function pickNode(event) {
    const rect=scene.getBoundingClientRect(),x=event.clientX-rect.left,y=event.clientY-rect.top;let best=null;
    for(let i=0;i<(graph?.points?.length||0);i++) {
      const point=project(graph.points[i]);if(!point||!point.every(Number.isFinite))continue;
      const distance=Math.hypot(point[0]-x,point[1]-y);
      if(distance<=16&&(!best||distance<best.distance-.5||Math.abs(distance-best.distance)<.5&&point[2]<best.depth))best={node:i,distance,depth:point[2]};
    }
    if(!best){$('#ph_pick_note').textContent='附近没有可选轴节点。请放大、旋转模型后重试，也可手动输入编号。';return;}
    $('#ph_base_node').value=String(best.node);$('#ph_endpoint_node').value=String(best.node);
    const point=graph.points[best.node],factor=report.scale_mm_per_model_unit??1;
    $('#ph_pick_note').textContent='已选节点 '+best.node+'；坐标 '+point.map(v=>format(v*factor)).join(', ')+' '+unit()+'；连接 '+graph.adjacency[best.node].length+' 条边。点击“添加根基”才会更改标注。';
    const bounds=report.model_bounds;let extent=1;if(bounds)extent=Math.max(...bounds.maximum.map((v,i)=>v-bounds.minimum[i]),1e-8);
    const cross=[];for(let k=0;k<3;k++){const a=point.slice(),b=point.slice();a[k]-=extent*.008;b[k]+=extent*.008;cross.push(a,b,null);}highlight?.(cross);
  }
  if(scene) {
    // Capture selection clicks without replacing the existing rotate/measure handlers.
    scene.addEventListener('pointerdown',event=>{if(!pickMode||event.button!==0)return;event.preventDefault();event.stopImmediatePropagation();pickDown={x:event.clientX,y:event.clientY,id:event.pointerId};scene.setPointerCapture?.(event.pointerId);},true);
    scene.addEventListener('pointermove',event=>{if(!pickMode||!pickDown||event.pointerId!==pickDown.id)return;event.preventDefault();event.stopImmediatePropagation();},true);
    scene.addEventListener('pointerup',event=>{if(!pickMode||!pickDown||event.pointerId!==pickDown.id)return;event.preventDefault();event.stopImmediatePropagation();const start=pickDown;pickDown=null;scene.releasePointerCapture?.(event.pointerId);if(Math.hypot(event.clientX-start.x,event.clientY-start.y)<8)guard(()=>pickNode(event));},true);
    scene.addEventListener('pointercancel',()=>{pickDown=null;},true);
  }
  setOptions(defaults);
  for(const key of Object.keys(defaults))$('#ph_'+key).addEventListener('input',()=>changed());
  $('#ph_analyze').onclick=()=>guard(()=>operate('analyze'));
  $('#ph_save').onclick=()=>guard(()=>operate('save'));
  $('#ph_reload').onclick=()=>guard(()=>load());
  $('#ph_origin_from_base').onclick=()=>guard(()=>setOriginFromBase(nodeInput('#ph_base_node')));
  $('#ph_pick').onclick=()=>{pickMode=!pickMode;pickDown=null;$('#ph_pick').textContent=pickMode?'退出选点模式':'在模型上选节点';$('#ph_pick').setAttribute('aria-pressed',String(pickMode));$('#ph_pick_note').textContent=pickMode?'选点模式：左键选择轴节点，右键拖动旋转，滚轮缩放。选中后点击“添加根基”，或在端点类型处标注。':'已退出选点模式，模型恢复原有操作方式。';};
  $('#ph_add_base').onclick=()=>guard(()=>{const node=nodeInput('#ph_base_node'),type=$('#ph_base_type').value.trim()||'unknown';return annotationEdit(()=>{
    annotations.root_bases=annotations.root_bases.filter(b=>b.node_id!==node);annotations.root_bases.push({node_id:node,root_type:type});
    annotations.endpoint_types=annotations.endpoint_types.filter(e=>e.node_id!==node);annotations.endpoint_types.push({node_id:node,type:'root_base'});
  },{topology:true,base:true});});
  $('#ph_axis_select').onchange=()=>{selectedAxis=$('#ph_axis_select').value||null;showAxis();};
  $('#ph_set_type').onclick=()=>guard(()=>{if(!selectedAxis)throw Error('请先选择根轴。');const id=selectedAxis,type=$('#ph_axis_type').value.trim()||'unknown';return annotationEdit(()=>{annotations.root_types=annotations.root_types.filter(r=>r.axis_id!==id);annotations.root_types.push({axis_id:id,type});});});
  $('#ph_use_apex').onclick=()=>guard(()=>{const axis=axisRecord();if(!axis)throw Error('请先选择根轴。');const node=axis.tip_node_id??axis.node_ids.at(-1);$('#ph_endpoint_node').value=String(node);$('#ph_endpoint_type').value=annotations.endpoint_types.find(e=>e.node_id===node)?.type||'unknown';});
  $('#ph_set_endpoint').onclick=()=>guard(()=>{const node=nodeInput('#ph_endpoint_node'),type=$('#ph_endpoint_type').value;const base=annotations.root_bases.some(b=>b.node_id===node);
    if(type==='root_base'&&!base)throw Error('请先将此节点添加为根基。');
    if(['root_tip','cut_end'].includes(type)&&(base||graph.adjacency[node].length!==1))throw Error('根尖和截断端必须是非根基的图端点。');
    return annotationEdit(()=>{annotations.endpoint_types=annotations.endpoint_types.filter(e=>e.node_id!==node);annotations.endpoint_types.push({node_id:node,type});},{confirmation:true});});
  $('#ph_clear_highlight').onclick=()=>highlight?.([]);
  $('#ph_download_json').onclick=()=>guard(()=>{if(!result||dirty)throw Error('请先更新分析。');jsonDownload(result,'root_phenotyping_'+measurementHash.slice(0,10)+'.json');});
  $('#ph_download_annotations').onclick=()=>guard(()=>{if(!result||dirty)throw Error('请先更新分析。');jsonDownload(profile(),'root_annotations_'+measurementHash.slice(0,10)+'.json');});
  $('#ph_import').onclick=()=>$('#ph_import_file').click();
  $('#ph_import_file').onchange=()=>{const file=$('#ph_import_file').files[0];$('#ph_import_file').value='';guard(()=>importFile(file));};
  for(const [kind,label] of [['distributions','分布 CSV'],['axes','根轴 CSV'],['branches','父子分支 CSV'],['profiles','沿线粗度 CSV'],['summaries','分组汇总 CSV']])button($('#ph_csv_buttons'),label,()=>operate('export',kind));
  updateButtons();
  try {const response=await fetch('/api/interaction-session');if(response.ok)token=(await response.json()).token||null;}catch{/* Read-only analysis can still be shown. */}
  await load();
  if(!token)message('分析可查看；本地交互会话不可用，标注、分析更新与保存暂不可用。请从测量工作台打开该结果。','warning');
  return {card,reload:load,getAnalysis:()=>result};
}
