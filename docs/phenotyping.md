# 根系表型分析与根轴标注

`phenotype` 在已生成的 `*_measurement.json` 上增加直径及空间分布、根轴身份标注、母侧根关系和沿轴粗度分析。输入图的节点、边、原实体边界和原截面半径保持不变；每次保存生成新的分析快照。

本模块区分三种信息：**观测几何、用户确认的根轴解释、独立的生物学准确率验证**。确认一条根轴，表示用户接受它的连续路径与身份；它不会自动把漏检、误连接或半径误差变成已经通过验证的结果。

## 1. 使用方式

### 命令行

先对已有测量文件计算无需根轴身份的分布：

```powershell
root-measurer phenotype "runs/default/results/model_delivery_local/run001/model_measurement.json" --output "runs/traits/first/model"
```

查看原测量导出的 `*_graph.json` 或工作台中的节点后，指定真实根基/根冠所在节点，生成有方向的候选根轴：

```powershell
root-measurer phenotype "runs/default/results/model_delivery_local/run001/model_measurement.json" --root-base 0 --diameter-thresholds 0.2 0.5 1.0 --tangent-window 1.0 --diameter-offset 1.0 --diameter-window 1.0 --output "runs/traits/rooted/model"
```

`0` 只是示例，必须替换为该测量图中的节点编号。每个需要定向的树分量应有一个明确根基。根基可以是根冠上的公共起点；没有证据时不要为了产生数字而任意指定。

检查候选轴、分叉续接与端点类型后，可重用上一轮导出的标注配置，再确认指定分量：

```powershell
root-measurer phenotype "runs/default/results/model_delivery_local/run001/model_measurement.json" --annotations "runs/traits/rooted/model_annotations.json" --confirm-component 0 --output "runs/traits/reviewed/model"
```

`--confirm-component` 的编号来自该结果的 `architecture.components`。仍有未解决续接歧义的分量不能确认。每次使用新的输出前缀；已有目标文件会被拒绝覆盖。

其它选项包括 `--vertical-axis`、`--down-sign`、`--depth-bins`、`--origin X Y Z`、`--radial-thresholds`，以及接收完整设置对象的 `--settings settings.json`。位置方位分箱通过设置对象中的 `bearing_edges_deg` 修改。优先级为：内置默认值、导入标注配置中的设置、此次显式设置。

### Python API

```python
from root_measurer import phenotype

result = phenotype(
    "runs/default/results/model_delivery_local/run001/model_measurement.json",
    options={
        "diameter_edges": [0.2, 0.5, 1.0],
        "vertical_axis": "z",
        "down_sign": -1,
        "tangent_window": 1.0,
        "diameter_offset": 1.0,
        "diameter_window": 1.0,
    },
    annotation_profile="runs/traits/rooted/model_annotations.json",
)

for axis in result["axis_traits"]["axes"]:
    print(axis["axis_id"], axis["confirmation_status"], axis["length"])
```

`annotation_profile` 可以省略，也可传入已经解析的配置字典。API 返回数据；命令行和本地工作台提供保存及CSV导出。工作台保存同样生成新的分析快照，重新进入结果时可恢复与该测量文件匹配的最近一次标注。

## 2. 单位、默认设置与尺度

若测量文件含 `scale_mm_per_model_unit`，新增性状和分布的长度均换算为 `mm`；否则沿用测量文件声明的坐标单位。直径阈值、径向阈值、原点坐标、角度拟合窗口、粗度偏移和粗度窗口均使用这个**输出长度单位**，不会按节点数解释。

完整JSON同时保存用于溯源的原始坐标字段，不能对整份文件统一假定毫米：

| 字段范围 | 长度/坐标约定 |
| --- | --- |
| `distributions`、`axis_traits`及对应性状CSV | 使用结果声明的 `unit` / `output_unit`；有毫米尺度时已换算 |
| `architecture.axes[].points`及名称以 `_model_units` 结尾的架构字段 | 保留原图坐标，包括 `length_model_units`、`arc_positions_model_units` 和 `attachment_arc_model_units`；需要时乘 `scale_mm_per_model_unit` 转为毫米 |
| `measurement_coverage` 内沿用原测量表的长度字段 | 保留原模型坐标长度，使用该对象的 `length_unit` / `coordinate_units`；换算因子另列为 `scale_mm_per_model_unit` |

例如原坐标长度为2且比例为3 mm/模型单位，轴性状的 `length` 为6 mm，架构的 `length_model_units` 仍为2。界面高亮继续使用原坐标；“用当前根基节点设参考点”会将该节点换算为设置所需的输出单位。

| 设置 | 默认值 | 含义 |
| --- | --- | --- |
| `vertical_axis` / `down_sign` | `z` / `-1` | 竖直坐标轴与向下符号 |
| `diameter_edges` | `[0.5, 1.0, 2.0]` | 直径分级内部边界；自动补0与无穷上界 |
| `depth_bins` | `10` | 沿观测模型竖直跨度等分的层数 |
| `depth_edges` | 未指定 | 可显式给出非等宽深度边界，替代等分 |
| `inclination_edges_deg` | `0,15,…,90` | 无向线段与竖直轴的倾角区间 |
| `azimuth_edges_deg` | `0,15,…,180` | 水平无向轴线方位区间 |
| `bearing_edges_deg` | `0,45,…,360` | 相对显式原点的位置方位区间 |
| `radial_edges` | `[1.0, 2.0, 5.0, 10.0]` | 水平径向长度分布的内部边界 |
| `origin` | `null` | 位置方位及径向分布的原点；未提供时两类分布不可用 |
| `tangent_window` | `1.0` | 局部父/子切线各自使用的弧长窗口 |
| `diameter_offset` | `1.0` | 粗度窗口距分叉连接点的弧长偏移 |
| `diameter_window` | `1.0` | 母根/侧根粗度窗口长度 |
| `minimum_window_fraction` | `0.5` | 实际可用几何窗口相对请求窗口的最低比例 |
| `minimum_radius_coverage` | `0.5` | 粗度窗口内合格半径支持弧长的最低比例 |

这些默认直径区间是可修改的数值设置，**不把某个直径阈值自动解释为“吸收根”与“运输根”的边界**。跨样品比较应使用相同单位、同一组阈值与窗口，并记录输入分辨率。未做物理标定时，输出也不会被自动冠以毫米或真实田间根深。

## 3. 不依赖根轴身份的分布

所有长度以原图的唯一边为分母，重复经过同一条边不会增加总长度。分布同时给出已分配长度、未分配长度、相对全部观测轴长的比例、相对已分配长度的比例和覆盖率。

| 分布 | 计算对象与定义 |
| --- | --- |
| 直径级根长 | 合格半径样本的 `d = 2r`，按样本支持弧长 `arc_weight` 分组 |
| 直径 × 深度 | 按样本记录位置所在深度层分配其支持弧长；属于样本支持估计 |
| 倾角 | 无向三维线段与竖直轴的夹角，0°平行竖直轴，90°水平 |
| 水平轴线方位 | 投影方向按180°取模；与其反向的线段属于同一轴线方位 |
| 深度 × 方位 | 先将直线边在深度层边界精确分段，再按其水平轴线方位累计长度 |
| 相对位置方位 | 从显式原点到根段所在位置的水平角，分布覆盖0–360° |
| 深度 × 相对位置方位 | 在深度平面和位置方位边界射线处共同切分边，按各层空间扇区累计长度 |
| 水平径向长度 | 以显式原点所经过的竖直线为中心，将三维线段在水平环带内的弧长积分 |

区间为左闭右开；有限的最后一个端点包含在最后一档，无穷上界以 `null` 表示。恰落在直径内部阈值上的值进入较高一档。

直径只使用现有 `radius`、`available` 和截面资格状态。混合分支截面、带孔截面以及缺失/不合格截面均不进入单根直径统计；`connected_section_radius` 不会作为替代值。缺失弧长保留在分母中，也不会被插值补齐。

默认深度零点取观测实体边界的上侧；没有实体范围时使用已提取轴的范围。它描述当前模型坐标中的深度。洗根后的姿态变化、样本截取及原位土壤信息都不会由本模块自动恢复。竖直线段没有确定的水平方位，其长度单独列为未分配；未提供径向原点时，径向结果为不可用，不能把空值当成0。

`azimuth` 描述线段**朝哪个轴线方向延伸**，`bearing` 描述根段**位于参考点的哪个方向**，`radial` 描述根段**距参考竖直轴有多远**。例如，以Z为竖直轴时，位置方位为 `atan2(y − origin_y, x − origin_x)` 按360°取模；0°为参考点的+X方向，90°为+Y方向。沿+Y延伸但位于参考点东侧的根段，其切线方向与位置方向不相同。

位置方位对整条直线边按角度边界射线切分，联合深度时同时在深度平面处切分；每个已切开的区间计入相应空间扇区。不会用整条长边的中点代表它穿过的全部扇区。恰好位于参考竖直轴上的正长度根段没有位置方位，其长度作为未分配保留；其它位置的竖直根段有确定的位置方位，即使它没有水平切线方向。没有显式 `origin` 时，`bearing` 与 `depth_bearing` 均为空，不会自动把坐标原点或猜测的根基作为参考点。

原有实体体积和边界表面积继续遵循原测量器定义。这里的直径分布输出长度及比例，不从局部半径拼接管体来替代原实体体积或表面积。

## 4. 根轴身份与标注

### 从几何路径到有方向的根轴

一条连续生物根可能跨越多个图分叉。`architecture` 将根轴作为原边的一个完整分区：所有原边各归属一次，不删环、不接合原本分离的分量、不改坐标。

在已标注唯一根基的树分量中，自动方法根据局部分叉的延续方向提出续接候选。用户可以指定续接规则或确认轴身份。自动方法本身不推断物种、根类型、真实根尖，也不把每条图路径当成一条生物根。

| 状态 | 含义与可用性 |
| --- | --- |
| `CANDIDATE` | 已有明确候选轴解释；可查看候选性状，尚未人工确认 |
| `CONFIRMED` | 用户确认该轴；仍需独立准确率验证 |
| `UNRESOLVED` | 续接选择尚未解决；根数、母侧根性状、整轴方向、渐细和发育序汇总为空 |
| `GEOMETRY_ONLY` | 缺少根基、存在环或存在竞争根基；保留几何描述，不给出生物轴关系 |

对称Y形分叉在无法判断哪条支路延续母根时，不会先把母根截断再把两支都作为已经确定的侧根。临时图分段只用于显示与保留几何。零长度ID边也被保留，其长度贡献为0；不足以确定方向的轴返回明确不可用状态。

### 配置与源测量绑定

导出的 `*_annotations.json` 外层格式为：

```json
{
  "schema_version": "root-phenotype-annotations-v001",
  "measurement_sha256": "该测量文件的实际SHA256",
  "settings": {},
  "annotations": {
    "schema_version": "root-axis-annotations-v001",
    "root_bases": [],
    "junction_continuations": [],
    "endpoint_types": [],
    "root_types": [],
    "confirmed_axes": [],
    "confirmed_components": []
  }
}
```

建议从一次实际导出开始修改，保留正确的测量SHA。各列表记录格式：

| 列表 | 单条记录 |
| --- | --- |
| `root_bases` | `{"node_id": 0, "root_type": "primary"}`，类型字段可省略 |
| `junction_continuations` | `{"node_id": 12, "incoming_node_id": 11, "outgoing_node_id": 13}` |
| `endpoint_types` | `{"node_id": 20, "type": "root_tip"}` |
| `root_types` | `{"axis_id": "实际轴ID", "type": "用户定义类型"}` |
| `confirmed_axes` | 实际轴ID字符串 |
| `confirmed_components` | 分量编号 |

续接的 `outgoing_node_id: null` 表示用户明确判断原轴在此结束；它与没有选择续接不同。端点类型包括 `root_tip`、`cut_end`、`unknown`、`root_base`。几何表面端点不自动成为真实根尖。

轴ID绑定其原边集合与有向起点；改变基部、续接或测量文件后，旧确认不能无条件转移。SHA不匹配的导入配置会被拒绝。改变根基时应重新检查所有轴身份与端点标注。

## 5. 根轴性状的公式与空值规则

设母根轴观测弧长为 `L`，从基部开始的弧长坐标为 `s`，侧根连接位置排序后为 `s1,…,sN`。

### 分枝数量、密度与未分枝区

| 指标 | 定义 | 不可用情形 |
| --- | --- | --- |
| 侧根数 | 当前母轴的子轴数量 `N` | 未定向、环、续接未决 |
| 整轴LRBD | `N / L` | `L = 0` 或身份不可用 |
| 分枝区LRBD | `N / (sN − s1)`，首末连接点均计入N | 无侧根、仅一条侧根或首末位置相同 |
| 分枝间距IBD | 相邻连接点弧长差 `s(i+1) − si` | 少于两条侧根时为空列表 |
| 基部未分枝区/DTB | 已识别根轴基部到首次侧根的弧长 `s1` | 没有侧根或自然基部不能确定 |
| 根尖未分枝区 | `L − sN`，且末端明确标注为 `root_tip` | 没有侧根，末端为切口或未知 |
| 观测末端未分枝区 | 最后观测侧根到当前观测轴末端的弧长 | 没有侧根或身份不可用 |

两种LRBD明确输出 `numerator`、`denominator`、公式和单位。这里分枝区密度使用 `N/(sN−s1)`，并非间距均值的倒数；后者在常规情况下是 `(N−1)/(sN−s1)`。论文或实验采用固定长度窗口时，应使用一致的采样约定，不能直接把这些不同分母的结果混在一起比较。

没有侧根的明确有向轴，整轴LRBD为0，首末分枝位置及未分枝区保持NA。只有一条侧根时没有IBD观测，不能输出“间距0”。同一点连接多条侧根时可以出现真实的0间距，但不会除以0产生无穷密度。

IBD输出均值、中位数、分位数、总体标准差和CV；均值为0时CV为空。长度分布按轴等权，直径分布按合格支持弧长加权。分位数采用经验分布的逆累计函数，标准差为描述当前观测集合的总体标准差，不代表测量不确定度或置信区间。

### 局部萌生角、整轴方向和生长角

局部萌生角使用母轴连接点前方的入射段 `[s − w, s]` 与侧轴起始段 `[0,w]`，各自用连续线段弧长加权PCA拟合方向，再计算向量夹角。整轴PCA方向使用整条轴；母子整轴夹角和局部萌生角分开输出。

PCA使用连续线段的积分矩，不按顶点数量加权，因此仅把某段加密采样不会改变该段权重。向量符号依照基部到尖端的方向确定。拟合主方向不唯一或符号不确定时，结果为空。

每个局部窗口报告请求长度、实际有效长度和实际弧长范围。短轴允许使用达到最低窗口比例的实际范围，并标记窗口缩短；更短时返回 `INSUFFICIENT_TANGENT_WINDOW`。这些几何窗口设置不等同于已经验证足够的图像分辨率。

整轴生长角使用**有方向的**PCA向量与向下向量：0°向下、90°水平、180°向上。它与无根基也能计算的无向线段倾角不同，局部母侧根萌生角也不能用竖直角替代。

### 沿轴直径、渐细与母子粗度

`diameter_profile` 记录每个原始合格样本的弧长位置、支持区间、直径、资格状态和源边。半径缺失的区间保持缺失。沿轴回归为：

```text
d(s) = intercept + slope × s
taper_rate = −slope
```

回归以样本支持弧长为权重，至少需要三个不同的合格弧长站位。斜率与 `taper_rate` 的单位为直径单位/弧长单位，当前同单位时为无量纲；正的 `taper_rate` 表示变细。输出R²、有效站位范围、有效弧长覆盖率和缺失状态；恒定直径的回归R²为空，不能用1来冒充解释了不存在的变异。未确认的明确候选轴可查看候选回归，续接未决轴不回归整轴渐细。

母子粗度在分叉旁留出偏移 `o` 后测量：

```text
母根窗口：[s_attach − o − w, s_attach − o]
侧根窗口：[o, o + w]
粗度比：侧根窗口弧长加权平均直径 / 母根窗口弧长加权平均直径
```

使用单根合格截面，排除混合截面；报告实际窗口、合格支持弧长及覆盖率。窗口不足、覆盖率不足、母根平均直径为0时，粗度比为空。该比值是形态描述，未被解释为水力导度比例。

## 6. 根类型、发育级次与Strahler序

根类型只来自显式标注，未标注为 `unknown`。类型和发育级次的轴数、总长、轴长分布、合格直径分布按确认状态分组；候选与确认结果不合并成“全部已确认”。

发育级次采用本模块明确约定：基部轴为0，侧轴为母轴级次加1。这是从基部向外的层级，与“胚根起源的primary root”这一植物器官类型不能混为一谈。只确认某条子轴、但尚未确认其祖先关系时，已确认发育级次仍为空。

Strahler序从末端边的一阶向基部递推：最大子序仅出现一次时取该最大值；多个子序共享最大值时加1。**一条连续轴可以在不同位置具有不同的Strahler序**，因此本模块保留逐边序，并单独提供明确标记的轴基部边序。

Strahler统计按边所属序累计长度和直径支持。`geometric_segment_count` 只是图线段数；细分一条线段可改变此数，不能把它叫作该根序的生物根数。原图有明确根基且无环时，Strahler候选边序不依赖哪条支路被选为母根续接；确认字段仍与候选字段分开。

## 7. 输出与复现

一个输出前缀生成七个文件：

| 后缀 | 内容 |
| --- | --- |
| `_phenotyping.json` | 完整结果、来源、设置、架构与所有性状 |
| `_annotations.json` | 与测量SHA绑定、可再次导入的设置和标注 |
| `_distributions.csv` | 直径及空间分布、未分配长度与覆盖率 |
| `_axes.csv` | 每轴身份、状态、长度、分枝、角度、渐细等 |
| `_branches.csv` | 母侧轴关系、局部角、有效窗口与粗度比 |
| `_profiles.csv` | 沿轴合格/缺失半径与直径样本 |
| `_summaries.csv` | 根类型、发育序、逐边Strahler序及确认状态分组 |

CSV使用UTF-8 BOM，保存方法版本、测量SHA、设置SHA、标注SHA、输出单位、方向设置和提取完整性标识。源测量与原实体/网格绑定在分析前后核验。候选结构、源测量状态、数值覆盖率和参考准确率是不同字段。

同一样品比较设置时应保留两个分析快照；跨样品比较时应统一单位、方向、分组边界与窗口，并保留生育阶段、物种、处理与取样完整性等实验记录。本模块不会从几何自动补出这些实验事实。

## 8. 当前范围与论文依据

本次扩展不计算需要额外数据的SRL、SRA、RTD等干质量指标，不生成土体距离场或土壤体积归一化的真实RLD，不推断生长时间序列、伸长速率，也不自动识别玉米节根轮次。吸水率、养分获取效率、水力导度、抗倒伏力及吸收根功能比例均未由这些几何量直接推定。

相关原始研究用于选择有解释价值的参数、澄清定义；它们不构成本程序已通过同类实验验证的证据：

1. **York & Lynch, 2015.** *Intensive field phenotyping of maize (Zea mays L.) root crowns identifies phenes and phene integration associated with plant growth and nitrogen acquisition.* Journal of Experimental Botany 66:5493–5505. [DOI:10.1093/jxb/erv241](https://doi.org/10.1093/jxb/erv241)。研究使用侧根密度、侧根长度、首次分枝距离等性状；本文档的整轴及分枝区分母是另外明确记录的统计约定。
2. **Zhan, Schneider & Lynch, 2015.** *Reduced Lateral Root Branching Density Improves Drought Tolerance in Maize.* Plant Physiology 168:1603–1615. [DOI:10.1104/pp.15.00187](https://doi.org/10.1104/pp.15.00187)。特定玉米材料在缺水条件下，少而长侧根表型显示优势。
3. **Zhan & Lynch, 2015.** *Reduced frequency of lateral root branching improves N capture from low-N soils in maize.* Journal of Experimental Botany 66:2055–2065. [DOI:10.1093/jxb/erv007](https://doi.org/10.1093/jxb/erv007)。低氮试验也支持联合考察侧根数量与长度。
4. **Jia, Liu & Lynch, 2018.** *Greater lateral root branching density in maize improves phosphorus acquisition from low phosphorus soil.* Journal of Experimental Botany 69:4961–4970. [DOI:10.1093/jxb/ery252](https://doi.org/10.1093/jxb/ery252)。低磷条件下，多而短侧根表型更有利，说明“分枝越多越好”或“越少越好”都不是通用评分规则。
5. **Uga et al., 2013.** *Control of root system architecture by DEEPER ROOTING 1 increases rice yield under drought conditions.* Nature Genetics 45:1097–1102. [DOI:10.1038/ng.2725](https://doi.org/10.1038/ng.2725)。稻DRO1研究支持根生长方向和深根的农学价值；其生长角与本模块的局部母侧根夹角是不同量。
6. **Guo et al., 2008.** *Anatomical traits associated with absorption and mycorrhizal colonization are linked to root branch order in twenty-three Chinese temperate tree species.* New Phytologist 180:673–683. [DOI:10.1111/j.1469-8137.2008.02573.x](https://doi.org/10.1111/j.1469-8137.2008.02573.x)。根序与解剖功能的关系支持分组研究，同时提醒统一直径阈值不能替代功能分类；木本研究的具体功能范围不能直接照搬到所有作物。
7. **Zeng et al., 2021.** *TopoRoot: a method for computing hierarchy and fine-grained traits of maize roots from 3D imaging.* Plant Methods 17:127. [DOI:10.1186/s13007-021-00829-z](https://doi.org/10.1186/s13007-021-00829-z)。按层级报告根数、长度、粗度、子根数及不同位置的角度，支持本模块按已辨别轴和层级统计的参数选择。本模块保留原图全部边，不复制其拓扑简化、断环或重建修复流程；本程序的候选轴规则及准确率不能等同于该论文方法。
8. **Herrero-Huerta, Raumonen & Gonzalez-Aguilera, 2022.** *4DRoot: Root phenotyping software for temporal 3D scans by X-ray computed tomography.* Frontiers in Plant Science 13:986856. [DOI:10.3389/fpls.2022.986856](https://doi.org/10.3389/fpls.2022.986856)。其静态根型、层级及直径/深度/方向分布为参数类别提供参考；本模块使用唯一原边长度与合格截面半径，不采用拟合圆柱的体积/表面积去替代原实体测量，也未实现该论文的时间序列分析。本文的无向切线方位和相对原点的位置方位各有明确约定，不能仅因都称“方位”就直接对比数值。
9. **Pagès, 2019.** *Analysis and Modeling of the Variations of Root Branching Density Within Individual Plants and Among Species.* Frontiers in Plant Science 10:1020. [DOI:10.3389/fpls.2019.01020](https://doi.org/10.3389/fpls.2019.01020)。研究分析IBD分布，并在已出现可见侧根的年轻分枝区测量，排除尚在出根的远端区。它支持保留IBD及其变异而非只报一个密度；本模块明确区分观测整轴的 `N/L`、首末侧根区间的 `N/(sN−s1)` 和平均IBD倒数，避免把不同采样区间、端点计数约定混为同一性状。这里不拟合其潜在分枝位点/侧根萌发概率模型。
10. **Symonova, Topp & Edelsbrunner, 2015.** *DynamicRoots: A Software Platform for the Reconstruction and Analysis of Growing Plant Roots.* PLOS ONE 10:e0127657. [DOI:10.1371/journal.pone.0127657](https://doi.org/10.1371/journal.pone.0127657)。论文用整条分枝体素集合的主PCA轴描述对重力角和母子夹角，明确讨论了“取根的哪一段”会影响角度定义。本模块分别输出局部窗口角与整轴PCA角，并使用连续轴线弧长加权而非体素集合加权；二者描述不同几何对象，不能把整轴角当作分叉处的萌生角。本次实现不包含其时间函数、逐时配准、生长分析或随时间修正根层级的流程。

本模块的方法验证包括解析几何、尺度换算、节点细分不变性、重复边与缺失样本、半边样本反向映射、续接歧义、零长度边、截断端和候选/确认语义。它验证计算约定与软件行为；真实材料上的准确率与不同生物处理的区分能力仍应由独立参考数据及实验设计评估。
