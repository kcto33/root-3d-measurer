# Root 3D Measurer · 根系三维测量器

在本机输入三维根系实体，输出唯一中心线总长、分支路径长度、沿线截面等效半径/直径、体积、边界面积和包围盒。浏览器中可点选端点分支、在同一分支选择两点、描绘有序相连根段，以及在连续轴位置现场查询原实体截面。

本仓库是独立迭代的工具。运行不需要原研究项目、历史实验目录、数据集、GPU、云账户或外部网络服务。当前版本为 `0.3.2`，既有数值内核保留可追溯来源；便携运行层和公共 API 单独维护。

大型连通原生实体使用占据邻接存储执行相同Lee细化顺序与判据，避免分配整个空白包围盒；小实体保留原稠密实现。原始格距、材料、图解释、截面及端点规则不变，资源上限按实际分配核验。完整参考与测量精度仍独立评价，计算可运行不等于全部字段精度合格。

## 安装

Python 3.11 或更新版本。推荐新建虚拟环境，再在仓库目录运行：

```powershell
python -m venv .venv
.venv\Scripts\Activate.ps1
python -m pip install -e ".[dev]"
root-measurer --version
```

Linux/macOS 使用 `source .venv/bin/activate`。依赖为 NumPy、SciPy、scikit-image、Shapely 和 psutil；无需终端以外的编译步骤或 Node.js。scikit-image 限定在已验证的 0.26 系列，因为部分现有内核需要其可检查的 Python 实现。

安装完成后，Windows 也可双击 `scripts/start_windows.cmd` 启动工作台。

## 快速开始

```powershell
# 生成一根可重复的小型合成直根；不是独立验收集。
root-measurer demo --output runs/demo_inputs --case straight

# 测量。尺度由用户确认；留空时全部使用模型单位。大实体可用 --section-workers N 启用多进程截面查询（结果与单进程严格一致）。
root-measurer measure runs/demo_inputs/straight --workspace runs/demo_workspace --scale-mm 100 --section-workers 1

# 浏览器工作台；输入/输出目录均由用户指定，安装目录不保存运行结果。
root-measurer gui --workspace runs/demo_workspace --input runs/demo_inputs/straight --open-browser
```

工作台只监听 `127.0.0.1`，默认由系统分配可用端口。命令会打印实际地址。测量进程禁用自身网络连接，浏览器界面通过本机 HTTP 访问。每次运行新建 `run001`、`run002` 等目录；已有结果不覆盖。

## 输入

**主要输入：原生实体包。** 选择 `entity_prediction.json` 或所在目录，旁边保留其 SHA 绑定的 NPZ 裁块。支持保存的原始占据/实体化占据及原格距；坐标约定与字段说明见 [输入规范](docs/input_format.md)。不需要把原实体转换成 PLY。

实体化阶段负责统一清理伪影、空腔、碎片、薄片和异常连接，确定最终测量对象。测量器只计算该输入实体的参数，不推测真实根系是否完整。轴长的可用状态表示算法对当前输入的计算覆盖；未解轴段或截面按字段报告，参见[职责与输出范围](docs/input_format.md#实体化与测量的职责)。

**通用三维二值体素输入（`.npy` / `.npz` / 3D `.tif` / `.tiff`）：** 适用于 X 射线 CT、MRI 或任意外部三维分割结果，直接导入或测量而无需预转换：

```powershell
# 直接测量三维体素文件（CT/TIFF 堆叠常用 --voxel-axis-order zyx；pitch 为体素边长）
root-measurer measure root_ct.npy --voxel-pitch 0.05 --voxel-axis-order zyx --scale-mm 1 --section-workers 4

# 或先转为 SHA-256 绑定的原生实体包再复用/批量测量
root-measurer import-volume root_ct.tif --output runs/ct_solid --pitch 0.05 --axis-order zyx
```

可指定分割标签（如 `--voxel-labels 1 3`）；标签模式优先于阈值。OME-TIFF 的物理间距／单位和 NPZ 的 `spacing`、`axis_order`、`unit`、`origin` 可自动读取，也可用 `--voxel-metadata acquisition.json` 提供。间距向量始终为 XYZ，与数组存储轴序分别设置。

```powershell
# 非等距体素：显式选择目标格距；按最近邻转换，保留原文件及体积变化记录。
root-measurer measure root_ct.npy --voxel-spacing 0.05 0.05 0.20 --voxel-axis-order zyx --voxel-unit mm --voxel-labels 1 --voxel-resample-pitch 0.05
```

未指定目标格距时，各向异性体素会明确拒绝转换。有可信单位时自动换算为毫米；显式 `--scale-mm` 优先。先 `import-volume` 再测量其实体包也保留同样的尺度。重采样后的参数描述派生实体，源原点保存在元数据中，测量坐标使用局部 XYZ。格式与约定见[输入规范](docs/input_format.md)。

TIFF 使用原始堆栈轴序读取，不按 RGB/RGBA 图像自动换轴。通过 `measure` 或工作台测量体素时，转换与后续测量均在受监督工作进程内完成；内存上限统计整个工作进程树（含并行截面进程），取消或超限会一起停止这些进程。独立 `import-volume` 只执行格式转换，不启动测量监督器。

**兼容输入：闭合三角网格 PLY。** 使用 `measure model.ply --mesh-target-pitch ... --mesh-target-reason "依据输入特征尺度"`。原三角面用于半径截面，派生占据仅用于轴提取。分辨率和相位挑战未通过时明确阻塞；没有宣称所有外部 PLY 均能完整测量。

## 批量测量

创建 `inputs.json`：相对路径以该清单所在目录为基准。

```json
{"inputs": ["demo_inputs/straight", "demo_inputs/curve"]}
```

```powershell
root-measurer demo --output runs/demo_inputs --case curve
root-measurer batch runs/inputs.json --workspace runs/batch_workspace
```

工作台、CLI 和 Python API 支持混合原生实体、体素和闭合 PLY。清单既兼容路径字符串，也接受逐样本设置：

```json
{
  "defaults": {"section_workers": 2},
  "inputs": [
    {"id": "root_a", "path": "root_a.npy", "voxel_pitch": 0.05, "voxel_unit": "mm", "voxel_labels": [1]},
    {"id": "root_b", "path": "root_b.npz", "vertical_axis": "y", "down_sign": 1},
    {"id": "root_c", "path": "root_c.ply", "mesh_target_pitch": 0.05, "mesh_target_reason": "输入最小特征尺度"}
  ]
}
```

路径及 `voxel_metadata` 相对清单文件解析；逐项设置优先于清单 `defaults`，清单 `defaults` 优先于命令行默认值。同一文件可用不同标签重复列入。输入失败后继续其余样本，取消项保留完整分母；每次使用新运行目录。`batch_summary.json/csv` 每模型一行，`batch_parameters.csv` 保存全部参数、单位、状态和覆盖率，失败模型也保留缺项行。工作台重启后可查看历史汇总、下载参数；没有结束汇总的批次显示待确认状态，已有进度记录可以导出，不会自动重启。

Python：`from root_measurer import batch`，调用 `result = batch("inputs.json", workspace="runs/batch_workspace")`，通过 `result.summary` 查看完整结果。

## 植物形态指标

新增17项几何形态指标：三维中心线有限尺度盒计数维数、深度、最大水平展幅、展深比、凸包体积、体积充实度、水平材料投影面积／凸包面积／充实度、凸包内根长密度、几何分叉密度、路径曲折度、竖直夹角、面积体积比、球形度、根长深度D50／D95；另提供10层根长和原实体体积深度分布。

默认 **Z轴竖直、负方向向下**。工作台、结果页、单模型和批量命令均可切换。深度零点为实体最上边界，未假定土壤表面。分形维数显示所有固定二倍尺度、8种网格偏移计数、拟合R²及适用状态；尺度不足不输出数值。形态描述与独立精度验收分开。定义、单位和限制见 [形态指标](docs/morphology.md)。

```powershell
root-measurer measure runs/demo_inputs/straight --workspace runs/demo_workspace --vertical-axis z --down-sign -1
# 对已有测量增量分析；不重新提取轴、不改写历史文件，输出前缀必须未使用。
root-measurer morphology runs/demo_workspace/results/model_delivery_local/run001/model_measurement.json --output runs/demo_workspace/analysis/straight_z_down --vertical-axis z --down-sign -1
```

Python：`from root_measurer import morphology`，再调用 `morphology(measurement_path, vertical_axis="z", down_sign=-1)`。原参数JSON/CSV记录测量时方向；结果页的“形态指标CSV／完整JSON／深度分布CSV”使用当前页面所选方向。批量汇总增加主要形态列、单位、方向和每项资格状态。

## 根系表型与根轴标注

对已有测量增加以下六组参数；无需重新提取轴。默认沿用已记录的物理尺度，未标定时明确使用模型坐标单位。

| 参数组 | 新增结果 |
| --- | --- |
| 直径分级 | 各直径区间的根长、比例、有效覆盖率及直径×深度分布 |
| 根轴与根序 | 指定根基、检查分叉续接，形成连续根轴与母侧根关系；按人工根类型、发育级次和逐边Strahler序汇总 |
| 侧根配置 | 侧根数、全母根/分枝区两种密度、沿母根的分枝间距、首次分枝距离和未分枝区 |
| 角度 | 局部母侧根夹角、整轴PCA夹角、整轴相对向下方向的生长角 |
| 沿轴粗度 | 直径曲线、带覆盖率的渐细回归、避开分叉的母侧根直径比 |
| 空间分布 | 无向根段倾角/朝向、深度×朝向，以及相对显式参考点的位置方位、深度×位置方位与水平径向根长 |

工作台结果页新增“根系表型与根轴标注”：先查看几何分布，再点选或输入根基节点，核对续接和根尖/截断端，最后确认相应分量。可以修改分箱与物理窗口、在模型上高亮根轴、导入标注、导出JSON/CSV及保存分析。重新打开结果时恢复最近一次与原测量SHA匹配的已保存标注。

复核清单显示待核对根基、未决分叉、未人工确认根轴和未解几何端点的数量；`review_summary` 随完整 JSON 及 CSV 输出。短枝、闭环、分离分量和缺失半径的可用部分继续导出，缺项不插补。改根基或续接后需要重新确认相应根轴。

**自动续接只提供候选。** 对称分叉等未决关系保持不可用；闭环和分离分量保留原状。人工确认表达根轴身份判断，独立精度验收另行记录。直径只使用合格单分支半径，混合/带孔/缺失截面不补值。各分布保留未分配弧长及完整分母。

```powershell
# 先计算无需根轴身份的分布，全部文件另存。
root-measurer phenotype runs/demo_workspace/results/model_delivery_local/run001/model_measurement.json --output runs/demo_workspace/traits/first
# 在结果页或图导出中核对节点后指定根基；0仅为示例。
root-measurer phenotype runs/demo_workspace/results/model_delivery_local/run001/model_measurement.json --root-base 0 --output runs/demo_workspace/traits/rooted
```

Python入口：`from root_measurer import phenotype`，调用 `phenotype(measurement_path, options={...}, annotation_profile=...)`。标注可省略，或传入导出的标注JSON路径/字典。CLI和工作台每次保存完整分析、标注配置以及分布、根轴、母侧根关系、沿线粗度、分组汇总五类CSV，共七个文件；已有文件拒绝覆盖。

定义、分母、单位、标注格式和论文依据见 [表型分析说明](docs/phenotyping.md)。当前范围不包括需要干质量的SRL/SRA/RTD、原位土体距离场、时间序列生长以及玉米节根轮层；也不从几何直接推断吸水或养分获取效率。

## 交互和 Python API

浏览器中选择操作方式后点选或描绘。粉色叠加表示实际计长范围或新计算的原截面。遮挡处先旋转/放大，或限定已选择的路径。系统只计量用户明确选择的根段，不搜索替代连接路径。

```python
from root_measurer import measure, query_length, query_radius

result = measure("runs/demo_inputs/straight", workspace="runs/api_workspace")
print(result.parameters)
segment = query_length(result.measurement_path, [
    {"path_id": "path00000", "start": 0.2, "end": 0.8}
])
section = query_radius(result.measurement_path, "path00000", 0.5317)
```

API/命令行区间与 `arc` 使用**模型单位的路径弧长**；结果根据记录的尺度输出 mm 或模型单位。浏览器输入框使用页面显示单位。命令行亦提供 `length`、`radius` 和独立参考 `score`，可用 `--help` 查看。

## 输出和验收

- `_measurement.json`：原始轴、端点、截面和异常记录。
- `_parameters.json/csv`：整模型参数；路径、沿线、图及原边界显示文件同时生成。
- `_supervision.json`、进度与日志：用时、内存预算、取消/阻塞状态。
- 现场查询 JSON：在对应运行的 `interactive_queries/` 保存，绑定测量和输入哈希。
- `EXPERIMENT_LOG.md` 与 `snapshot/`：可复核的操作和配置记录。

主交互半径是**原实体整个连通材料截面的面积等效半径** `sqrt(A/pi)`，孔洞扣除，混合/分叉截面标注。严格单分支半径作为单独字段保留。预测输出有效率和完整参考覆盖率使用不同分母，不能把覆盖率当成精度。

`score` 使用独立参考轴和原有限实体截面，预测不读取参考文件：

```powershell
root-measurer score runs/demo_workspace/results/model_delivery_local/run001/model_measurement.json --reference runs/demo_inputs/straight/reference_graph.json --output runs/demo_workspace/analysis/straight_score.json
```

功能完成、几何端点可用和独立精度通过分别记录。没有参考时显示“未评分”。接触/非唯一截面、未解决端点和无唯一轴的分量不补值、不改成通过。历史批次结论及适用范围见 [验证基线](docs/validation.md) 和 [已知限制](docs/known_limitations.md)。

## 开发

```powershell
python -m pytest
python -m build
```

公共接口位于 `src/root_measurer/api.py`、`cli.py` 和 `runtime.py`。现有内核的文件身份/源码读取机制保留在包内 `src/`、`scripts/`，属于内部实现。源码及网页资源都进入 wheel，安装后无需原仓库。来源映射见 [source_provenance.json](docs/source_provenance.json)。

首版迭代任务见 [路线图](docs/roadmap.md) 和 [变更记录](CHANGELOG.md)。仓库当前为私有，尚未选择开源许可证。
