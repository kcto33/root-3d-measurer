# 原生实体输入规范

一个输入目录包含 `entity_prediction.json` 和该文件 `fields[].cache` 所绑定的 NPZ 文件。多个裁块必须共享格距，重合占据不能重复计数。

JSON 至少包含 `root`、`fields`、`duplicate_occupied_cells`（false）、`all_source_negative_cells_accounted`。每个 cache 记录文件名/路径、SHA-256 和字节数，工具从 JSON 旁按文件名读取并复核内容。`storage.pitch` 用于页面显示；实际格距以校验后的 NPZ 为准。

每个 NPZ 的字段恰好为：

| 字段 | 约定 |
|---|---|
| `source_negative_indices` | N×3、唯一的整数局部占据索引 |
| `regularized_negative_indices` | M×3、唯一整数索引，保留全部 source 材料 |
| `pitch` | 正的有限标量，所有裁块一致 |
| `origin` | 三维第一格中心，满足 `(integer_offset + 0.5) * pitch` |
| `shape` | 三个正整数，局部索引在范围内 |

完整实体是索引对应的闭合立方体并集，全球格中心为 `(global_index + 0.5) * pitch`。输入不提供中心线/半径答案；这些由工具计算。`reference_graph.json` 仅用于可选的独立评分，预测不会读取它。

物理尺度不由归一化坐标推断。每模型单位对应毫米的比例必须由用户确认。面积按比例平方、体积按比例立方换算。

## 外部体素与采集元数据

支持三维 NPY、NPZ、TIFF/OME-TIFF。数组保持文件原轴序；`auto` 使用采集元数据，普通 TIFF 的 QYX/IYX 页堆栈按 ZYX 解释，其余无轴序信息时使用 XYZ。也可显式指定 `xyz` 或 `zyx`。多数组 NPZ 必须能唯一确定三维数组，或给出 `array_key`。

NPZ 可保存标量字符串 `axis_order`、`unit` 及三元素 `spacing`、`origin`；这些键不参与三维数组候选选择。OME-TIFF 读取 PhysicalSizeX/Y/Z 及其已声明单位。可选采集 JSON 示例如下，显式 JSON 覆盖内嵌元数据，显式参数再覆盖 JSON：

```json
{"spacing": [0.05, 0.05, 0.20], "axis_order": "zyx", "unit": "mm", "origin": [10, 20, 30]}
```

`spacing` 和 `origin` 始终为 XYZ，单位支持 mm、cm、m、um/µm/μm、nm。原点仅作为来源坐标信息保存；测量采用以体素零号角为原点的局部 XYZ，不自动施加扫描仪旋转或平移。未声明单位时保留模型单位；声明单位可自动换算为毫米，显式毫米比例优先。该信息保存在派生实体包中，后续直接测量此包时仍生效。

不指定标签时取大于阈值的前景；指定 `labels` 时仅取这些整数标签的并集，忽略阈值。同一文件的不同标签可分别列入批次。原始文件只读。

内部测量域为等距立方体。输入三轴间距不等时，必须显式给出 `resample_pitch`；按目标格中心最近邻采样，不默认猜测目标分辨率。JSON 的 `resampling` 保存源标签体积、输出体积、方法及原文件保留状态；参数和测量报告的 `input_preparation`/`volume_import` 保留这些信息。重采样可能改变体积、薄枝和接触关系，结果对应派生实体。

运行 `root-measurer demo --output runs/demo_inputs --case straight` 可生成完整可读示例。示例不会覆盖已有目录。

## 实体化与测量的职责

实体化阶段负责按同一套已登记规则输出最终测量对象：去除伪影、内部空腔、离散碎片和薄片，处理异常连接，并保存保留、剔除、填补、断开或连接位置及理由。原始重建与处理账本在实体化侧保留；它们不作为待测材料自动混入最终对象。若形成派生实体包，`source_negative_indices`表示该包的源材料，`regularized_negative_indices`仍须保留该包全部源材料；跨包的剔除另存账本，不能破坏现有字段契约。

测量器接收并测量冻结后的输入实体。它不判断该实体是否包含真实根系的全部根段，不从原根缺失推断长度，也不自行删碎片、补洞或连接断根。一个截断但边界明确的输入实体，可以测量其自身几何参数。

参数输出中的`measurement_scope`明确对象为`supplied_finite_entity`，并标记`original_root_completeness_assessed=false`。兼容字段`complete_model_length_available`只表示当前算法是否覆盖输入实体的分量和几何端点；它不是原根系完整性判断。若算法尚未求得输入实体某段的轴或截面，应输出具体计算限制和已测区间，不把未解值填零，也不把部分轴长标成输入实体总长。体积、面积等独立可计算字段不因轴失败而停止。

重建是否完整、清理是否误删真根、输出是否接近原始根系，属于实体化验证或独立重建评价；这些问题不作为测量器对当前输入实体计算参数的前提。
