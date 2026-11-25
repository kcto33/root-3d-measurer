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

运行 `root-measurer demo --output runs/demo_inputs --case straight` 可生成完整可读示例。示例不会覆盖已有目录。

## 实体化与测量的职责

实体化阶段负责按同一套已登记规则输出最终测量对象：去除伪影、内部空腔、离散碎片和薄片，处理异常连接，并保存保留、剔除、填补、断开或连接位置及理由。原始重建与处理账本在实体化侧保留；它们不作为待测材料自动混入最终对象。若形成派生实体包，`source_negative_indices`表示该包的源材料，`regularized_negative_indices`仍须保留该包全部源材料；跨包的剔除另存账本，不能破坏现有字段契约。

测量器接收并测量冻结后的输入实体。它不判断该实体是否包含真实根系的全部根段，不从原根缺失推断长度，也不自行删碎片、补洞或连接断根。一个截断但边界明确的输入实体，可以测量其自身几何参数。

参数输出中的`measurement_scope`明确对象为`supplied_finite_entity`，并标记`original_root_completeness_assessed=false`。兼容字段`complete_model_length_available`只表示当前算法是否覆盖输入实体的分量和几何端点；它不是原根系完整性判断。若算法尚未求得输入实体某段的轴或截面，应输出具体计算限制和已测区间，不把未解值填零，也不把部分轴长标成输入实体总长。体积、面积等独立可计算字段不因轴失败而停止。

重建是否完整、清理是否误删真根、输出是否接近原始根系，属于实体化验证或独立重建评价；这些问题不作为测量器对当前输入实体计算参数的前提。
