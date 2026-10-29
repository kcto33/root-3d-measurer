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
