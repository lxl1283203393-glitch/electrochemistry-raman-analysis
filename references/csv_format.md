# Optosky CSV 与实验状态

导入器扫描目录内全部 CSV，自动寻找实际光谱表头，不固定起始行。支持编码差异、空白行、大小写与列顺序变化。规范化通道为 `raw`、`dark`、`dark_subtracted`、`baseline_subtracted`，波数为 `raman_shift_cm1`。可用原通道保留；缺失所需通道要报告，不生成虚构列。

元数据原文和原始行保留，并提取内部采集时间、积分时间、激光功率与平均次数。文件系统修改时间不能自动当作采集时间。内部时间缺失时，时间序列有证据缺口。共享波数轴、重复波数、非数字值、缺失值、异常范围和 Raw−Dark 对应关系均需核查。不要假定仪器的 Dark Subtracted 恰好等于逐列相减。

每条记录包含路径、批次、稳定 ID、SHA256、编码、波数轴、原始通道、元数据、实验状态、电位、内部时间、警告、处理阶段与历史。`to_ramanspy(channel=...)` 生成 RamanSPy 光谱对象；依赖不可用应明确失败。

| 文件名例子 | 解析状态 | 电位含义 |
| --- | --- | --- |
| `-1.023cv` | `applied_potential` | 原符号 −1.023 V |
| `-1.353ing cv` | `applied_potential` | −1.353 V；`ing` 操作含义未知 |
| `0cv` | `applied_potential` | 0 V，不能称无外加偏压 |
| `OCP` / 明确开路标记 | `open_circuit` | 非指定施加电位 |
| 明确无外加偏压标记 | `unbiased` | 不能假设等于 0 V |
| `no sample` | `blank` | 空白，类型 `no_sample` |
| 不明确名称 | `unknown` | 不猜测 |

`blank` 可由配置指定 `bare_ni_foam/electrolyte/optical_window/other`，保留来源。电位名称缺少负号时不根据 HER 背景自动补负号。只有用户或原始实验记录证实时才能在 `filename_overrides` 给出 `potential_v/state/reason`；原标签及覆盖来源均保留。

RHE 换算不替换原始电位。参比、内充液、pH、温度与 iR 状态不足时保留原始参比；明确启用近似且具备可解释条件时，输出假设而非伪装为精确校准。不能从 1 M KOH 的文字默认补成已实测 pH。
