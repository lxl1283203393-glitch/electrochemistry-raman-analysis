---
name: electrochemistry-raman-analysis
description: Analyze Optosky Raman Microscope CSV exports with executable RamanSPy preprocessing, quality audits, candidate bands, constrained mathematical fits, and traceable electrochemical potential or batch comparisons. Use for folder-based in situ Raman analysis; chemical assignments and mechanisms require separate evidence review.
---

# Electrochemistry Raman analysis

使用本目录的 Python 工具完成计算；安装和四种命令见 [README.md](README.md)。不得只给出提示词或凭图估算结果。输入目录、样品、参比及处理参数由用户数据和配置确定，不使用内置样品结论。

1. 检查输入文件和 `assets/experiment_config.yaml`。未知实验信息保持 `null`；首先阅读 [CSV 与状态规则](references/csv_format.md)，保留文件名电位符号。明确用户要求 `audit`、`analyze`、`compare` 或已有结果的 `report`。
2. 首次遇到一个目录时执行 `audit`，核对解析失败、元数据、内部采集顺序、Raw−Dark 一致性与异常标记。原始 CSV 只读，输出到另一个目录；可疑文件保留。
3. 执行 `analyze` 或 `compare` 前按 [预处理方法](references/preprocessing_methods.md) 核对 A/B/C 路径、基线和平滑敏感性、设置可比性及空白条件。拟合只用于配置窗口，`unidentifiable/failed/no_peak` 不进入可靠峰趋势。PCA/NMF 仅在质量审计后按需开启。
4. 查看输出表、处理历史、图表清单及报告；实际打开关键 PNG，检查原始全谱、可疑光谱和拟合残差。图表清单的 `skipped` 要解释原因，`error` 属于失败。候选峰、基线及拟合的合并面板受配置数量上限控制，完整候选与拟合结果保留在表中。
5. 按 [科研规则](references/scientific_rules.md) 和 [解释边界](references/interpretation_guidelines.md) 区分直接观察、支持性解释、待验证假设与无法判断。每项计算发现关联文件、处理方法、参数、表和图；没有独立重复只做描述性分析。

可移植入口：`python scripts/main.py --mode analyze --input INPUT --config assets/experiment_config.yaml --output OUTPUT`。安装后也可用 `electrochemistry-raman`，兼容命令 `optosky-raman` 保留。`report --results EXISTING --output REPORT` 重建报告，不重新处理光谱。配置字段详见 [配置说明](references/configuration.md)。

遇到错误时保留已有输出，报告具体文件和原因；不悄悄跳过 RamanSPy、改写数据、猜测电位或强行得到拟合。检查退出码及错误表后修正输入/配置并重跑。未取得真实原始文件时，将验证明确标为模拟数据的软件测试。缺参比、内充液、pH、温度或 iR 信息时，保留原始参比；近似 RHE 换算必须注明假设。

软件依据和 API 核查记录见 [软件来源](references/software_sources.md)。实际测试环境与交付边界见 VALIDATION.md。对相关残差导致的不确定度抑制和 `model_adequacy: not_established`，即使数学拟合状态为 ok 也不能宣称峰参数已经充分验证。
