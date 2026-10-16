# v1.0.1 — Electrochemistry Raman analysis

本补丁版本统一项目命名，不改变光谱处理算法。

- GitHub 仓库和 Skill 标识改为 `electrochemistry-raman-analysis`。
- 显示名、README、报告标题、版本记录与安装包同步更新。
- 新 CLI 为 `electrochemistry-raman`；保留 `optosky-raman` 兼容入口，Python 内部模块保持 `optosky_raman`。
- 108 项自动测试通过，完整 ZIP 包含新 Skill 根目录及 SHA256 校验。
- v1.0.0 保留作历史版本。

## v1.0.0 — 初始版本

首个公开发布版本，提供可安装的 Codex Skill 和基于 RamanSPy 的可执行 Python 分析工具。

- 四种 CLI 模式：audit、analyze、compare、report。
- 动态识别 Optosky CSV 表头、编码、原始通道、元数据、采集时间及实验状态。
- RamanSPy 光谱对象与预处理；保留仪器 Dark Subtracted、BaseLine Subtracted 以及统一处理路径。
- AsLS/arPLS/airPLS、平滑、尖峰处理、条件核查后的可选空白扣除与参数敏感性分析。
- 保守候选峰检测；Gaussian/Lorentzian/pseudo-Voigt 模型比较、残差与不可辨识诊断。
- 电位/时间混杂、同电位和跨批次描述性比较，可选 PCA/NMF。
- 16 类按数据条件生成的科研图，导出 PNG/PDF/SVG 和对应 CSV；Markdown/HTML 报告。
- SHA256、原始文本、数组、实际参数、软件版本和运行日志用于追溯。

验证：108 项自动化测试通过；实际测试环境为 Windows x64、Python 3.14.5、RamanSPy 0.2.10。详见 VALIDATION.md。

本发布仅包含工具源代码、配置模板、文档和软件测试，不包含私人原始实验数据、处理结果或本机虚拟环境。自动结果只计算光谱特征；化学归属与机制需要独立证据。报告中的数学拟合和仪器平均次数不能替代独立实验重复。

下载安装：使用本 Release 的完整 Skill ZIP，或克隆仓库后按 README.md 安装。ZIP 附带 SHA256 校验文件与逐文件 MANIFEST.sha256。
