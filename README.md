# Electrochemistry Raman analysis

可执行的 Codex Skill 和 Python 分析包，处理 Optosky Raman Microscope CSV。RamanSPy 负责光谱对象及预处理，NumPy/Pandas/SciPy 补充导入、质量审计、候选峰和数学峰形拟合；Matplotlib 导出 PNG/PDF/SVG 和逐图数据。程序不自动给出物相或 HER 机制结论。

## Windows 安装

GitHub 仓库：[lxl1283203393-glitch/electrochemistry-raman-analysis](https://github.com/lxl1283203393-glitch/electrochemistry-raman-analysis)。完整 Skill 压缩包见 [Releases](https://github.com/lxl1283203393-glitch/electrochemistry-raman-analysis/releases)。下载后可用附带的 `.sha256` 和 `MANIFEST.sha256` 核对完整性。仓库不包含私人原始实验数据或处理报告。

也可直接克隆到 Codex 用户技能目录：

```powershell
git clone https://github.com/lxl1283203393-glitch/electrochemistry-raman-analysis.git "$env:USERPROFILE\.agents\skills\electrochemistry-raman-analysis"
```

需要 Python 3.11 或更新版本。建议使用 Python 3.11–3.14 与独立虚拟环境；具体已测试环境及结果见交付的验证记录。路径中有中文或空格时保留 PowerShell 引号。

将完整项目解压到用户技能目录，例如 `$env:USERPROFILE\.agents\skills\electrochemistry-raman-analysis`。当前 Codex 文档使用 `.agents/skills`；已有 `.codex/skills` 安装也可保留，避免两个目录同时放同名旧版本。重启 Codex 或新建会话使技能被发现。也可从任意项目目录使用命令行工具，无需安装到技能目录。

```powershell
Set-Location "$env:USERPROFILE\.agents\skills\electrochemistry-raman-analysis"
py -3.11 -m venv .venv
& .\.venv\Scripts\python.exe -m pip install --upgrade pip
& .\.venv\Scripts\python.exe -m pip install -e '.[test]'
& .\.venv\Scripts\python.exe -m pytest -q
& .\.venv\Scripts\python.exe scripts/main.py --help
```

若本机只有另一受支持版本，将 `py -3.11` 换为实际解释器。全程直接调用虚拟环境解释器，不依赖 PowerShell 激活脚本。生产运行也可安装 `requirements.txt` 后执行 `scripts/main.py`。安装失败时保留具体依赖错误并修复，不跳过 RamanSPy。

## 四种模式

先将 `assets/experiment_config.yaml` 复制为自己的配置，填写有证据的元数据。下例中 `D:\Raman\batch-a`、`D:\Raman\batch-b` 仅为示例，需替换为实际目录。输出目录必须与原始数据目录分开。

```powershell
# 只导入、审计与质量控制
& .\.venv\Scripts\python.exe scripts/main.py --mode audit --input 'D:\Raman\batch-a' --config '.\experiment.yaml' --output '.\results\audit-a'

# 单批完整处理、候选峰、配置窗口拟合、电位分析及报告
& .\.venv\Scripts\python.exe scripts/main.py --mode analyze --input 'D:\Raman\batch-a' --config '.\experiment.yaml' --output '.\results\analyze-a'

# 多批联合处理与条件核查
& .\.venv\Scripts\python.exe scripts/main.py --mode compare --input 'D:\Raman\batch-a' 'D:\Raman\batch-b' --config '.\experiment.yaml' --output '.\results\compare'

# 读取已有结果重建报告，不重新拟合或处理光谱
& .\.venv\Scripts\python.exe scripts/main.py --mode report --results '.\results\analyze-a' --output '.\results\report-a'
```

安装后可用 `.\.venv\Scripts\electrochemistry-raman.exe` 替代 `python scripts/main.py`。在 Codex 中可说：“使用 `$electrochemistry-raman-analysis` 分析这个目录中的全部原位 Raman CSV，先审计，再生成报告；未知实验条件保持空缺。”首次使用时如 Codex 未找到虚拟环境，提供项目目录及 `.venv\Scripts\python.exe` 位置。

## 工程结构

```text
electrochemistry-raman-analysis/
  SKILL.md                    Codex 调用与证据审查
  agents/openai.yaml          Codex 界面元数据
  pyproject.toml              Python 包和命令行入口
  requirements.txt            运行依赖
  optosky_raman/              可复用的数据模型、导入、审计、处理、峰、比较、绘图与报告模块
  scripts/main.py             可移植命令行入口
  references/                科研规则、格式、算法和配置说明
  assets/experiment_config.yaml
  tests/                     可执行行为与算法测试
```

与最初建议的多个独立脚本相比，公共实现集中在 `optosky_raman` 包中，避免重复解析器、参数漂移和脚本间循环导入。`scripts/main.py` 保留易于直接运行的入口，安装包提供相同 CLI。每个谱图记录同时保留原数组、元数据、SHA256、处理阶段和处理历史。

## 输出和图表

结果目录以 `results.json` 为重建索引，另有 `errors.json`、`config_used.yaml`、`run_log.jsonl`、`audit_report.md` 和 `report.md/report.html`。`provenance/source_snapshots/` 保存原始文本快照，`arrays/` 保存处理数组。`tables/` 中包括 `audit.csv`、`peaks.csv`、`preprocessing_sensitivity.csv`、`detection_sensitivity.csv`、`fit_components.csv`、`fit_models.csv`、`potential_metrics.csv`、`adjacent.csv`、`difference_spectra.csv`、`local_correlations.csv` 和 `batch_comparison.csv`；实际适用表以运行清单为准。`report` 模式读取保存的结果重建 Markdown/HTML，不执行新的科学计算。

`figures/` 默认导出 16 类图的适用项目，`figures/data/` 保存对应长表 CSV：原始全谱、统一处理全谱、电位排序堆叠、时间排序堆叠、分区放大、空白比较、基线方法、候选峰、拟合与残差、峰位/面积/FWHM 对电位、相邻差谱、二维热图、异常对比、跨批次对比。`figures_manifest.json` 记录每类的 `generated/skipped/error`、原因、路径与选取记录；缺少空白、第二批、采集时间或可靠拟合时，相应图不会伪造。

默认 300 dpi PNG 与矢量 PDF/SVG。CSV 包含来源记录、处理方式、波数、原强度及绘图偏移等；标题和图注保留原始电位参比。原强度可为负值，显示偏移不改写分析数据。二维热图按谱图逐行展示，保留重复电位；行间距离不表示电位步长。不同波数轴不会隐式插值成热图。

全谱叠加、堆叠、热图和表保留全部适用记录。基线比较和候选峰图默认各显示前 6 条符合条件的谱图；拟合总览默认前 6 个有有限拟合曲线的结果，每个含残差面板。清单列出选择记录及数量，全部结果仍在候选峰/拟合表中；可调 `plotting.max_annotated_records`、`plotting.max_fit_panels`。避免为每条谱、每个方法与窗口无条件生成数百张图。基线比较不等于参数敏感性全部图形，完整参数变体见敏感性表。

拟合趋势只使用 `status=ok` 的结果，散点不自动连线，也不声称不同电位的同一组件序号对应同一真实峰。未设置独立重复，不输出实验误差线或 p 值；拟合参数的不确定度属于模型内估计，不能替代独立重复。

## 科研边界与扩展

依次审阅 [科学规则](references/scientific_rules.md)、[格式与状态](references/csv_format.md)、[预处理](references/preprocessing_methods.md)、[解释](references/interpretation_guidelines.md) 和 [配置](references/configuration.md)。空白默认仅比较；只有显式核查兼容条件后才扣除。电位文件名缺失负号时保留原标记，只有有原因记录的显式覆盖才更改。

其他仪器可扩展导入适配器，输出相同记录对象；其他样品调整配置窗口、空白类别和文献审阅，无需修改通用峰检测代码；更多批次通过多个 `--input` 目录传入。新光学条件优先补充元数据和对照，不能用归一化掩盖不可比条件。

RamanSPy 来源：[项目](https://github.com/barahona-research-group/RamanSPy)、[官方文档](https://ramanspy.readthedocs.io/)、[论文 DOI](https://doi.org/10.1021/acs.analchem.4c00383)。这些是软件来源，不是任何 Raman 物相归属证据。

实际环境、108项软件测试及科学限制见 [VALIDATION.md](VALIDATION.md)。依赖版本快照见 requirements-tested.txt；该快照对应本次Windows/Python3.14验证，其他环境由pyproject.toml解析兼容版本。私人实验数据未随仓库发布，外部用户不能据此独立复现实验结论。

对保存结果可执行 `python scripts/verify_run.py --results RESULTS --compare OTHER_RUN`，检查原始文件哈希、文本/数组保真、逐图文件及重复计算一致性。`--compare` 可省略。报告与源码压缩包分开保存；源码包不含用户原始CSV及虚拟环境。

## 名称迁移

自 v1.0.1 起，项目与 Skill 名为 `electrochemistry-raman-analysis`，显示名为 Electrochemistry Raman analysis。命令行入口为 `electrochemistry-raman`，保留 `optosky-raman` 命令兼容已有脚本。Python 内部模块仍为 `optosky_raman`，当前原始数据适配范围仍是 Optosky CSV。

旧版本用户请在新虚拟环境安装新项目，避免新旧分发同时管理同一模块；在 Codex 中使用 `$electrochemistry-raman-analysis`。v1.0.0 保留作历史版本。
