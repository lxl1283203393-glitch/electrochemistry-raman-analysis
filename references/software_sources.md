# 已核查的软件来源

核查日期：2026-10-16。本项目调用官方库 API，不把软件参考资料当作化学峰归属文献。

- [Codex 官方技能规范](https://learn.chatgpt.com/docs/build-skills)：技能目录、name/description、渐进加载、用户级 `.agents/skills` 与自动/显式调用。已结合本机 skill-creator 规范验证。
- [RamanSPy 官方仓库](https://github.com/barahona-research-group/RamanSPy)：实际安装测试版本 0.2.10；使用 Spectrum、SpectralContainer、Pipeline、WhitakerHayes、SavGol、ASLS、ARPLS、AIRPLS 和分析接口。
- [RamanSPy 官方 API](https://ramanspy.readthedocs.io/en/latest/)：预处理、数据容器、PCA/NMF。具体签名同时通过安装版本运行检查；网页文档与本机版本可能不同。
- [软件论文](https://doi.org/10.1021/acs.analchem.4c00383)：Georgiev 等，Analytical Chemistry (2024), 96, 8492–8500。此 DOI 已由官方仓库/文档核验，不提供任何样品物相归属。
- [Gamry 参比电位换算说明](https://www.gamry.com/resources-2/electrochemical-calculators-tools/convert-potentials-to-another-reference-electrode/)：参比内充液不同需区别处理；代码中的常见 Ag/AgCl 标称偏移只用于注明条件的换算。优先使用实际电极校准的 rhe_offset_v。无填充液/pH/温度时默认不换算。

完整依赖和版本见 requirements-tested.txt 及每次运行的 results.json。安装依赖失败时不得悄悄换成未使用 RamanSPy 的核心路径。
