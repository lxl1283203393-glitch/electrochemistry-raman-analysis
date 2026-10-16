# 配置说明

入口读取 YAML 并与程序默认值递归合并；未知字段和不合法参数报错。模板见 `assets/experiment_config.yaml`，该文件是运行参数示例，不是已知的实验事实。除单位、算法及默认开关外，未知实验信息使用 `null`。

| 分组 | 字段与用途 |
| --- | --- |
| `experiment` | `sample_name, preparation_notes, electrolyte, ph, reference_electrode, reference_fill, temperature_c, potential_unit, scan_mode, hold_time_s, laser_wavelength_nm, laser_power_mw, objective_magnification, measurement_position, optical_path_id, ir_corrected` |
| RHE 选择 | `rhe_offset_v` 由有依据的外部校准提供；`allow_approximate_rhe` 默认 false。原始电位仍保留 |
| `filename_overrides` | 按完整文件名覆盖 `potential_v/state/blank_type`，必须有非空 `reason`，不用于无依据自动修正 |
| `blank` | `files` 与 `type` 标识空白；不自动启用扣除 |
| `preprocessing` | `channel, despike, smooth, baseline, normalization, background, sensitivity` |
| `peaks` | prominence 或噪声倍数、最小 SNR、距离/宽度阈值及敏感性比例；距离和宽度使用 cm⁻¹ |
| `fitting` | `windows, shapes, max_components, multistart, seed, criterion` 及约束、稳定性阈值 |
| `plotting` | `dpi, formats, regions, max_annotated_records, max_fit_panels, stack_offset`；`stack_offset: null` 自动选择仅用于显示的偏移 |
| `output` | `formats: [markdown, html]`，`save_arrays` 控制保存处理数组 |
| `qc` | 波数范围、Raw−Dark 容差、尖峰和饱和线索阈值，不自动删除数据 |
| `multivariate` | 默认 `enabled: false`；可选 PCA/NMF、组分数、归一化和种子 |

模板默认拟合窗口 `[200,800]` 仅演示配置格式，不预设物相或固定峰位。应依研究问题和谱图可用范围修改；无须拟合时设 `fitting.windows: []`。拟合数量上限是算法搜索限制，不是物理组分数结论。

`peaks.noise_sigma_floor: null` 由管线根据每条谱的平滑前噪声计算下限；显式数值只能提高该下限。独立调用 detect_peaks 函数时需要自行提供同单位下限，特别是输入已平滑的数组。方法和来源见 preprocessing_methods.md。

示例显式纠正一个经实验记录核实的标签：

```yaml
filename_overrides:
  '1.123cv.csv':
    potential_v: -1.123
    state: applied_potential
    reason: '原始电化学记录确认此文件漏写负号；记录编号由用户填写'
```

不要复制此覆盖到无关实验。`0cv` 默认就是 0 V；若用户证实另一操作含义，单独记录覆盖证据。空白路径及输入路径由本次运行指定，分发包不包含用户机器绝对路径。
