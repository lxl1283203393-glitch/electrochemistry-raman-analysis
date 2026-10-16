"""Reports from saved results only; no refitting or chemical assignment."""
from __future__ import annotations
from collections import Counter
from html import escape
import json
from pathlib import Path
import pandas as pd
import markdown


def _table(frame: pd.DataFrame,columns: list[str]|None=None,n: int=12) -> str:
    if frame.empty: return '无可用记录；见下方数据缺口与跳过原因。'
    if columns: frame=frame[[c for c in columns if c in frame.columns]]
    frame=frame.head(n).fillna('null')
    cols=list(frame.columns)
    text=['| '+' | '.join(cols)+' |','| '+' | '.join(['---']*len(cols))+' |']
    for row in frame.itertuples(index=False,name=None):
        text.append('| '+' | '.join(str(v).replace('|','/').replace('\n',' ') for v in row)+' |')
    return '\n'.join(text)


def generate_report(results_dir: Path,output_dir: Path|None=None) -> dict:
    results_dir=Path(results_dir).resolve(); output_dir=Path(output_dir or results_dir).resolve()
    saved=results_dir/'results.json'
    if not saved.exists(): raise FileNotFoundError(f'saved results not found: {saved}')
    data=json.loads(saved.read_text(encoding='utf-8'))
    output_dir.mkdir(parents=True,exist_ok=True)
    # If report is emitted elsewhere, retain relative artifact links to the saved analysis.
    import os
    try: prefix=Path(os.path.relpath(results_dir,output_dir)).as_posix()
    except ValueError: prefix=results_dir.as_uri()  # Windows cross-volume report destinations.
    def link(path): return path if prefix=='.' else prefix+'/'+path
    def read(name):
        path=results_dir/'tables'/f'{name}.csv'
        return pd.read_csv(path) if path.exists() else pd.DataFrame()
    audit,peaks,components,adjacent,batch=[read(n) for n in ('audit','peaks','fit_components','adjacent','batch_comparison')]
    counts=data['counts']; exp=data['config']['experiment']
    missing=[k for k,v in exp.items() if v is None]
    mode=data['mode']; figure_rows=data.get('figures',[])
    scale_rows=adjacent[adjacent['common_scale_suspected'].astype(str).str.lower()=='true'] if 'common_scale_suspected' in adjacent else pd.DataFrame()
    if scale_rows.empty:
        scale_interpretation='本次没有相邻谱满足共同倍率变化筛查判据。不能仅凭整体下降确定气泡/焦点影响，也不能排除这些竞争解释；请同时检查谱形残差。'
    else:
        scale_interpretation=f'二级：支持性解释——{len(scale_rows)} 对相邻谱满足共同倍率筛查判据（见表中 common_scale_suspected）。这些具体比较与共同光学倍率变化相容；气泡、焦点和背景改变仍须独立记录确认。'
    flagged_change=('candidate_change' in adjacent and adjacent['candidate_change'].astype(str).str.lower().eq('true').any())
    if counts['processed']<2:
        cautious_statement='“本次可用处理光谱不足两份，无法形成电位依赖或谱形演变结论；需要补充数据后再分析。”'
    elif flagged_change:
        cautious_statement='“本批次部分相邻光谱达到描述性变化筛查判据；电位、时间和采集条件的贡献尚未分离，目前不足以确认特定工作态物相或重构机制。”'
    else:
        cautious_statement='“本次获得了可追溯光谱与数学特征；相邻谱未触发配置中的变化筛查判据，这不构成结构不变或不存在弱谱带的证明。”'
    lines=['# Optosky 电化学原位 Raman 科研报告','',
        f"运行模式：`{mode}`；状态：`{data['status']}`；生成数据时间：{data['created_at']}。",'',
        '**本报告的数值来自程序计算；化学归属、结构演变和因果机制没有经过自动确认。**','',
        '## 1. 实验信息及缺失参数','',
        '实验信息来自配置与文件元数据。配置中的 null 保留为未知，仪器平均次数不视为独立重复。','',
        _table(pd.DataFrame([{'parameter':k,'value':str(v)} for k,v in exp.items()]),n=100),'',
        '缺失参数：'+(', '.join(missing) if missing else '未发现 null 配置字段，但仍需核验采集记录。')+'。','',
        '## 2. 数据完整性审计','',
        f"导入 {counts['imported']} 份，统一处理 {counts['processed']} 份，错误事件 {counts['failed']} 项。",'',
        f"[完整审计表]({link('tables/audit.csv')})；[独立审计报告]({link('audit_report.md')})；[错误明细]({link('errors.json')})。",'',
        _table(audit,n=8),'',
        '文件名日期不替代内部采集时间；电位覆盖规则和原始文件哈希保存在 results.json 与 source_snapshots。','',
        '## 3. 预处理方法和参数','',
        'A：仪器 Dark Subtracted；B：仪器 BaseLine Subtracted；C：RamanSPy 统一预处理。A/B 原数组保持不变。',
        '负值保留；基线敏感性分支不归一化，C 的归一化设置单独记录。SavGol 窗口按数据点计，在非均匀波数轴上对应的物理宽度会变化。','',
        '```json',json.dumps(data['config']['preprocessing'],ensure_ascii=False,indent=2),'```','',
        f"[完整配置]({link('config_used.yaml')})；[方法/参数敏感性]({link('tables/preprocessing_sensitivity.csv')})；软件版本见 results.json。",'',
        '## 4. 光谱总体特征','',
        '一级：直接观察仅限原始/处理数据中可复核的光谱特征。整体积分和谱形相似度属于描述性计算值。',
        f"[全谱数值表]({link('tables/potential_metrics.csv')})。积分、峰高受背景、光路及采集设置影响。",'',
        _table(read('potential_metrics'),['filename','batch','potential_v','integrated_intensity','raw_integrated_intensity'],n=12),'',
        '## 5. 候选谱峰与数学拟合','',
        f"候选峰 {len(peaks)} 条。峰拟合状态：{dict(Counter(f.get('status') for f in data.get('fits',[])))}。",'',
        f"[候选峰及原始证据核对]({link('tables/peaks.csv')})；[检测敏感性]({link('tables/detection_sensitivity.csv')})；[模型比较]({link('tables/fit_models.csv')})；[分量参数]({link('tables/fit_components.csv')})。",'',
        _table(components,['record_id','window','status','center_cm1','fwhm_cm1','height','area','shape'],n=12),'',
        '分量是数学峰形分解，面积为该峰形在全实数轴的解析积分，不代表物相含量。局部协方差给出的条件参数不确定度不代表实验重复误差；不可辨识拟合不用于趋势图。',
        '当残差明显相关时，局部协方差不确定度将抑制；status=ok 仍只是条件数学拟合，不保证完整谱形模型或化学归属成立。',
        '检测到宽带不能直接等同于识别出多个独立峰。跨电位的不同分量不自动视为同一化学振动。','',
        '## 6. 电位与时间响应','',
        json.dumps(data.get('comparisons',{}).get('confounding',{}),ensure_ascii=False,indent=2),'',
        f"[相邻采集比较]({link('tables/adjacent.csv')})；[局部谱形相关性]({link('tables/local_correlations.csv')})；差谱见 figures。",'',
        '电位响应与时间漂移可能混杂。候选突变只定位到采样间隔，扫描步长不是结构转变电位精度。',
        '同一标称电位的多次采集单列如下，不当作独立重复：','',
        _table(read('same_potential_comparison'),['record_id','other_id','potential_v','raw_correlation','raw_shape_residual','correlation','shape_residual'],n=15),'',
        '标称电位始终保留原参比与符号；RHE 缺必要参数时保留 null。0 V、无外加电位和开路状态分别记录。','',
        '## 7. 异常采集与共同光学倍率','',
        _table(adjacent,['record_id','previous_id','scale_ratio','shape_residual','raw_scale_ratio','common_scale_suspected','candidate_change'],n=20),'',
        scale_interpretation+' 异常文件保留。','',
        '## 8. 空白对照','',
        f"识别空白 {sum(r['state']=='blank' for r in data.get('records',[]))} 份。空白首先比较，不默认扣除。",'',
        '仅显式确认条件匹配、波数轴匹配且设置不冲突时执行扣除；比例因子、扣除前后数组和负值诊断写入处理历史。',
        '无样品空白不能代替裸 Ni foam、电解液及光学窗口各自的对照；缺哪一种就保留对应缺口。','',
        '## 9. 跨批次比较','',_table(batch,n=15),'',
        '不同激光功率、积分时间、采样位置或光学条件下不作未经论证的绝对强度比较。同电极连续采集与仪器平均次数都不视为独立样本重复。','',
        '## 10. 峰归属候选与证据边界','',
        '四级：无法判断——程序不建立化学峰归属库，本次没有自动核验的归属文献。后续人工候选应逐条注明 DOI、原始参比谱、测试环境和竞争归属。',
        '本报告未由单峰证明 Zn 溶出、空位、CoP/Ni₂P 异质界面、电子转移、d 带中心变化或 HER 活性位。','',
        '## 11. 可能的结构演变解释','',
        '三级：待验证假设——若同一光谱特征在原始数据、合理预处理及独立实验中重复出现，可进一步检验其是否与表面结构或吸附状态变化有关。本次程序结果不确认重构路径。','',
        '## 12. 竞争解释','',
        '液态水/界面水、Ni foam、窗口、荧光背景、HER 气泡、光学焦点及激光作用均可能影响信号；首先对照原始全谱和共同倍率诊断。','',
        '## 13. 局限性与执行状态','',
        '未提供独立样品重复，因此只报告描述性结果，无 p 值、重复误差线或文献 benchmark。谱仪分辨率和校准信息缺失时，不把拟合小数位视为真实精度。',
        f"错误明细：{json.dumps(data.get('errors',[]),ensure_ascii=False)}",'',
        '图像覆盖：','',_table(pd.DataFrame(figure_rows),['category','status','reason'],n=100),'',
        '## 14. 最小补充实验建议','',
        '补充同光路的 1 M KOH、裸 Ni foam 和窗口空白；在可疑电位前后返回同一电位复测以分离时间漂移；记录气泡/焦点状态，并用独立电极重复确认关键谱形变化。',
        '如需 RHE 电位，补齐参比内充液、温度、pH 及 iR 处理记录。若要支持物相/溶出/机制，需另加对应独立表征。','',
        '## 15. 可用于组会的谨慎表述','',
        cautious_statement,'',
        '## 证据索引与图像','',
        f"[运行记录]({link('run_log.jsonl')})；[源文件和软件版本索引]({link('results.json')})；[图表清单]({link('figures_manifest.json')})。",'']
    for fig in figure_rows:
        if fig.get('status')=='generated' and fig.get('path'):
            lines.extend([f"### {fig['title']}",'',f"![{fig['title']}]({link(fig['path'])})",'',
                f"[绘图数据]({link(fig['data_path'])})；处理：{fig.get('processing','见图表清单')}。",''])
    text='\n'.join(lines)
    formats=data['config'].get('output',{}).get('formats',['markdown','html'])
    outputs={}
    if 'markdown' in formats:
        (output_dir/'report.md').write_text(text,encoding='utf-8'); outputs['markdown']=str(output_dir/'report.md')
    if 'html' in formats:
        body=markdown.markdown(text,extensions=['tables','fenced_code'])
        html='<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>Optosky Raman 科研报告</title><style>body{font-family:system-ui,"Microsoft YaHei",sans-serif;max-width:1160px;margin:32px auto;padding:0 24px;line-height:1.7;color:#17212d}img{max-width:100%}table{border-collapse:collapse;display:block;overflow:auto;font-size:13px}td,th{padding:6px 10px;border:1px solid #ccd3db}pre{overflow:auto;background:#f1f4f7;padding:16px}h2{margin-top:2em;border-bottom:1px solid #ccd3db}a{color:#145b96}</style>'+body+'</html>'
        (output_dir/'report.html').write_text(html,encoding='utf-8'); outputs['html']=str(output_dir/'report.html')
    return outputs
