"""Four-mode pipeline with source hashing, per-file failures and saved artifacts."""
from __future__ import annotations
from collections import Counter
from datetime import datetime,timezone
from hashlib import sha256
from importlib.metadata import version,PackageNotFoundError
from pathlib import Path
from typing import Any
import json
import platform
import warnings
import numpy as np
import pandas as pd
import yaml
from . import __version__
from .config import load_config
from .importer import load_directory,parse_file
from .audit import audit_records
from .preprocessing import process_record
from .peaks import detect_peaks,noise_sigma
from .fitting import fit_window
from .comparison import compare_records
from .storage import write_json,write_table,serializable
from .visualize import create_figures
from .report import generate_report


def _versions() -> dict:
    packages=['ramanspy','numpy','scipy','pandas','matplotlib','pybaselines','PyYAML','scikit-learn']
    versions={'python':platform.python_version(),'platform':platform.platform(),'electrochemistry-raman-analysis':__version__}
    for name in packages:
        try: versions[name]=version(name)
        except PackageNotFoundError: versions[name]='not_installed'
    code=sha256()
    for module in sorted(Path(__file__).parent.glob('*.py')):
        code.update(module.name.encode('utf-8')); code.update(module.read_bytes())
    versions['implementation_sha256']=code.hexdigest()
    return versions


def _safe_name(record_id: str) -> str:
    return sha256(record_id.encode('utf-8')).hexdigest()[:20]


def _multivariate(records: list,config: dict,output: Path) -> dict:
    import ramanspy as rp
    cfg=config['multivariate']
    if not cfg['enabled']: return {'status':'disabled'}
    usable=[r for r in records if 'unified' in r.stages and r.state!='blank']
    n=int(cfg['n_components'])
    if len(usable)<max(3,n+1):
        return {'status':'skipped','reason':'too few spectra; require at least max(3,n_components+1)'}
    x=usable[0].x[np.argsort(usable[0].x)]
    if any(r.x.shape!=x.shape or not np.allclose(np.sort(r.x),x,rtol=0,atol=1e-6) for r in usable):
        return {'status':'skipped','reason':'multivariate analysis requires common axes; no automatic extrapolation'}
    matrix=np.stack([r.stages['unified'][np.argsort(r.x)] for r in usable])
    normalization=cfg['normalization']
    if normalization=='vector':
        norms=np.linalg.norm(matrix,axis=1,keepdims=True)
        if (norms<=1e-15).any(): return {'status':'skipped','reason':'zero spectrum norm'}
        matrix=matrix/norms
    elif normalization in ('area','auc','max','minmax'):
        normalizers={'area':rp.preprocessing.normalise.AUC,'auc':rp.preprocessing.normalise.AUC,
            'max':rp.preprocessing.normalise.MaxIntensity,'minmax':rp.preprocessing.normalise.MinMax}
        if (np.ptp(matrix,axis=1)<=1e-15).any(): return {'status':'skipped','reason':'normalization undefined for a flat spectrum'}
        matrix=normalizers[normalization]().apply(rp.SpectralContainer(matrix,x)).spectral_data
        if not np.isfinite(matrix).all(): return {'status':'skipped','reason':'normalization produced nonfinite values'}
    elif normalization!='none': raise ValueError('unsupported multivariate normalization')
    method=cfg['method'].lower()
    if method=='nmf' and (matrix<0).any():
        return {'status':'skipped','reason':'NMF requires nonnegative input; negative data not clipped or shifted'}
    container=rp.SpectralContainer(matrix,x)
    model=rp.analysis.decompose.PCA(n_components=n,random_state=cfg['seed']) if method=='pca' else rp.analysis.decompose.NMF(n_components=n,random_state=cfg['seed'],max_iter=1000)
    with warnings.catch_warnings(record=True) as captured:
        warnings.simplefilter('always'); scores,components=model.apply(container)
    scores=np.asarray(scores); components=np.asarray(components)
    write_table(output/'tables'/'multivariate_components.csv',pd.DataFrame(components.T,index=None).assign(raman_shift_cm1=x))
    # RamanSPy returns one score map per component (components x spectra).
    sf=pd.DataFrame(scores.reshape(n,len(usable)).T,columns=[f'component_{i+1}' for i in range(n)])
    sf.insert(0,'record_id',[r.id for r in usable]); write_table(output/'tables'/'multivariate_scores.csv',sf)
    return {'status':'computed','method':method,'n_spectra':len(usable),'n_components':n,
        'normalization':normalization,'warnings':[str(w.message) for w in captured],
        'interpretation':'mathematical components, not chemical species; same-electrode spectra are not independent replicates'}


def run_pipeline(inputs: list[Path],output: Path,config_path: Path|None=None,mode: str='analyze') -> dict:
    if mode not in ('audit','analyze','compare'): raise ValueError(f'pipeline mode must be audit/analyze/compare, got {mode!r}')
    if not inputs: raise ValueError('at least one input directory is required')
    if mode=='analyze' and len(inputs)!=1: raise ValueError('analyze mode requires exactly one input directory; use compare for batches')
    if mode=='compare' and len(inputs)<2: raise ValueError('compare mode requires at least two input directories')
    inputs=[Path(p).resolve() for p in inputs]; output=Path(output).resolve()
    if len(set(inputs))!=len(inputs): raise ValueError('duplicate input directories are not independent batches')
    if any(output==p or output.is_relative_to(p) for p in inputs):
        raise ValueError('output must be outside input directories to protect originals and prevent recursive re-import')
    if len({p.name for p in inputs})!=len(inputs):
        raise ValueError('input batch directory names must be distinct for unambiguous record IDs')
    config=load_config(config_path)
    output.mkdir(parents=True,exist_ok=True)
    (output/'tables').mkdir(exist_ok=True); (output/'arrays').mkdir(exist_ok=True)
    (output/'provenance'/'source_snapshots').mkdir(parents=True,exist_ok=True)
    (output/'config_used.yaml').write_text(yaml.safe_dump(config,allow_unicode=True,sort_keys=False),encoding='utf-8')
    log=output/'run_log.jsonl'
    log.write_text('',encoding='utf-8')
    def event(operation,**fields):
        with log.open('a',encoding='utf-8') as f:
            f.write(json.dumps(serializable({'time':datetime.now(timezone.utc).isoformat(),'operation':operation,**fields}),ensure_ascii=False)+'\n')
    event('start',mode=mode,inputs=[str(p) for p in inputs])
    records=[]; errors=[]
    for directory in inputs:
        loaded,failed=load_directory(directory,config=config)
        records.extend(loaded); errors.extend(failed)
        event('import_batch',input=directory,imported=len(loaded),errors=failed)
    for rec in records:
        # CSV metadata remain source evidence; missing contextual fields may be provided explicitly.
        for key in ('measurement_position','optical_path_id'):
            if rec.metadata.get(key) is None and config['experiment'].get(key) is not None:
                rec.metadata[key]=config['experiment'][key]
    audit=audit_records(records,config)
    write_table(output/'tables'/'audit.csv',audit)
    audit_text=['# 独立数据质量审计','',f'读取 {len(records)} 份光谱，导入错误 {len(errors)} 项。',
        '', '[完整逐文件审计表](tables/audit.csv)。内部采集时间用于排序，文件名时间不代替采集时间。',
        '', 'Raw−Dark 与 Dark Subtracted 的差异逐文件保存；保持仪器输出，不自行替换。',
        '', '异常标记是可疑采集状态提示，不能单独确定气泡、焦点或化学变化。原文件均保留。',
        '', '负值、零值、疑似尖峰/平台、相似度和采集设置见表。饱和在无ADC上限时只报告迹象。',
        '', '仪器平均次数不等于独立重复；不存在自动显著性检验。',
        '', '导入错误：', '```json',json.dumps(errors,ensure_ascii=False,indent=2),'```']
    (output/'audit_report.md').write_text('\n'.join(audit_text),encoding='utf-8')
    event('audit',rows=len(audit))
    peak_rows=[]; fits=[]; sensitivity=[]; detection_sensitivity=[]; processed=0; background_provenance=None
    if mode!='audit':
        bg=config['preprocessing']['background']; blank=None
        if bg['enabled']:
            if not bg['file']: raise ValueError('background.file required when subtraction enabled')
            blank_path=Path(bg['file'])
            if not blank_path.is_absolute():
                base=Path(config_path).resolve().parent if config_path else inputs[0]
                blank_path=base/blank_path
            blank=parse_file(blank_path,config=config)
            key=_safe_name('background/'+str(blank.path))
            np.savez_compressed(output/'arrays'/f'{key}.npz',x=blank.x,**blank.channels)
            (output/'provenance'/'source_snapshots'/f'{key}.txt').write_text(blank.metadata['original_text'],encoding='utf-8',newline='')
            background_provenance={'path':str(blank.path),'sha256':blank.sha256,'arrays':f'arrays/{key}.npz',
                'snapshot':f'provenance/source_snapshots/{key}.txt','parameters':bg}
        for rec in records:
            try:
                details=process_record(rec,config['preprocessing'],blank)
                sensitivity.extend({'record_id':rec.id,**d} for d in details)
                processed+=1
                order=np.argsort(rec.x); x=rec.x[order]; y=rec.stages['unified'][order]
                # Smoothing creates correlated samples; a second-difference MAD of C
                # alone underestimates noise. Carry forward the pre-smoothing noise.
                source=rec.channels[config['preprocessing']['channel']][order]
                raw_noise=noise_sigma(x,source)
                corrected=rec.stages['smoothed'][order]-rec.stages['baseline'][order]
                gain=float(np.ptp(y)/np.ptp(corrected)) if np.ptp(corrected)>1e-15 else 1.0
                configured_floor=config['peaks'].get('noise_sigma_floor') or 0
                effective_peak_config={**config['peaks'],'noise_sigma_floor':max(configured_floor,raw_noise*gain)}
                raw_config={**config['peaks'],'noise_sigma_floor':max(configured_floor,raw_noise)}
                original=detect_peaks(x,rec.channels['dark_subtracted'][order],raw_config)
                detected=detect_peaks(x,y,effective_peak_config)
                rec.history.append({'operation':'candidate_detection','parameters':effective_peak_config,
                    'noise_source':'pre-smoothing selected channel; floor prevents correlated-noise underestimation',
                    'pre_smoothing_noise_sigma':raw_noise,'normalization_gain':gain})
                for peak in detected:
                    tolerance=max(5.0,float(peak['fwhm_cm1'])*.25)
                    raw_support=any(abs(p['position_cm1']-peak['position_cm1'])<=tolerance for p in original)
                    peak_rows.append({'record_id':rec.id,'batch':rec.batch,'filename':rec.path.name,
                        'potential_v':rec.potential_v,'state':rec.state,'processing':'unified',
                        'raw_supported':raw_support,'evidence_status':'candidate_with_raw_support' if raw_support else 'processed_only_unconfirmed',**peak})
                for factor in config['peaks']['sensitivity_factors']:
                    pc={**effective_peak_config,'prominence_sigma':config['peaks']['prominence_sigma']*factor}
                    if pc['prominence'] is not None: pc['prominence']*=factor
                    varied=detect_peaks(x,y,pc)
                    for p in varied: detection_sensitivity.append({'record_id':rec.id,'variant':f'detection_{factor:g}',
                        'factor':factor,'position_cm1':p['position_cm1'],'prominence':p['prominence'],'snr':p['snr']})
                # Re-detect baseline/smoothing variants so sensitivity is about bands, not only totals.
                for detail in details:
                    variant=detail['variant']
                    for p in detect_peaks(x,rec.stages[variant][order],raw_config):
                        detection_sensitivity.append({'record_id':rec.id,'variant':variant,'factor':None,
                            'position_cm1':p['position_cm1'],'prominence':p['prominence'],'snr':p['snr']})
                fitting_options={k:v for k,v in config['fitting'].items() if k!='windows'}
                for window in config['fitting']['windows']:
                    fit=fit_window(x,y,{**fitting_options,'range':window})
                    fit.update(record_id=rec.id,window=window,processing='unified')
                    fits.append(fit)
                event('processed',record_id=rec.id,candidate_peaks=len(detected),warnings=rec.warnings)
            except Exception as exc:
                failure={'file':str(rec.path),'record_id':rec.id,'stage':'analysis','error_type':type(exc).__name__,'error':str(exc)}
                errors.append(failure); event('error',**failure)
    peak_frame=pd.DataFrame(peak_rows)
    if not peak_frame.empty:
        varied=pd.DataFrame(detection_sensitivity)
        for index,row in peak_frame.iterrows():
            candidates=varied[(varied.record_id==row.record_id)&(~varied.variant.str.startswith('detection_'))] if not varied.empty else pd.DataFrame()
            tolerance=max(5.0,float(row.fwhm_cm1)*.25)
            matched=candidates[np.abs(candidates.position_cm1-row.position_cm1)<=tolerance] if not candidates.empty else pd.DataFrame()
            variants=sorted(matched.variant.unique()) if not matched.empty else []
            methods=sorted({v.split('_')[0] for v in variants if not v.startswith('smooth_')})
            peak_frame.loc[index,'supporting_preprocessing_variants']=';'.join(variants)
            peak_frame.loc[index,'supporting_baseline_methods']=';'.join(methods)
            peak_frame.loc[index,'preprocessing_variant_count']=len(variants)
            peak_frame.loc[index,'evidence_status']='raw_and_processed_candidate' if row.raw_supported else 'processed_only_unconfirmed'
    write_table(output/'tables'/'peaks.csv',peak_frame)
    write_table(output/'tables'/'preprocessing_sensitivity.csv',pd.DataFrame(sensitivity))
    write_table(output/'tables'/'detection_sensitivity.csv',pd.DataFrame(detection_sensitivity))
    component_rows=[]; model_rows=[]
    for fit in fits:
        for comp in fit['components']:
            component_rows.append({'record_id':fit['record_id'],'window':str(fit['window']),'status':fit['status'],
                **{k:v for k,v in comp.items() if k!='fitted'}})
        for model in fit['models']:
            model_rows.append({'record_id':fit['record_id'],'window':str(fit['window']),**model})
    write_table(output/'tables'/'fit_components.csv',pd.DataFrame(component_rows))
    write_table(output/'tables'/'fit_models.csv',pd.DataFrame(model_rows))
    comparisons=compare_records(records,fits,config)
    for name,frame in comparisons.items():
        if isinstance(frame,pd.DataFrame): write_table(output/'tables'/f'{name}.csv',frame)
    figures=create_figures(records,audit,peak_frame,fits,comparisons,output,config) if mode!='audit' else []
    write_json(output/'figures_manifest.json',figures)
    for figure in figures:
        if figure['status']=='error': errors.append({'stage':'figure','category':figure['category'],'error':figure['reason']})
    try: multivariate=_multivariate(records,config,output) if mode!='audit' else {'status':'not_requested_in_audit'}
    except Exception as exc:
        multivariate={'status':'failed','error':str(exc)}
        errors.append({'stage':'multivariate','error_type':type(exc).__name__,'error':str(exc)})
    summaries=[]
    for rec in records:
        key=_safe_name(rec.id); snapshot=output/'provenance'/'source_snapshots'/f'{key}.txt'
        snapshot.write_text(rec.metadata.get('original_text',''),encoding='utf-8',newline='')
        metadata={k:v for k,v in rec.metadata.items() if k not in ('original_text','original_rows')}
        arr_path=output/'arrays'/f'{key}.npz'
        # Saved arrays are never a substitute for original file hashes and preserved raw text.
        if config['output']['save_arrays']:
            np.savez_compressed(arr_path,x=rec.x,**{f'channel_{k}':v for k,v in rec.channels.items()},
                **{f'stage_{k}':v for k,v in rec.stages.items()})
        unchanged=sha256(rec.path.read_bytes()).hexdigest()==rec.sha256
        if not unchanged: errors.append({'file':str(rec.path),'stage':'source_integrity','error':'Source file changed during run'})
        summaries.append({'record_id':rec.id,'path':str(rec.path),'batch':rec.batch,'filename':rec.path.name,
            'sha256':rec.sha256,'source_unchanged':unchanged,'encoding':rec.encoding,
            'state':rec.state,'potential_v':rec.potential_v,'acquired_at':rec.acquired_at,
            'metadata':metadata,'warnings':rec.warnings,'history':rec.history,
            'snapshot':snapshot.relative_to(output).as_posix(),
            'arrays':arr_path.relative_to(output).as_posix() if config['output']['save_arrays'] else None})
    if background_provenance:
        unchanged=sha256(blank.path.read_bytes()).hexdigest()==blank.sha256
        background_provenance['source_unchanged']=unchanged
        if not unchanged: errors.append({'file':str(blank.path),'stage':'source_integrity','error':'Background source changed during run'})
    result={'schema_version':'1.0','created_at':datetime.now(timezone.utc).isoformat(),
        'mode':mode,'status':'partial_failure' if errors else 'complete','config':config,
        'inputs':[str(p) for p in inputs],'counts':{'imported':len(records),'processed':processed,'failed':len(errors)},
        'versions':_versions(),'random_seed':config['fitting']['seed'],
        'random_seeds':{'fitting':config['fitting']['seed'],'multivariate':config['multivariate']['seed']},'records':summaries,
        'fits':fits,'figures':figures,'errors':errors,'multivariate':multivariate,'background_source':background_provenance,
        'comparisons':{'confounding':comparisons['confounding'],'batch_pairs':len(comparisons['batch_comparison'])}}
    write_json(output/'errors.json',errors); write_json(output/'results.json',result)
    generate_report(output)
    event('complete',status=result['status'],counts=result['counts'])
    return result
