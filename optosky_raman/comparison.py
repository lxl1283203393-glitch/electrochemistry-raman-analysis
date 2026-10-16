"""Descriptive potential/time and batch comparisons, never independent repeats."""
from __future__ import annotations
from itertools import combinations
from typing import Any
import numpy as np
import pandas as pd


def _corr(a: np.ndarray,b: np.ndarray) -> float|None:
    if np.std(a)<1e-12 or np.std(b)<1e-12: return None
    return float(np.corrcoef(a,b)[0,1])


def _pair(a: Any,b: Any,channel: str='unified') -> dict:
    # Comparison interpolation only on the shared measured interval, never extrapolation.
    lo,hi=max(a.x.min(),b.x.min()),min(a.x.max(),b.x.max())
    ai=np.argsort(a.x); bi=np.argsort(b.x)
    ax=a.x[ai]; bx=b.x[bi]
    mask=(ax>=lo)&(ax<=hi)
    x=ax[mask]
    if len(x)<5:
        return {'status':'insufficient_overlap','scale_ratio':None,'shape_residual':None,
                'correlation':None,'common_scale_suspected':False,'axis_interpolation':None}
    ya=a.stages.get(channel,a.channels['dark_subtracted'])[ai][mask]
    yb=np.interp(x,bx,b.stages.get(channel,b.channels['dark_subtracted'])[bi])
    denom=float(np.dot(ya,ya))
    scale=float(np.dot(ya,yb)/denom) if denom>1e-20 else None
    norm=float(np.linalg.norm(yb))
    residual=float(np.linalg.norm(yb-scale*ya)/norm) if scale is not None and norm>1e-20 else None
    return {'status':'ok','x':x,'a':ya,'b':yb,'scale_ratio':scale,
        'shape_residual':residual,'correlation':_corr(ya,yb),
        'common_scale_suspected':bool(scale is not None and residual is not None and scale<.5 and residual<.15),
        'axis_interpolation':not (a.x.shape==b.x.shape and np.allclose(a.x,b.x,atol=1e-6,rtol=0))}


def compare_records(records: list, fits: list[dict], config: dict) -> dict:
    rows=[]; adjacent=[]; differences=[]; regional=[]; batches=[]; repeated=[]; confounding={}
    for rec in records:
        if 'unified' not in rec.stages: continue
        y=rec.stages['unified']; raw=rec.channels['dark_subtracted']; order=np.argsort(rec.x)
        rows.append({'record_id':rec.id,'batch':rec.batch,'filename':rec.path.name,
            'state':rec.state,'potential_v':rec.potential_v,
            'time':rec.acquired_at.isoformat() if rec.acquired_at else None,
            'integrated_intensity':float(np.trapezoid(y[order],rec.x[order])),
            'raw_integrated_intensity':float(np.trapezoid(raw[order],rec.x[order])),
            'maximum_intensity':float(y.max()),'processing':'unified',
            'integration_time_s':rec.metadata.get('integration_time_s'),
            'laser_power_mw':rec.metadata.get('laser_power_mw')})
    available=[r for r in records if 'unified' in r.stages]
    for batch in sorted(set(r.batch for r in available)):
        applied=[r for r in available if r.batch==batch and r.state=='applied_potential' and r.potential_v is not None]
        timed=sorted([r for r in applied if r.acquired_at],key=lambda r:(r.acquired_at,r.id))
        if len(timed)>=3:
            t=np.array([(r.acquired_at-timed[0].acquired_at).total_seconds() for r in timed]); v=np.array([r.potential_v for r in timed])
            corr=_corr(t,v)
            rank_corr=_corr(pd.Series(t).rank().to_numpy(),pd.Series(v).rank().to_numpy())
            confounding[batch]={'n':len(timed),'pearson_time_potential':corr,'rank_time_potential':rank_corr,
                'time_potential_confounded':bool(rank_corr is not None and abs(rank_corr)>.8),
                'interpretation':'sequential observations on one electrode; causal potential response is not isolated'}
        else:
            confounding[batch]={'n':len(timed),'time_potential_confounded':None,'interpretation':'too few timestamped applied-potential observations'}
        for a,b in zip(timed,timed[1:]):
            pair=_pair(a,b); rawpair=_pair(a,b,'dark_subtracted')
            row={'record_id':b.id,'previous_id':a.id,'batch':batch,
                'potential_v':b.potential_v,'previous_potential_v':a.potential_v,
                'interval_s':(b.acquired_at-a.acquired_at).total_seconds(),
                **{k:v for k,v in pair.items() if k not in ('x','a','b')}}
            row['raw_scale_ratio']=rawpair.get('scale_ratio')
            row['raw_shape_residual']=rawpair.get('shape_residual')
            setting_keys=('integration_time_s','laser_power_mw','averages')
            changed=[k for k in setting_keys if a.metadata.get(k)!=b.metadata.get(k)]
            known=all(a.metadata.get(k) is not None and b.metadata.get(k) is not None for k in ('integration_time_s','laser_power_mw'))
            row['acquisition_settings_changed']=';'.join(changed)
            row['absolute_intensity_comparable']=not changed and known
            row['common_scale_suspected']=bool(row.get('common_scale_suspected',False) and not changed and known)
            row['candidate_change']=bool(not changed and pair.get('correlation') is not None and
                (pair['correlation']<.8 or pair.get('scale_ratio',1)<.7 or pair.get('scale_ratio',1)>1.4))
            row['screening_rule']='correlation<0.8 or scale<0.7 or scale>1.4; settings changes excluded; not a transition determination'
            adjacent.append(row)
            if pair['status']=='ok' and not changed and known:
                differences.extend({'record_id':b.id,'previous_id':a.id,'batch':batch,
                    'raman_shift_cm1':float(x),'difference':float(d),'potential_v':b.potential_v}
                    for x,d in zip(pair['x'],pair['b']-pair['a']))
                for region in config.get('plotting',{}).get('regions',[[200,800],[2800,3800]]):
                    mask=(pair['x']>=region[0])&(pair['x']<=region[1])
                    regional.append({'record_id':b.id,'previous_id':a.id,'batch':batch,
                        'range_low':region[0],'range_high':region[1],
                        'correlation':_corr(pair['a'][mask],pair['b'][mask]) if mask.sum()>=5 else None})
    for a,b in combinations(available,2):
        if a.state!='applied_potential' or b.state!='applied_potential': continue
        if a.potential_v is None or b.potential_v is None or not np.isclose(a.potential_v,b.potential_v,atol=1e-6,rtol=0): continue
        pair=_pair(a,b)
        if a.batch==b.batch:
            rawpair=_pair(a,b,'dark_subtracted')
            repeated.append({'record_id':a.id,'other_id':b.id,'batch':a.batch,'potential_v':a.potential_v,
                'time':str(a.acquired_at),'other_time':str(b.acquired_at),
                'correlation':pair.get('correlation'),'scale_ratio':pair.get('scale_ratio'),
                'shape_residual':pair.get('shape_residual'),'raw_correlation':rawpair.get('correlation'),
                'raw_scale_ratio':rawpair.get('scale_ratio'),'raw_shape_residual':rawpair.get('shape_residual'),
                'interpretation':'same nominal potential; acquisition condition labels are not assumed equivalent; not independent replicates'})
            continue
        fields=['integration_time_s','laser_power_mw','measurement_position','optical_path_id']
        mismatched=[key for key in fields if a.metadata.get(key)!=b.metadata.get(key)]
        missing=[key for key in fields if a.metadata.get(key) is None or b.metadata.get(key) is None]
        batches.append({'record_id':a.id,'other_id':b.id,'batch':a.batch,'other_batch':b.batch,
            'potential_v':a.potential_v,'correlation':pair.get('correlation'),
            'shape_residual':pair.get('shape_residual'),
            'absolute_intensity_comparable':not mismatched and not missing,
            'setting_mismatches':';'.join(mismatched),'missing_conditions':';'.join(missing),
            'interpretation':'descriptive spectrum-shape comparison; no independent replicate assumption'})
    fit_rows=[]
    lookup={r.id:r for r in available}
    for fit in fits:
        if fit.get('status')!='ok' or fit.get('record_id') not in lookup: continue
        rec=lookup[fit['record_id']]
        for i,c in enumerate(fit.get('components',[])):
            fit_rows.append({'record_id':rec.id,'batch':rec.batch,'filename':rec.path.name,
                'potential_v':rec.potential_v,'state':rec.state,'window':str(fit.get('window')),
                'component_index':i,'position_cm1':c['center_cm1'],'height':c['height'],
                'area':c['area'],'fwhm_cm1':c['fwhm_cm1'],'shape':c['shape'],
                'processing':'unified','tracking':'unlinked mathematical components; identity across potentials not established'})
    return {'potential_metrics':pd.DataFrame(rows),'adjacent':pd.DataFrame(adjacent),
        'difference_spectra':pd.DataFrame(differences),'local_correlations':pd.DataFrame(regional),
        'batch_comparison':pd.DataFrame(batches),'fit_metrics':pd.DataFrame(fit_rows),
        'same_potential_comparison':pd.DataFrame(repeated),
        'confounding':confounding}
