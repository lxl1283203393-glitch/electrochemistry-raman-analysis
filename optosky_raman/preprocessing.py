"""RamanSPy preprocessing; originals and every numerical stage are retained."""
from __future__ import annotations
from datetime import datetime, timezone
from importlib.metadata import version
from typing import Any
from copy import copy, deepcopy
import warnings
import numpy as np
import ramanspy as rp


def subtract_background(record: Any, blank: Any, options: dict) -> tuple[np.ndarray, dict]:
    if not options.get('conditions_verified',False):
        raise ValueError('background conditions must be explicitly verified before subtraction')
    if record.x.shape != blank.x.shape or not np.allclose(record.x,blank.x,atol=1e-6,rtol=0):
        raise ValueError('background axis mismatch; automatic resampling is not permitted')
    for key in ('integration_time_s','laser_power_mw','optical_path_id','measurement_position',
                'objective_magnification','laser_wavelength_nm','averages'):
        a,b=record.metadata.get(key),blank.metadata.get(key)
        if a is not None and b is not None and a!=b:
            raise ValueError(f'background conditions mismatch: {key} ({a!r} vs {b!r})')
    scale=float(options.get('scale',1))
    if not np.isfinite(scale) or scale<0:
        raise ValueError('background scale must be finite and nonnegative')
    y=record.channels['dark_subtracted']-scale*blank.channels['dark_subtracted']
    details={'scale':scale,'blank_file':str(blank.path),'conditions_verified':True,
        'blank_sha256':getattr(blank,'sha256',None),
        'negative_fraction':float(np.mean(y<0)),
        'negative_min':float(y.min()),'interpretation':'negative dips may be subtraction artifacts; no chemical attribution'}
    return y,details


def _baseline(spectrum: Any, settings: dict) -> Any:
    method=str(settings.get('method','asls')).lower()
    if method=='none':
        return spectrum
    classes={'asls':rp.preprocessing.baseline.ASLS,
             'arpls':rp.preprocessing.baseline.ARPLS,
             'airpls':rp.preprocessing.baseline.AIRPLS}
    if method not in classes:
        raise ValueError(f'unsupported baseline method: {method}')
    kwargs={k:settings[k] for k in ('lam','max_iter','tol') if k in settings}
    if method=='asls': kwargs['p']=settings.get('p',.01)
    return rp.preprocessing.Pipeline([classes[method](**kwargs)]).apply(spectrum)


def process_record(record: Any, options: dict, blank: Any|None=None) -> list[dict]:
    """Publish stages only after the complete preprocessing operation succeeds.

    A failed file retains its earlier state; a fresh import therefore cannot expose
    an incomplete C spectrum to downstream comparisons or multivariate analysis.
    """
    pending=copy(record)
    pending.stages={}
    pending.history=deepcopy(record.history)
    pending.warnings=list(record.warnings)
    sensitivity=_process_record_inplace(pending,deepcopy(options),blank)
    record.stages=pending.stages
    record.history=pending.history
    record.warnings=pending.warnings
    return sensitivity


def _process_record_inplace(record: Any, options: dict, blank: Any|None=None) -> list[dict]:
    """Compute temporary A/B/C arrays; public process_record commits atomically."""
    order=np.argsort(record.x)
    inverse=np.argsort(order)
    x=record.x[order]
    for key in ('dark_subtracted','baseline_subtracted'):
        if key in record.channels: record.stages[key]=record.channels[key].copy()
    channel=options.get('channel','dark_subtracted')
    if channel not in record.channels:
        raise ValueError(f'{record.path}: requested preprocessing channel {channel!r} unavailable')
    y=record.channels[channel].copy()
    bg=options.get('background',{})
    if bg.get('enabled',False) and getattr(record,'state',None)=='blank':
        bg={**bg,'enabled':False}
        record.history.append({'operation':'background_subtraction','status':'not_applied_to_blank_control'})
    if bg.get('enabled',False):
        if blank is None: raise ValueError('background enabled but configured blank file was not loaded')
        if channel!='dark_subtracted': raise ValueError('background subtraction requires dark_subtracted channel')
        y,details=subtract_background(record,blank,bg)
        record.stages['background_subtracted']=y.copy()
        record.history.append({'operation':'background_subtraction',**details})
    y=y[order]
    spectrum=rp.Spectrum(y,x)
    caught=[]
    with warnings.catch_warnings(record=True) as captured:
        warnings.simplefilter('always')
        despike=options.get('despike',{})
        if despike.get('enabled',True):
            spectrum=rp.preprocessing.Pipeline([rp.preprocessing.despike.WhitakerHayes(
                kernel_size=int(despike.get('kernel_size',5)),threshold=float(despike.get('threshold',10)))]).apply(spectrum)
        record.stages['despiked']=spectrum.spectral_data.copy()
        smooth=options.get('smooth',{})
        if smooth.get('enabled',True):
            window=int(smooth.get('window_length',7)); poly=int(smooth.get('polyorder',2))
            if window%2!=1 or window<=poly or window>len(x):
                raise ValueError('SavGol window must be odd, exceed polyorder and not exceed point count')
            spectrum=rp.preprocessing.Pipeline([rp.preprocessing.denoise.SavGol(window_length=window,polyorder=poly)]).apply(spectrum)
        record.stages['smoothed']=spectrum.spectral_data.copy()
        corrected=_baseline(spectrum,options.get('baseline',{}))
        record.stages['baseline']=spectrum.spectral_data-corrected.spectral_data
        norm=options.get('normalization','none')
        norms={'vector':rp.preprocessing.normalise.Vector,'minmax':rp.preprocessing.normalise.MinMax,
               'max':rp.preprocessing.normalise.MaxIntensity,'auc':rp.preprocessing.normalise.AUC,
               'area':rp.preprocessing.normalise.AUC}
        if norm not in (None,'none'):
            if norm not in norms: raise ValueError(f'unsupported normalization: {norm}')
            if np.ptp(corrected.spectral_data)<=np.finfo(float).eps:
                raise ValueError('cannot normalize flat spectrum')
            corrected=rp.preprocessing.Pipeline([norms[norm]()]).apply(corrected)
        record.stages['unified']=corrected.spectral_data.copy()
        sensitivity=[]
        sens=options.get('sensitivity',{})
        if sens.get('enabled',True):
            base=dict(options.get('baseline',{})); base.setdefault('lam',1e6)
            for method in sens.get('methods',['asls','arpls','airpls']):
                for factor in sens.get('lam_factors',[.3,1,3]):
                    settings={**base,'method':method,'lam':base['lam']*factor}
                    variant=_baseline(spectrum,settings).spectral_data.copy()
                    key=f'{method}_lam_{float(factor):g}'
                    record.stages[key]=variant
                    if float(factor)==1: record.stages[method]=variant
                    sensitivity.append({'variant':key,'method':method,'lam':settings['lam'],
                        'window_length':smooth.get('window_length',7) if smooth.get('enabled',True) else None,
                        'maximum_position_cm1':float(x[np.argmax(variant)]),
                        'integrated_intensity':float(np.trapezoid(variant,x)),
                        'negative_fraction':float(np.mean(variant<0))})
            for window in sens.get('smooth_windows',[5,7,9]):
                if not smooth.get('enabled',True) or window==smooth.get('window_length',7): continue
                if window%2!=1 or window<=smooth.get('polyorder',2) or window>len(x):
                    raise ValueError(f'invalid sensitivity SavGol window: {window}')
                varied=rp.preprocessing.Pipeline([rp.preprocessing.denoise.SavGol(
                    window_length=window,polyorder=smooth.get('polyorder',2))]).apply(rp.Spectrum(record.stages['despiked'],x))
                variant=_baseline(varied,options.get('baseline',{})).spectral_data.copy()
                key=f'smooth_{window}'
                record.stages[key]=variant
                sensitivity.append({'variant':key,'method':options.get('baseline',{}).get('method','asls'),
                    'lam':options.get('baseline',{}).get('lam',1e6),'window_length':window,
                    'maximum_position_cm1':float(x[np.argmax(variant)]),
                    'integrated_intensity':float(np.trapezoid(variant,x)),
                    'negative_fraction':float(np.mean(variant<0))})
        caught=[f'{w.category.__name__}: {w.message}' for w in captured]
    for key in list(record.stages):
        if key not in ('dark_subtracted','baseline_subtracted','background_subtracted'):
            record.stages[key]=record.stages[key][inverse]
    for key,arr in record.stages.items():
        if not np.isfinite(arr).all(): raise ValueError(f'nonfinite preprocessing stage: {key}')
    record.warnings.extend(caught)
    record.history.append({'framework':'RamanSPy','version':version('ramanspy'),
        'operation':'preprocessing','parameters':options,'channel':channel,
        'processed_at':datetime.now(timezone.utc).isoformat(),'warnings':caught,
        'axis_note':'SavGol windows are in sample points; irregular shift spacing changes physical smoothing width',
        'normalization_applies_to':'unified only; sensitivity stages preserve unnormalized intensity'})
    return sensitivity
