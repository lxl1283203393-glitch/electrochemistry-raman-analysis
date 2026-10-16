"""Strict configuration schema, recursive defaults, and scientific guards."""
from __future__ import annotations

from copy import deepcopy
import math
from pathlib import Path
from typing import Any
import yaml


DEFAULT_CONFIG = {
    'experiment': {
        'sample_name': None, 'preparation_notes': None, 'electrolyte': None,
        'ph': None, 'reference_electrode': None, 'reference_fill': None,
        'temperature_c': None, 'potential_unit': 'V', 'scan_mode': None,
        'hold_time_s': None, 'laser_wavelength_nm': None, 'laser_power_mw': None,
        'objective_magnification': None, 'measurement_position': None,
        'ir_corrected': None, 'rhe_offset_v': None, 'allow_approximate_rhe': False,
        'optical_path_id': None,
    },
    'filename_overrides': {},
    'blank': {'files': [], 'type': None},
    'preprocessing': {
        'channel': 'dark_subtracted',
        'despike': {'enabled': True, 'threshold': 10, 'kernel_size': 5},
        'smooth': {'enabled': True, 'window_length': 7, 'polyorder': 2},
        'baseline': {'method': 'asls', 'lam': 1e6, 'p': 0.01, 'max_iter': 100, 'tol': 0.001},
        'normalization': 'none',
        'background': {'enabled': False, 'file': None, 'scale': 1.0, 'conditions_verified': False},
        'sensitivity': {'enabled': True, 'methods': ['asls', 'arpls', 'airpls'],
                        'lam_factors': [0.3, 1, 3], 'smooth_windows': [5, 7, 9]},
    },
    'peaks': {'prominence': None, 'prominence_sigma': 5, 'min_height': None,
              'min_snr': 5, 'noise_sigma_floor': None, 'distance_cm1': 8, 'min_width_cm1': 1, 'max_width_cm1': None,
              'local_noise_window_cm1': 80, 'broad_band_fwhm_cm1': 100,
              'sensitivity_factors': [0.75, 1, 1.25]},
    'fitting': {'windows': [[200, 800]], 'shapes': ['gaussian', 'lorentzian', 'pseudo_voigt'],
                'max_components': 3, 'multistart': 6, 'seed': 0, 'criterion': 'bic',
                'min_fwhm_cm1': None, 'max_fwhm_cm1': None, 'min_component_snr': 3,
                'max_parameter_correlation': 0.98, 'boundary_tolerance': 0.01,
                'stability_rtol': 0.1, 'max_nfev': 3000, 'initial_components': None},
    'plotting': {'dpi': 300, 'figsize': [10, 6], 'regions': [[200, 800], [2800, 3800]],
                 'formats': ['png', 'pdf', 'svg'], 'max_annotated_records': 6,
                 'max_fit_panels': 6, 'stack_offset': None},
    'output': {'formats': ['markdown', 'html'], 'save_arrays': True},
    'qc': {'axis_min_cm1': -1000, 'axis_max_cm1': 10000, 'raw_dark_atol': 0.001,
           'raw_dark_rtol': 0.00001, 'spike_zscore': 12, 'saturation_level': None},
    'multivariate': {'enabled': False, 'method': 'pca', 'n_components': 2,
                     'normalization': 'vector', 'seed': 42},
}


class ConfigError(ValueError):
    pass


def _number(value: Any, key: str, minimum=None, maximum=None, integer=False):
    if value is None:
        return
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ConfigError(f'{key}: expected finite number')
    if integer and not isinstance(value, int):
        raise ConfigError(f'{key}: expected integer')
    if minimum is not None and value < minimum or maximum is not None and value > maximum:
        raise ConfigError(f'{key}: outside permitted range [{minimum}, {maximum}]')


def _merge(base: dict, incoming: dict, prefix='') -> dict:
    if not isinstance(incoming, dict):
        raise ConfigError(f'{prefix or "config"}: expected mapping')
    for key, value in incoming.items():
        name = f'{prefix}.{key}' if prefix else key
        if key not in base:
            raise ConfigError(f'Unknown configuration field: {name}')
        if key == 'filename_overrides' and not prefix:
            if not isinstance(value, dict):
                raise ConfigError(f'{name}: expected filename mapping')
            base[key] = deepcopy(value)
        elif isinstance(base[key], dict):
            base[key] = _merge(base[key], value, name)
        else:
            base[key] = deepcopy(value)
    return base


def _validate(config: dict) -> None:
    def walk(obj, default, prefix=''):
        for key, value in obj.items():
            name = f'{prefix}.{key}' if prefix else key
            expected = default.get(key)
            if name == 'filename_overrides':
                continue
            if isinstance(expected, dict):
                walk(value, expected, name)
            elif isinstance(expected, bool) and type(value) is not bool:
                raise ConfigError(f'{name}: expected boolean')
            elif isinstance(expected, list) and not isinstance(value, list):
                raise ConfigError(f'{name}: expected list')
            elif isinstance(expected, (int, float)) and not isinstance(expected, bool):
                if value is None:
                    raise ConfigError(f'{name}: null is not permitted')
                _number(value, name, integer=isinstance(expected, int))
            elif isinstance(expected, str) and not isinstance(value, str):
                raise ConfigError(f'{name}: expected string')
    walk(config, DEFAULT_CONFIG)
    exp = config['experiment']
    for key in ['sample_name', 'preparation_notes', 'electrolyte', 'reference_electrode', 'reference_fill',
                'scan_mode', 'measurement_position', 'optical_path_id']:
        if exp[key] is not None and not isinstance(exp[key], str):
            raise ConfigError(f'experiment.{key}: expected string or null')
    _number(exp['ph'], 'experiment.ph', -2, 16)
    _number(exp['temperature_c'], 'experiment.temperature_c', -20, 100)
    for key in ['hold_time_s', 'laser_wavelength_nm', 'laser_power_mw', 'objective_magnification']:
        _number(exp[key], f'experiment.{key}', 0)
    _number(exp['rhe_offset_v'], 'experiment.rhe_offset_v')
    if exp['ir_corrected'] is not None and type(exp['ir_corrected']) is not bool:
        raise ConfigError('experiment.ir_corrected: expected boolean or null')
    if exp['potential_unit'] not in ['V', 'mV']:
        raise ConfigError('experiment.potential_unit: expected V or mV')
    states = {'applied_potential', 'unbiased', 'open_circuit', 'blank', 'unknown'}
    for filename, values in config['filename_overrides'].items():
        if not isinstance(filename, str) or not isinstance(values, dict):
            raise ConfigError('filename_overrides: keys must be filenames and values mappings')
        unknown = set(values) - {'potential_v', 'state', 'reason', 'blank_type'}
        if unknown:
            raise ConfigError(f'filename_overrides.{filename}: unknown fields {sorted(unknown)}')
        if not isinstance(values.get('reason'), str) or not values['reason'].strip():
            raise ConfigError(f'filename_overrides.{filename}: explicit nonempty reason required')
        _number(values.get('potential_v'), f'filename_overrides.{filename}.potential_v')
        if 'state' in values and values['state'] not in states:
            raise ConfigError(f'filename_overrides.{filename}.state: invalid state')
        if values.get('state') in {'unbiased', 'open_circuit', 'blank', 'unknown'} and values.get('potential_v') is not None:
            raise ConfigError(f'filename_overrides.{filename}: non-applied state cannot contain applied potential')
    prep = config['preprocessing']
    if prep['channel'] not in ['dark_subtracted', 'baseline_subtracted', 'raw', 'dark']:
        raise ConfigError('preprocessing.channel: unknown channel')
    smooth = prep['smooth']
    _number(smooth['window_length'], 'smooth.window_length', 3, integer=True)
    _number(smooth['polyorder'], 'smooth.polyorder', 0, integer=True)
    if smooth['window_length'] % 2 != 1 or smooth['polyorder'] >= smooth['window_length']:
        raise ConfigError('smooth: odd window_length greater than polyorder required')
    for val in prep['sensitivity']['smooth_windows']:
        _number(val, 'sensitivity.smooth_windows', 3, integer=True)
        if val % 2 != 1:
            raise ConfigError('sensitivity.smooth_windows: odd windows required')
    if prep['baseline']['method'] not in {'asls', 'arpls', 'airpls', 'none'}:
        raise ConfigError('baseline.method: use asls/arpls/airpls/none')
    if not set(prep['sensitivity']['methods']) <= {'asls', 'arpls', 'airpls'}:
        raise ConfigError('sensitivity.methods: unknown baseline method')
    _number(prep['baseline']['lam'], 'baseline.lam', 1e-12)
    _number(prep['baseline']['p'], 'baseline.p', 1e-12, 1-1e-12)
    _number(prep['baseline']['max_iter'], 'baseline.max_iter', 1, integer=True)
    _number(prep['baseline']['tol'], 'baseline.tol', 0)
    _number(prep['despike']['threshold'], 'despike.threshold', 0)
    _number(prep['despike']['kernel_size'], 'despike.kernel_size', 3, integer=True)
    if prep['despike']['kernel_size'] % 2 != 1:
        raise ConfigError('despike.kernel_size: expected odd integer')
    _number(prep['background']['scale'], 'background.scale', 0)
    for name in ['lam_factors']:
        for value in prep['sensitivity'][name]:
            _number(value, f'sensitivity.{name}', 1e-12)
    if prep['normalization'] not in {'none', 'vector', 'area', 'auc', 'max', 'minmax'}:
        raise ConfigError('preprocessing.normalization: unsupported normalization')
    fit = config['fitting']
    _number(fit['max_components'], 'fitting.max_components', 0, integer=True)
    _number(fit['seed'], 'fitting.seed', 0, integer=True)
    for key in ['multistart', 'max_nfev']:
        _number(fit[key], f'fitting.{key}', 1, integer=True)
    if not set(fit['shapes']) <= {'gaussian', 'lorentzian', 'pseudo_voigt'} or not fit['shapes']:
        raise ConfigError('fitting.shapes: expected supported nonempty list')
    if fit['criterion'] not in {'bic', 'aicc'}:
        raise ConfigError('fitting.criterion: expected bic/aicc')
    for group, low, high in [('fitting', 'min_fwhm_cm1', 'max_fwhm_cm1'), ('peaks', 'min_width_cm1', 'max_width_cm1')]:
        if config[group][low] is not None and config[group][high] is not None and config[group][low] >= config[group][high]:
            raise ConfigError(f'{group}.{low}: must be below {high}')
    if fit['initial_components'] is not None:
        if not isinstance(fit['initial_components'], list):
            raise ConfigError('fitting.initial_components: expected list or null')
        for component in fit['initial_components']:
            if not isinstance(component, dict) or set(component) - {'center_cm1', 'fwhm_cm1', 'height', 'eta'} or not {'center_cm1', 'fwhm_cm1', 'height'} <= set(component):
                raise ConfigError('fitting.initial_components: expected center_cm1/fwhm_cm1/height and optional eta')
            _number(component['center_cm1'], 'initial_components.center_cm1')
            _number(component['fwhm_cm1'], 'initial_components.fwhm_cm1', 1e-12)
            _number(component['height'], 'initial_components.height', 0)
            _number(component.get('eta'), 'initial_components.eta', 0, 1)
    for key, ranges in [('fitting.windows', fit['windows']), ('plotting.regions', config['plotting']['regions'])]:
        for bounds in ranges:
            if not isinstance(bounds, list) or len(bounds) != 2:
                raise ConfigError(f'{key}: expected [low, high] pairs')
            for value in bounds:
                _number(value, key)
            if bounds[0] >= bounds[1]:
                raise ConfigError(f'{key}: low must be below high')
    for group in ['peaks', 'fitting', 'qc']:
        for key, value in config[group].items():
            if key in {'prominence', 'prominence_sigma', 'min_snr', 'noise_sigma_floor', 'distance_cm1', 'min_width_cm1', 'max_width_cm1',
                       'local_noise_window_cm1', 'broad_band_fwhm_cm1', 'min_fwhm_cm1', 'max_fwhm_cm1',
                       'min_component_snr', 'stability_rtol', 'boundary_tolerance', 'raw_dark_atol', 'raw_dark_rtol', 'spike_zscore'}:
                _number(value, f'{group}.{key}', 0)
    _number(fit['max_parameter_correlation'], 'fitting.max_parameter_correlation', 0, 1)
    for value in config['peaks']['sensitivity_factors']:
        _number(value, 'peaks.sensitivity_factors', 1e-12)
    if not config['plotting']['formats'] or not set(config['plotting']['formats']) <= {'png', 'pdf', 'svg'}:
        raise ConfigError('plotting.formats: expected png/pdf/svg')
    if not set(config['output']['formats']) <= {'markdown', 'html'}:
        raise ConfigError('output.formats: expected markdown/html')
    if len(config['plotting']['figsize']) != 2:
        raise ConfigError('plotting.figsize: expected [width,height]')
    for value in config['plotting']['figsize']:
        _number(value, 'plotting.figsize', 0.1)
    for value in config['blank']['files']:
        if not isinstance(value, str):
            raise ConfigError('blank.files: expected filenames')
    if config['qc']['axis_min_cm1'] >= config['qc']['axis_max_cm1']:
        raise ConfigError('qc.axis_min_cm1: must be below axis_max_cm1')
    if config['multivariate']['method'] not in {'pca', 'nmf'}:
        raise ConfigError('multivariate.method: expected pca/nmf')
    if config['multivariate']['normalization'] not in {'none', 'vector', 'area', 'auc', 'max', 'minmax'}:
        raise ConfigError('multivariate.normalization: unsupported method')
    _number(config['multivariate']['n_components'], 'multivariate.n_components', 1, integer=True)


def merge_config(values: dict | None = None) -> dict:
    result = _merge(deepcopy(DEFAULT_CONFIG), {} if values is None else values)
    _validate(result)
    return result


def load_config(path: Path | None = None) -> dict:
    if path is None:
        return merge_config()
    try:
        with Path(path).open('r', encoding='utf-8-sig') as handle:
            values = yaml.safe_load(handle)
    except (OSError, yaml.YAMLError) as exc:
        raise ConfigError(f'{path}: {exc}') from exc
    if values is None:
        values = {}
    return merge_config(values)
