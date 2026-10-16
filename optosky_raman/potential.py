"""Filename state parsing and explicitly qualified reference conversion."""
from __future__ import annotations

from pathlib import Path
import re
import math


def parse_potential(filename: str, override: dict | None = None, potential_unit: str = 'V') -> dict:
    label = Path(filename).stem.replace('−', '-').strip()
    label = re.split(r'_20\d{2}[_-]\d{1,2}[_-]\d{1,2}', label, maxsplit=1)[0].strip()
    lower = label.casefold()
    state, value, suffix, blank_type = 'unknown', None, None, None
    warnings = []
    if re.search(r'\bno[ _-]*sample\b|无样品', lower):
        state, blank_type = 'blank', 'no_sample'
    elif re.search(r'\bbare[ _-]*ni|裸.*镍|裸.*ni', lower):
        state, blank_type = 'blank', 'bare_ni_foam'
    elif lower in {'electrolyte', '电解液', 'optical window', 'window', '光学窗口'}:
        state = 'blank'
        blank_type = 'electrolyte' if lower in {'electrolyte', '电解液'} else 'optical_window'
    elif re.search(r'open\s*(?:circuit|cv)|\bocp\b|开路', lower):
        state = 'open_circuit'
    elif re.search(r'\bunbiased\b|\bno\s+(?:bias|applied\s+potential)\b|未施加|无外加', lower):
        state = 'unbiased'
    else:
        match = re.search(r'(?<![\d.])([+-]?(?:\d+(?:\.\d+)?|\.\d+))\s*([a-z]*?)\s*cv\b', lower)
        if match is None:
            match = re.fullmatch(r'(?:\d+(?:\.\d+)?\s*(?:ms|s)\s+)?([+-]?(?:\d+\.\d+|\.\d+))(?:\s*(v|mv))?\s*', lower)
        if match:
            value = float(match.group(1))
            suffix = match.group(2) or None
            if suffix == 'mv' or suffix != 'v' and potential_unit == 'mV':
                value /= 1000
            if suffix in {'v', 'mv'}:
                suffix = None
            state = 'applied_potential'
            if suffix:
                warnings.append('unknown_suffix_meaning')
    result = {'state': state, 'potential_v': value, 'filename_potential_v': value,
              'filename_state': state, 'suffix': suffix, 'blank_type': blank_type,
              'warnings': warnings, 'provenance': {'source': 'filename', 'label': label}}
    if override:
        result['provenance']['override'] = dict(override)
        if 'potential_v' in override and override['potential_v'] is not None and 'state' not in override:
            result['state'] = 'applied_potential'
        for key in ['state', 'potential_v', 'blank_type']:
            if key in override:
                result[key] = override[key]
        if result['state'] != 'applied_potential':
            result['potential_v'] = None
        result['warnings'].append('explicit_filename_override')
    return result


def convert_to_rhe(potential_v: float | None, experiment: dict) -> dict:
    result = {'potential_rhe_v': None, 'approximate': False, 'assumptions': [],
              'missing': [], 'reference_offset_v': None, 'ir_corrected': experiment.get('ir_corrected')}
    if potential_v is None:
        result['missing'] = ['applied_potential']
        return result
    if experiment.get('rhe_offset_v') is not None:
        result['potential_rhe_v'] = potential_v + float(experiment['rhe_offset_v'])
        result['reference_offset_v'] = float(experiment['rhe_offset_v'])
        result['assumptions'] = ['RHE offset explicitly supplied by configuration; calibration validity not inferred']
        return result
    values = dict(experiment)
    approximate = bool(values.get('allow_approximate_rhe', False))
    required = ['reference_electrode', 'reference_fill', 'ph', 'temperature_c']
    result['missing'] = [key for key in required if values.get(key) is None]
    if result['missing'] and not approximate:
        result['assumptions'] = ['RHE conversion withheld: missing reference filling solution, pH or temperature']
        return result
    defaults = {'reference_electrode': 'Ag/AgCl', 'reference_fill': 'saturated KCl', 'ph': 14, 'temperature_c': 25}
    if approximate:
        for key in result['missing']:
            values[key] = defaults[key]
            result['assumptions'].append(f'Approximate conversion assumes {key}={values[key]}')
        result['approximate'] = bool(result['missing'])
    reference = re.sub(r'\s', '', str(values.get('reference_electrode', '')).lower())
    fill = re.sub(r'\s', '', str(values.get('reference_fill', '')).lower())
    if reference not in {'ag/agcl', 'agcl', '银/氯化银'}:
        result['assumptions'].append('Unsupported reference electrode: supply a calibrated rhe_offset_v')
        return result
    offsets = {'saturatedkcl': 0.197, 'satkcl': 0.197, 'sat.kcl': 0.197, '饱和kcl': 0.197,
               '3mkcl': 0.210, '3.0mkcl': 0.210, '1mkcl': 0.235, '1.0mkcl': 0.235}
    if fill not in offsets:
        result['assumptions'].append('Unknown reference filling solution: supply a calibrated rhe_offset_v')
        return result
    if float(values['temperature_c']) != 25:
        if not approximate:
            result['assumptions'].append('Reference potential table valid at 25 C only; provide calibrated rhe_offset_v at actual temperature')
            return result
        result['approximate'] = True
        result['assumptions'].append('Approximation uses 25 C reference potential at configured temperature')
    nernst = math.log(10) * 8.314462618 * (273.15 + float(values['temperature_c'])) / 96485.33212
    offset = offsets[fill] + nernst * float(values['ph'])
    result['reference_offset_v'] = offset
    result['potential_rhe_v'] = potential_v + offset
    result['assumptions'].append('Tabulated Ag/AgCl reference potential; pH supplied by configuration; liquid junction potential not measured')
    if experiment.get('ir_corrected') is None:
        result['assumptions'].append('iR correction status unknown; conversion does not apply iR correction')
    return result
