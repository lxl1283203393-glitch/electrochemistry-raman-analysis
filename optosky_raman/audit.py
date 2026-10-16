"""Descriptive numeric quality audit; no spectrum is deleted."""
from __future__ import annotations

from collections import defaultdict
from datetime import datetime
import numpy as np
import pandas as pd

from .config import merge_config
from .models import SpectrumRecord


def _noise(y: np.ndarray) -> float:
    d = np.diff(y)
    if not len(d):
        return 0.0
    return float(np.median(np.abs(d - np.median(d))) / (0.67448975 * np.sqrt(2)))


def _correlation(a: np.ndarray, b: np.ndarray) -> float:
    if np.std(a) == 0 or np.std(b) == 0:
        return float('nan')
    return float(np.corrcoef(a, b)[0, 1])


def _longest_run(values: np.ndarray) -> int:
    if len(values) == 0:
        return 0
    boundaries = np.flatnonzero(np.diff(np.r_[False, values, False].astype(int)))
    return int(np.max(boundaries[1::2] - boundaries[::2])) if len(boundaries) else 0


def audit_records(records: list[SpectrumRecord], config: dict | None = None) -> pd.DataFrame:
    cfg = merge_config(config)
    if not records:
        return pd.DataFrame(columns=['record_id', 'batch', 'filename', 'suspect_flags'])
    batches = defaultdict(list)
    for record in records:
        batches[record.batch].append(record)
    orders, intervals, previous = {}, {}, {}
    for group in batches.values():
        ordered = sorted(group, key=lambda r: (r.acquired_at is None, r.acquired_at or datetime.max, r.id))
        last = None
        for order, record in enumerate(ordered, start=1):
            orders[record.id] = order if record.acquired_at else None
            intervals[record.id] = (record.acquired_at-last.acquired_at).total_seconds() if last and record.acquired_at and last.acquired_at else None
            previous[record.id] = last
            last = record
    rows = []
    for record in records:
        x, y = record.x, record.channels['dark_subtracted']
        raw, dark = record.channels['raw'], record.channels['dark']
        noise = _noise(y)
        curvature = np.diff(y, n=2)
        curve_noise = float(np.median(np.abs(curvature-np.median(curvature))) / 0.67448975) if len(curvature) else 0
        spike_count = int(np.count_nonzero(np.abs(curvature-np.median(curvature)) > cfg['qc']['spike_zscore'] * max(curve_noise, np.finfo(float).eps))) if len(curvature) else 0
        level = cfg['qc']['saturation_level']
        max_plateau = _longest_run(np.isclose(raw, np.max(raw), rtol=0, atol=1e-8))
        saturation = bool(np.any(raw >= level)) if level is not None else bool(max_plateau >= max(5, int(np.ceil(len(raw)*0.01))))
        group = batches[record.batch]
        reference_x = group[0].x
        consistent = len(x) == len(reference_x) and np.array_equal(x, reference_x)
        correlations = []
        similar_id = None
        for other in records:
            if other is record or len(other.x) != len(x) or not np.array_equal(other.x, x):
                continue
            corr = _correlation(y, other.channels['dark_subtracted'])
            if np.isfinite(corr):
                correlations.append((corr, other.id))
        max_corr = max(correlations, default=(float('nan'), None))
        suspect = []
        if not consistent:
            suspect.append('axis_differs_within_batch')
        if not record.metadata.get('raw_dark_equal', False):
            suspect.append('raw_dark_mismatch')
        if spike_count:
            suspect.append('candidate_sharp_spikes')
        if saturation:
            suspect.append('possible_saturation_or_clipping')
        if record.acquired_at is None:
            suspect.append('missing_acquisition_time')
        if intervals[record.id] == 0:
            suspect.append('identical_acquisition_timestamp')
        if record.state == 'unknown':
            suspect.append('unknown_experiment_state')
        prev = previous.get(record.id)
        ratio, prev_corr = None, None
        if prev and np.array_equal(prev.x, x):
            py = prev.channels['dark_subtracted']
            denominator = float(np.dot(py, py))
            ratio = float(np.dot(py, y) / denominator) if denominator > 0 else None
            prev_corr = _correlation(y, py)
            settings_changed = any(prev.metadata.get(key) != record.metadata.get(key) for key in ['integration_time_s', 'laser_power_mw', 'averages'])
            if settings_changed:
                suspect.append('acquisition_settings_changed')
            if not settings_changed and ratio is not None and ratio > 0 and (ratio < 0.5 or ratio > 2) and prev_corr > 0.95:
                suspect.append('possible_common_intensity_scaling_check_optics_bubbles_focus')
        row = {
            'record_id': record.id, 'file': str(record.path), 'filename': record.path.name,
            'batch': record.batch, 'sha256': record.sha256, 'encoding': record.encoding,
            'state': record.state, 'blank_type': record.metadata.get('blank_type'),
            'potential_v': record.potential_v, 'potential_rhe_v': record.metadata.get('potential_rhe_v'),
            'filename_potential_v': record.metadata.get('potential', {}).get('filename_potential_v'),
            'potential_provenance': str(record.metadata.get('potential', {}).get('provenance', {})),
            'acquired_at': record.acquired_at, 'acquisition_order': orders[record.id],
            'interval_s': intervals[record.id], 'integration_time_s': record.metadata.get('integration_time_s'),
            'laser_power_mw': record.metadata.get('laser_power_mw'), 'averages': record.metadata.get('averages'),
            'n_points': len(x), 'x_min_cm1': float(np.min(x)), 'x_max_cm1': float(np.max(x)),
            'axis_consistent': consistent, 'duplicate_x_count': int(len(x)-len(np.unique(x))),
            'raw_min': float(np.min(raw)), 'raw_max': float(np.max(raw)),
            'dark_median': float(np.median(dark)), 'dark_min': float(np.min(dark)), 'dark_max': float(np.max(dark)),
            'signal_min': float(np.min(y)), 'signal_max': float(np.max(y)), 'signal_median': float(np.median(y)),
            'signal_integral': float(np.trapezoid(y[np.argsort(x)], np.sort(x))),
            'noise_mad': noise, 'snr_estimate': float((np.percentile(y, 95)-np.percentile(y, 5))/noise) if noise > 0 else None,
            'zero_count': int(np.count_nonzero(y == 0)), 'negative_count': int(np.count_nonzero(y < 0)),
            'zero_fraction': float(np.mean(y == 0)), 'negative_fraction': float(np.mean(y < 0)),
            'spike_count': spike_count, 'saturation_suspected': saturation, 'max_plateau_points': max_plateau,
            'raw_dark_equal': record.metadata.get('raw_dark_equal'),
            'raw_dark_max_abs_difference': record.metadata.get('raw_dark_max_abs_difference'),
            'raw_dark_rmse': record.metadata.get('raw_dark_rmse'),
            'similarity_max': max_corr[0], 'most_similar_record_id': max_corr[1],
            'previous_intensity_ratio': ratio, 'previous_similarity': prev_corr,
            'suspect_flags': ';'.join(suspect), 'warnings': ';'.join(record.warnings),
            'qc_interpretation': 'Descriptive screening only; suspect flags do not identify a physical or chemical cause',
        }
        rows.append(row)
    return pd.DataFrame(rows)
