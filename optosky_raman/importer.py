"""Read-only Optosky CSV importer with line-specific failures."""
from __future__ import annotations

import csv
from datetime import datetime
import hashlib
import io
from pathlib import Path
import re
from typing import Any
import numpy as np

from .config import merge_config
from .models import SpectrumRecord
from .potential import parse_potential, convert_to_rhe


class SpectrumImportError(ValueError):
    pass


_COLUMNS = {
    'ramanshift': 'x', 'ramanshiftcm1': 'x', 'ramanshiftcm': 'x',
    'raw': 'raw', 'dark': 'dark', 'darksubtracted': 'dark_subtracted',
    'baselinesubtracted': 'baseline_subtracted',
}


def _canonical(value: str) -> str:
    return re.sub(r'[^a-z0-9]', '', value.casefold())


def _decode(data: bytes, path: Path) -> tuple[str, str]:
    encodings = ['utf-16', 'utf-8-sig', 'gb18030', 'cp1252'] if data.startswith((b'\xff\xfe', b'\xfe\xff')) else ['utf-8-sig', 'gb18030', 'utf-16', 'cp1252']
    failures = []
    for encoding in encodings:
        try:
            decoded = data.decode(encoding)
            if '\x00' in decoded:
                continue
            return decoded, encoding
        except UnicodeError as exc:
            failures.append(f'{encoding}: {exc}')
    raise SpectrumImportError(f'{path.name}: text decoding failed ({"; ".join(failures)})')


def _float_metadata(value: str | None) -> float | None:
    if value is None:
        return None
    match = re.search(r'[+-]?(?:\d+(?:\.\d*)?|\.\d+)', value)
    return float(match.group()) if match else None


def _extract_metadata(rows: list[list[str]]) -> tuple[dict[str, Any], datetime | None, list[str]]:
    raw = {}
    for row in rows:
        if len(row) >= 2 and row[0].strip():
            key = row[0].strip()
            val = ','.join(row[1:]).strip()
            if key in raw:
                raw[key] = raw[key] + [val] if isinstance(raw[key], list) else [raw[key], val]
            else:
                raw[key] = val
    meta = {'raw_metadata': raw, 'metadata_rows': rows,
            'integration_time_s': None, 'laser_power_mw': None, 'averages': None}
    warnings = []
    acquired_at = None
    for key, value in raw.items():
        value = value[-1] if isinstance(value, list) else value
        norm = _canonical(key)
        if '积分时间' in key or norm in {'integrationtime', 'integrationtimems', 'integrationtimes', 'exposuretime', 'exposuretimems'}:
            val = _float_metadata(value)
            meta['integration_time_s'] = val / 1000 if val is not None and ('ms' in key.lower() or '毫秒' in key) else val
        elif '激光功率' in key or norm in {'laserpower', 'laserpowermw'}:
            meta['laser_power_mw'] = _float_metadata(value)
        elif '平均次数' in key or norm in {'averages', 'average', 'accumulations'}:
            val = _float_metadata(value)
            meta['averages'] = int(val) if val is not None and val.is_integer() else val
        elif '扫描时间' in key or '采集时间' in key or norm in {'scantime', 'acquisitiontime', 'acquiredat', 'datetime'}:
            for fmt in [None, '%Y/%m/%d %H:%M:%S', '%Y-%m-%d %H:%M:%S.%f', '%Y年%m月%d日 %H:%M:%S']:
                try:
                    acquired_at = datetime.fromisoformat(value) if fmt is None else datetime.strptime(value, fmt)
                    break
                except ValueError:
                    pass
            if acquired_at is None:
                warnings.append(f'unparsed_acquisition_time:{value}')
    for field in ['integration_time_s', 'laser_power_mw', 'averages']:
        if meta[field] is None:
            warnings.append(f'missing_metadata:{field}')
        elif meta[field] <= 0:
            warnings.append(f'nonpositive_metadata:{field}')
    if acquired_at is None:
        warnings.append('missing_metadata:acquisition_time')
    return meta, acquired_at, warnings


def parse_file(path: Path, batch: str | None = None, config: dict | None = None) -> SpectrumRecord:
    path = Path(path)
    cfg = merge_config(config)
    try:
        data = path.read_bytes()
    except OSError as exc:
        raise SpectrumImportError(f'{path.name}: cannot read file: {exc}') from exc
    original, encoding = _decode(data, path)
    parsed = None
    header_idx = None
    columns = None
    delimiter_used = None
    for delimiter in [',', ';', '\t']:
        try:
            rows = list(csv.reader(io.StringIO(original, newline=''), delimiter=delimiter))
        except csv.Error:
            continue
        for idx, row in enumerate(rows):
            names = [_COLUMNS.get(_canonical(v)) for v in row]
            if 'x' in names:
                duplicates = sorted({n for n in names if n and names.count(n) > 1})
                if duplicates:
                    raise SpectrumImportError(f'{path.name}: duplicate data columns on row {idx+1}: {duplicates}')
                missing = set(_COLUMNS.values()) - set(names)
                if missing:
                    raise SpectrumImportError(f'{path.name}: missing required columns on row {idx+1}: {sorted(missing)}')
                parsed, header_idx, columns, delimiter_used = rows, idx, names, delimiter
                break
        if parsed is not None:
            break
    if parsed is None:
        raise SpectrumImportError(f'{path.name}: data header containing Raman Shift and required signal channels not found')
    arrays = {name: [] for name in columns if name}
    for row_idx, row in enumerate(parsed[header_idx+1:], start=header_idx+2):
        if not row or all(not cell.strip() for cell in row):
            continue
        if len(row) != len(columns):
            raise SpectrumImportError(f'{path.name}: row {row_idx}: missing/extra values ({len(row)} cells; expected {len(columns)})')
        for col_idx, name in enumerate(columns):
            if name is None:
                continue
            cell = row[col_idx].strip()
            if not cell:
                raise SpectrumImportError(f'{path.name}: row {row_idx}, column {name}: missing numeric value')
            try:
                value = float(cell)
            except ValueError as exc:
                raise SpectrumImportError(f'{path.name}: row {row_idx}, column {name}: non-numeric value {cell!r}') from exc
            if not np.isfinite(value):
                raise SpectrumImportError(f'{path.name}: row {row_idx}, column {name}: non-finite value {cell!r}')
            arrays[name].append(value)
    arrays = {name: np.asarray(value, dtype=float) for name, value in arrays.items()}
    x = arrays.pop('x')
    if len(x) < 2:
        raise SpectrumImportError(f'{path.name}: fewer than two numeric spectrum points')
    values, counts = np.unique(x, return_counts=True)
    if np.any(counts > 1):
        raise SpectrumImportError(f'{path.name}: duplicate Raman Shift values {values[counts>1][:5].tolist()}')
    dx = np.diff(x)
    if not (np.all(dx > 0) or np.all(dx < 0)):
        raise SpectrumImportError(f'{path.name}: Raman Shift axis is not strictly monotonic')
    metadata, acquired_at, warnings = _extract_metadata(parsed[:header_idx])
    metadata.update(original_text=original, original_rows=parsed, delimiter=delimiter_used,
                    data_header_row=header_idx+1, data_columns=parsed[header_idx],
                    axis_direction='increasing' if dx[0] > 0 else 'decreasing',
                    configured_experiment=dict(cfg['experiment']))
    for key, value in cfg['experiment'].items():
        if key not in metadata or metadata[key] is None:
            metadata[key] = value
    batch = batch or path.parent.name
    override = cfg['filename_overrides'].get(f'{batch}/{path.name}', cfg['filename_overrides'].get(path.name))
    potential = parse_potential(path.name, override, potential_unit=cfg['experiment']['potential_unit'])
    # Filename values conventionally denote volts; explicit mV suffix is handled by parser.
    metadata['potential'] = potential
    metadata['blank_type'] = potential['blank_type']
    if path.name in cfg['blank']['files'] or str(path) in cfg['blank']['files']:
        metadata['blank_type'] = cfg['blank']['type'] or 'other'
        metadata['configured_blank'] = True
        potential['state'], potential['potential_v'] = 'blank', None
    metadata['rhe'] = convert_to_rhe(potential['potential_v'], cfg['experiment'])
    metadata['potential_rhe_v'] = metadata['rhe']['potential_rhe_v']
    warnings.extend(potential['warnings'])
    difference = arrays['raw'] - arrays['dark'] - arrays['dark_subtracted']
    metadata['raw_dark_equal'] = bool(np.allclose(arrays['raw'] - arrays['dark'], arrays['dark_subtracted'],
                                               atol=cfg['qc']['raw_dark_atol'], rtol=cfg['qc']['raw_dark_rtol']))
    metadata['raw_dark_max_abs_difference'] = float(np.max(np.abs(difference)))
    metadata['raw_dark_rmse'] = float(np.sqrt(np.mean(difference ** 2)))
    if not metadata['raw_dark_equal']:
        warnings.append('raw_dark_mismatch: instrument Dark Subtracted is retained; extra vendor processing unknown')
    if np.min(x) < cfg['qc']['axis_min_cm1'] or np.max(x) > cfg['qc']['axis_max_cm1']:
        warnings.append('axis_outside_configured_range')
    if dx[0] < 0:
        warnings.append('descending_axis_preserved')
    return SpectrumRecord(path=path.resolve(), batch=batch, x=x, channels=arrays, metadata=metadata,
                          state=potential['state'], potential_v=potential['potential_v'],
                          acquired_at=acquired_at, warnings=warnings,
                          sha256=hashlib.sha256(data).hexdigest(), encoding=encoding,
                          history=[{'operation': 'import', 'source_sha256': hashlib.sha256(data).hexdigest(),
                                    'encoding': encoding, 'header_row': header_idx+1, 'raw_data_modified': False}])


def load_directory(path: Path, batch: str | None = None, config: dict | None = None) -> tuple[list[SpectrumRecord], list[dict]]:
    path = Path(path)
    cfg = merge_config(config)
    if not path.is_dir():
        return [], [{'file': str(path), 'error': 'Input directory does not exist or is not a directory', 'error_type': 'DirectoryError'}]
    records, errors = [], []
    files = sorted((p for p in path.rglob('*') if p.is_file() and p.suffix.casefold() == '.csv'), key=lambda p: str(p).casefold())
    for file in files:
        try:
            record = parse_file(file, batch or path.name, cfg)
            record.metadata['source_relative_path'] = file.relative_to(path).as_posix()
            records.append(record)
        except (SpectrumImportError, OSError) as exc:
            errors.append({'file': str(file), 'error': str(exc), 'error_type': type(exc).__name__})
    if not files:
        errors.append({'file': str(path), 'error': 'No CSV files found recursively', 'error_type': 'NoDataError'})
    return records, errors
