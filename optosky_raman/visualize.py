"""Traceable Raman figures. Rendering never changes analytical arrays."""
from __future__ import annotations

from pathlib import Path
from typing import Any, Iterable
import json
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.gridspec import GridSpec
import numpy as np
import pandas as pd


CATEGORIES = (
    ("raw_full", "Original full spectra"),
    ("processed_full", "Unified preprocessing: full spectra"),
    ("potential_stack", "Spectra ordered by original potential"),
    ("time_stack", "Spectra ordered by acquisition time"),
    ("regional_zoom", "Configured spectral regions"),
    ("blank_comparison", "Blank and sample comparison"),
    ("baseline_methods", "Baseline method sensitivity"),
    ("candidate_peaks", "Detected candidate bands"),
    ("peak_fits", "Mathematical peak fits and residuals"),
    ("peak_position_potential", "Fitted position versus original potential"),
    ("peak_area_potential", "Fitted area versus original potential"),
    ("fwhm_potential", "Fitted FWHM versus original potential"),
    ("adjacent_differences", "Adjacent acquisition difference spectra"),
    ("potential_heatmap", "Potential-ordered spectral heatmap"),
    ("anomaly_comparison", "Flagged and other spectra"),
    ("batch_comparison", "Cross-batch descriptive comparison"),
)
_STAGE_NAMES = {"dark_subtracted": "A: instrument Dark Subtracted", "baseline_subtracted": "B: instrument BaseLine Subtracted",
                "unified": "C: unified preprocessing", "asls": "AsLS", "arpls": "arPLS", "airpls": "airPLS"}
_AVAILABLE_FONTS = {entry.name for entry in font_manager.fontManager.ttflist}
_FONT = next((name for name in ("Microsoft YaHei", "SimHei", "Noto Sans CJK SC", "DejaVu Sans")
              if name in _AVAILABLE_FONTS), "DejaVu Sans")
_STYLE = {"font.family": _FONT, "font.size": 9, "axes.labelsize": 10,
          "axes.titlesize": 11, "legend.fontsize": 7, "axes.spines.top": False,
          "axes.spines.right": False, "axes.linewidth": .7, "svg.fonttype": "none", "pdf.fonttype": 42}


def _array(record: Any, stage: str) -> np.ndarray | None:
    values = record.stages.get(stage) if stage in record.stages else record.channels.get(stage)
    if values is None:
        return None
    values = np.asarray(values, dtype=float)
    return values if values.shape == np.asarray(record.x).shape and np.isfinite(values).all() else None


def _label(record: Any, reference: str) -> str:
    potential = f"{record.potential_v:+.3f} V" if record.potential_v is not None else record.state
    return f"{record.batch} | {record.path.name} | {potential}"


def _panel_title(record: Any) -> str:
    """Compact panels; exact filenames remain in CSV and manifest record IDs."""
    potential = f"{record.potential_v:+.3f} V" if record.potential_v is not None else record.state
    if re.search(r"ing\b|ing[ _]", record.path.stem, flags=re.IGNORECASE):
        potential += " [ing]"
    clock = record.acquired_at.strftime("%H:%M:%S") if record.acquired_at is not None else "time unknown"
    return f"{record.batch} | {potential} | {clock}"


def _frame(record: Any, stage: str, y: np.ndarray, **columns: Any) -> pd.DataFrame:
    frame = pd.DataFrame({"record_id": record.id, "batch": record.batch, "filename": record.path.name,
                          "source_sha256": record.sha256, "potential_v": record.potential_v,
                          "state": record.state, "acquired_at": str(record.acquired_at or ""),
                          "raman_shift_cm1": record.x, "intensity": y,
                          "processing": _STAGE_NAMES.get(stage, stage)})
    for name, values in columns.items():
        frame[name] = values
    return frame


def _decorate(ax: Any, title: str, ylabel: str = "Intensity (instrument a.u.)") -> None:
    ax.set(title=title, xlabel=r"Raman shift (cm$^{-1}$)", ylabel=ylabel)
    ax.tick_params(direction="out", length=3)
    ax.grid(alpha=.15, linewidth=.4)


def _concat(frames: list[pd.DataFrame]) -> pd.DataFrame:
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def create_figures(records: Iterable[Any], audit_df: pd.DataFrame, peaks_df: pd.DataFrame,
                   fits: list[dict], comparisons: dict, output: Path, config: dict) -> list[dict]:
    """Write up to 16 figure groups and their exact source tables.

    Every category gets a manifest record, including skipped/error categories.
    ``path`` and ``data_path`` are relative to *output*. Limited multi-panel
    figures name selected/omitted record IDs; full spectra are never downsampled.
    """
    records = list(records)
    output = Path(output)
    folder = output / "figures"
    data_folder = folder / "data"
    data_folder.mkdir(parents=True, exist_ok=True)
    plotting = config.get("plotting", {})
    formats = plotting.get("formats", ["png", "pdf", "svg"])
    dpi = int(plotting.get("dpi", 300))
    reference = config.get("experiment", {}).get("reference_electrode") or "unspecified reference"
    reference = str(reference)
    max_panels = max(1, int(plotting.get("max_annotated_records", 6)))
    fit_limit = max(1, min(int(plotting.get("max_fit_panels", 6)), int(plotting.get("max_fit_plots", 12))))
    regions = plotting.get("regions", [[200, 800], [2800, 3800]])
    preprocessing = config.get("preprocessing", {})
    c_description = ""
    if preprocessing:
        smooth = preprocessing.get("smooth", {})
        base = preprocessing.get("baseline", {})
        smoothing = f"SavGol({smooth.get('window_length')},{smooth.get('polyorder')})" if smooth.get("enabled") else "no smoothing"
        c_description = (f"C settings: {smoothing}; baseline={base.get('method', 'unknown')}, "
                         f"lam={base.get('lam', 'unknown')}; normalization={preprocessing.get('normalization', 'none')}.")
    by_id = {record.id: record for record in records}
    colors = {record.id: plt.get_cmap("tab20")(index % 20) for index, record in enumerate(records)}
    manifest: list[dict] = []

    def skipped(category: str, title: str, reason: str) -> None:
        manifest.append({"category": category, "title": title, "status": "skipped", "reason": reason,
                         "path": None, "paths": {}, "data_path": None, "reference": reference})

    def save(category: str, title: str, fig: Any, frame: pd.DataFrame, notes: str = "",
             processing: str = "", selected: list[str] | None = None) -> None:
        frame = frame.copy()
        frame["potential_reference"] = reference
        data_path = data_folder / f"{category}.csv"
        frame.to_csv(data_path, index=False, encoding="utf-8-sig")
        caption = f"Original potentials: V vs {reference}. {processing}"
        if "C:" in processing or "unified" in processing.lower():
            caption += f"\n{c_description}" if c_description else ""
        fig.text(.01, .003, caption + (f"\n{notes}" if notes else ""), fontsize=7, va="bottom")
        # Reserve physical space for the caption plus tick/axis labels. A fixed
        # fractional margin fails for both tall panels and short heatmaps.
        fig.subplots_adjust(bottom=min(.3, .98 / fig.get_figheight()), hspace=.46, wspace=.30)
        paths = {}
        try:
            for fmt in formats:
                if fmt not in {"png", "pdf", "svg"}:
                    raise ValueError(f"Unsupported figure format: {fmt}")
                filename = folder / f"{category}.{fmt}"
                fig.savefig(filename, dpi=dpi, bbox_inches="tight", facecolor="white")
                paths[fmt] = filename.relative_to(output).as_posix()
        finally:
            plt.close(fig)
        manifest.append({"category": category, "title": title, "status": "generated", "reason": "",
                         "path": paths.get("png", next(iter(paths.values()), None)), "paths": paths,
                         "data_path": data_path.relative_to(output).as_posix(), "reference": reference,
                         "processing": processing, "parameters": json.loads(json.dumps({
                             "preprocessing": preprocessing, "peaks": config.get("peaks", {}),
                             "fitting": config.get("fitting", {}), "plotting": plotting}, default=str)),
                         "notes": notes, "selected_record_ids": selected or []})

    def overlay(category: str, title: str, selected: list[Any], stage: str = "unified", stack: bool = False,
                flagged: set[str] | None = None) -> None:
        selected = [r for r in selected if _array(r, stage) is not None]
        if not selected:
            skipped(category, title, f"No finite {stage} spectra available")
            return
        fig, ax = plt.subplots(figsize=(9.5, 5.3))
        frames = []
        span = max(float(np.ptp(_array(r, stage))) for r in selected)
        offset = plotting.get("stack_offset")
        offset = float(offset) if offset is not None else (span * .65 if span else 1.)
        for index, record in enumerate(selected):
            y = _array(record, stage)
            delta = index * offset if stack else 0.
            line = "--" if record.state == "blank" or (flagged is not None and record.id not in flagged) else "-"
            label = _label(record, reference)
            if category == "time_stack":
                label = f"{record.acquired_at.isoformat(sep=' ')} | {label}"
            ax.plot(record.x, y + delta, color=colors[record.id], lw=.9, ls=line, label=label)
            frames.append(_frame(record, stage, y, displayed_intensity=y + delta, stack_offset=delta))
        _decorate(ax, title, "Intensity + display offset (a.u.)" if stack else "Intensity (instrument a.u.)")
        ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1), ncols=max(1, len(selected) // 22 + 1), frameon=False)
        notes = (f"Display offset = {offset:.5g} a.u.; no normalization." if stack else
                 "Descriptive overlay; acquisition settings/optical geometry may differ. No intensity correction.")
        save(category, title, fig, _concat(frames), notes, _STAGE_NAMES.get(stage, stage), [r.id for r in selected])

    def renderer(category: str, title: str) -> None:
        if category == "raw_full":
            overlay(category, title, records, "dark_subtracted")
        elif category == "processed_full":
            overlay(category, title, records)
        elif category == "potential_stack":
            selected = sorted([r for r in records if r.potential_v is not None and r.state != "blank"],
                              key=lambda r: (r.batch, r.potential_v, str(r.acquired_at), r.id))
            overlay(category, title, selected, stack=True)
        elif category == "time_stack":
            selected = sorted([r for r in records if r.acquired_at is not None], key=lambda r: (r.acquired_at, r.id))
            overlay(category, title, selected, stack=True)
        elif category == "regional_zoom":
            valid = [(float(lo), float(hi)) for lo, hi in regions if any(np.any((r.x >= lo) & (r.x <= hi)) for r in records)]
            selected = [r for r in records if _array(r, "unified") is not None]
            if not selected or not valid:
                skipped(category, title, "No unified spectra inside configured zoom regions")
                return
            fig, axes = plt.subplots(1, len(valid), figsize=(5.5 * len(valid), 4.5), squeeze=False)
            frames = []
            for ax, (lo, hi) in zip(axes.flat, valid):
                for record in selected:
                    mask = (record.x >= lo) & (record.x <= hi)
                    y = _array(record, "unified")
                    ax.plot(record.x[mask], y[mask], lw=.8, color=colors[record.id], label=_label(record, reference))
                    frame = _frame(record, "unified", y, region=f"{lo:g}-{hi:g}")
                    frames.append(frame.loc[mask])
                _decorate(ax, f"{lo:g}–{hi:g} cm$^{{-1}}$")
            axes[0, -1].legend(loc="upper left", bbox_to_anchor=(1.01, 1), frameon=False)
            save(category, title, fig, _concat(frames), "Regions are configured, not chemical assignments.", "C: unified preprocessing")
        elif category == "blank_comparison":
            corrected = [r for r in records if _array(r, "background_subtracted") is not None]
            if not corrected and (not any(r.state == "blank" for r in records) or not any(r.state != "blank" for r in records)):
                skipped(category, title, "Requires both an explicitly identified blank and sample")
            else:
                if not corrected:
                    overlay(category, title + " (no subtraction)", records, "dark_subtracted")
                    return
                fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
                frames = []
                for record in records:
                    y = _array(record, "dark_subtracted")
                    if y is None:
                        continue
                    axes[0].plot(record.x, y, lw=.8, ls="--" if record.state == "blank" else "-",
                                 color=colors[record.id], label=_label(record, reference))
                    frames.append(_frame(record, "dark_subtracted", y, panel="before_subtraction"))
                for record in corrected:
                    before, after = _array(record, "dark_subtracted"), _array(record, "background_subtracted")
                    if before is not None:
                        axes[1].plot(record.x, before, lw=.7, ls="--", color=colors[record.id], label=f"{record.path.name}: before")
                    details = next((h for h in reversed(record.history) if h.get("operation") == "background_subtraction"), {})
                    scale = details.get("scale", config.get("preprocessing", {}).get("background", {}).get("scale", "unrecorded"))
                    axes[1].plot(record.x, after, lw=.9, color=colors[record.id], label=f"{record.path.name}: after (scale={scale})")
                    frames.append(_frame(record, "background_subtracted", after, panel="after_subtraction",
                                         background_scale=scale, blank_file=details.get("blank_file", "unrecorded")))
                _decorate(axes[0], "A before subtraction; blank source recorded in CSV")
                _decorate(axes[1], "Verified background subtraction: before/after")
                axes[1].axhline(0, color="black", lw=.5)
                for ax in axes:
                    ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1), frameon=False)
                save(category, title, fig, _concat(frames), "Subtraction is explicit; negative dips retained and require artifact review.",
                     "A: instrument Dark Subtracted; background_subtracted uses recorded scale and verified conditions.",
                     [r.id for r in corrected])
        elif category == "baseline_methods":
            selected = [r for r in records if sum(_array(r, key) is not None for key in ("asls", "arpls", "airpls")) >= 2][:max_panels]
            if not selected:
                skipped(category, title, "Fewer than two baseline-method stages are available")
                return
            fig, axes = plt.subplots(int(np.ceil(len(selected) / 2)), min(2, len(selected)),
                                     figsize=(11, 3.4 * int(np.ceil(len(selected) / 2))), squeeze=False)
            frames = []
            for ax, record in zip(axes.flat, selected):
                for stage in ("dark_subtracted", "baseline_subtracted", "asls", "arpls", "airpls"):
                    y = _array(record, stage)
                    if y is not None:
                        ax.plot(record.x, y, lw=.8, label=_STAGE_NAMES[stage])
                        frames.append(_frame(record, stage, y))
                _decorate(ax, _panel_title(record))
                ax.legend(frameon=False)
                if len(selected) > 2 and list(axes.flat).index(ax) < len(selected) - 2:
                    ax.set_xlabel("")
            for ax in list(axes.flat)[len(selected):]:
                ax.set_visible(False)
            save(category, title, fig, _concat(frames), f"First {len(selected)} eligible records shown; selection IDs in manifest.",
                 "A/B and available AsLS/arPLS/airPLS outputs; no extra smoothing.", [r.id for r in selected])
        elif category == "candidate_peaks":
            ids = set(peaks_df.get("record_id", pd.Series(dtype=str)).dropna())
            selected = [r for r in records if r.id in ids and _array(r, "unified") is not None][:max_panels]
            if not selected or "position_cm1" not in peaks_df:
                skipped(category, title, "No candidate peaks with traceable unified spectra")
                return
            fig, axes = plt.subplots(int(np.ceil(len(selected) / 2)), min(2, len(selected)),
                                     figsize=(11, 3.4 * int(np.ceil(len(selected) / 2))), squeeze=False)
            frames = []
            for ax, record in zip(axes.flat, selected):
                y = _array(record, "unified")
                ax.plot(record.x, y, color=colors[record.id], lw=.8, label="C: unified")
                peak_frame = peaks_df.loc[peaks_df.record_id == record.id].copy()
                peak_frame = peak_frame.loc[np.isfinite(peak_frame.position_cm1.to_numpy(float))]
                positions = peak_frame.position_cm1.to_numpy(float)
                order = np.argsort(record.x)
                marker_y = np.interp(positions, record.x[order], y[order])
                ax.scatter(positions, marker_y, marker="v", s=20,
                           color="#b13749", label="Candidate bands")
                _decorate(ax, _panel_title(record))
                ax.legend(frameon=False)
                if len(selected) > 2 and list(axes.flat).index(ax) < len(selected) - 2:
                    ax.set_xlabel("")
                frame = _frame(record, "unified", y, row_type="spectrum")
                peak_frame["row_type"] = "candidate_peak"
                peak_frame["plotted_intensity"] = marker_y
                peak_frame["filename"] = record.path.name
                peak_frame["batch"] = record.batch
                peak_frame["processing"] = "C: unified preprocessing"
                peak_frame["source_sha256"] = record.sha256
                frames.extend([frame, peak_frame])
            for ax in list(axes.flat)[len(selected):]:
                ax.set_visible(False)
            save(category, title, fig, _concat(frames), f"First {len(selected)} peak-bearing records; bands have no automatic phase assignment.",
                 "C: unified preprocessing", [r.id for r in selected])
        elif category == "peak_fits":
            eligible = [fit for fit in fits if fit.get("record_id") in by_id and len(fit.get("x", []))
                        and np.isfinite(np.asarray(fit.get("fitted", []), dtype=float)).any()]
            selected = eligible[:fit_limit]
            if not selected:
                skipped(category, title, "No finite fitted spectra are available; failed fits remain in fit tables")
                return
            rows = int(np.ceil(len(selected) / 2))
            fig = plt.figure(figsize=(11, 4.7 * rows))
            grid = GridSpec(rows * 2, min(2, len(selected)), figure=fig, height_ratios=[3, 1] * rows, hspace=.70, wspace=.32)
            frames = []
            for index, fit in enumerate(selected):
                row, col = divmod(index, min(2, len(selected)))
                ax = fig.add_subplot(grid[row * 2, col])
                residual_ax = fig.add_subplot(grid[row * 2 + 1, col], sharex=ax)
                record = by_id[fit["record_id"]]
                x, y, fitted = (np.asarray(fit[key]) for key in ("x", "y", "fitted"))
                residual = np.asarray(fit["residual"])
                ax.plot(x, y, lw=.8, color="#444444", label="Input C")
                ax.plot(x, fitted, lw=1.2, color="#b13749", label=f"Fit ({fit.get('status')})")
                frame = pd.DataFrame({"record_id": record.id, "filename": record.path.name, "source_sha256": record.sha256,
                                      "raman_shift_cm1": x, "intensity": y, "fitted": fitted, "residual": residual,
                                      "processing": "C: mathematical peak fit", "fit_status": fit.get("status"),
                                      "window": str(fit.get("window", ""))})
                baseline = np.asarray(fit.get("baseline", np.zeros_like(x)))
                frame["fitted_baseline"] = baseline
                for number, component in enumerate(fit.get("components", []), 1):
                    values = component.get("fitted")
                    if values is not None and len(values) == len(x):
                        values = np.asarray(values)
                        ax.plot(x, values + baseline, ls="--", lw=.7, label=f"Component {number} + baseline")
                        frame[f"component_{number}"] = values
                residual_ax.plot(x, residual, color="#2a788e", lw=.8)
                residual_ax.axhline(0, color="black", lw=.5)
                _decorate(ax, _panel_title(record))
                ax.set_xlabel("")
                _decorate(residual_ax, "", "Residual (a.u.)")
                if row < rows - 1:
                    residual_ax.set_xlabel("")
                ax.legend(frameon=False)
                frames.append(frame)
            save(category, title, fig, _concat(frames), f"Shown {len(selected)}/{len(eligible)} finite fits; status shown, not phase amounts.",
                 "C: unified preprocessing; fitted functions are mathematical components.", [f["record_id"] for f in selected])
        elif category in {"peak_position_potential", "peak_area_potential", "fwhm_potential"}:
            field, ylabel = {"peak_position_potential": ("center_cm1", r"Fitted position (cm$^{-1}$)"),
                             "peak_area_potential": ("area", r"Fitted area (a.u. cm$^{-1}$)"),
                             "fwhm_potential": ("fwhm_cm1", r"Fitted FWHM (cm$^{-1}$)")}[category]
            rows = []
            for fit in fits:
                record = by_id.get(fit.get("record_id"))
                if fit.get("status") != "ok" or record is None or record.potential_v is None or record.state == "blank":
                    continue
                for index, component in enumerate(fit.get("components", [])):
                    if field in component and np.isfinite(component[field]):
                        rows.append({"record_id": record.id, "filename": record.path.name, "source_sha256": record.sha256,
                                     "batch": record.batch, "potential_v": record.potential_v, "component_index": index,
                                     "window": str(fit.get("window", "")), field: component[field], "fit_status": "ok",
                                     "processing": "C: reliable mathematical fit"})
            frame = pd.DataFrame(rows)
            if frame.empty:
                skipped(category, title, "No reliable fitted components at known applied potentials")
                return
            fig, ax = plt.subplots(figsize=(7, 4.7))
            for batch, group in frame.groupby("batch", sort=False):
                ax.scatter(group.potential_v, group[field], s=22, alpha=.85, label=str(batch))
            ax.set(title=title, xlabel=f"Original potential (V vs {reference})", ylabel=ylabel)
            ax.legend(title="Batch", frameon=False)
            ax.grid(alpha=.15)
            save(category, title, fig, frame, "Unconnected components; no automatic band correspondence, replicate statistics or error bars.",
                 "Only status=ok mathematical fits; descriptive observations.", list(frame.record_id.unique()))
        elif category == "adjacent_differences":
            frame = comparisons.get("difference_spectra", pd.DataFrame())
            required = {"record_id", "previous_id", "raman_shift_cm1", "difference"}
            if frame.empty or not required.issubset(frame.columns):
                skipped(category, title, "No condition-compatible adjacent difference spectra")
                return
            fig, ax = plt.subplots(figsize=(9.5, 5.3))
            frame = frame.copy()
            for (record_id, previous_id), group in frame.groupby(["record_id", "previous_id"], sort=False):
                record = by_id.get(record_id)
                label = f"{record_id} minus {previous_id}"
                ax.plot(group.raman_shift_cm1, group.difference, lw=.8, label=label,
                        color=colors.get(record_id))
                if record is not None:
                    frame.loc[group.index, "source_sha256"] = record.sha256
            frame["processing"] = "C: current minus previous acquisition; compatible-settings pairs only"
            _decorate(ax, title, "Difference intensity (a.u.)")
            ax.axhline(0, color="black", lw=.5)
            ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1), frameon=False)
            save(category, title, fig, frame, "Acquisition order; adjacent potentials may also differ in time.",
                 "C: unified preprocessing; differences supplied by comparison module.")
        elif category == "potential_heatmap":
            selected = sorted([r for r in records if r.potential_v is not None and r.state != "blank" and _array(r, "unified") is not None],
                              key=lambda r: (r.batch, r.potential_v, str(r.acquired_at), r.id))
            if len(selected) < 2:
                skipped(category, title, "Requires at least two unified spectra with known potentials")
                return
            x = selected[0].x
            if not all(len(r.x) == len(x) and np.allclose(r.x, x, rtol=0, atol=1e-7) for r in selected):
                skipped(category, title, "Wave-number axes differ; heatmap resampling is not performed implicitly")
                return
            matrix = np.vstack([_array(r, "unified") for r in selected])
            fig, ax = plt.subplots(figsize=(9, max(4.5, len(selected) * .2)))
            mesh = ax.pcolormesh(x, np.arange(len(selected)), matrix, shading="nearest", cmap="cividis", rasterized=True)
            ax.set_yticks(np.arange(len(selected)), [_panel_title(r) for r in selected], fontsize=7)
            ax.set(title=title, xlabel=r"Raman shift (cm$^{-1}$)", ylabel=f"Original potential (V vs {reference})")
            fig.colorbar(mesh, ax=ax, label="C intensity (a.u.)")
            frames = [_frame(r, "unified", _array(r, "unified"), heatmap_row=i) for i, r in enumerate(selected)]
            save(category, title, fig, _concat(frames), "Record rows, not uniform potential spacing; duplicates retained; no interpolation or normalization.",
                 "C: unified preprocessing", [r.id for r in selected])
        elif category == "anomaly_comparison":
            flagged: set[str] = set()
            if "record_id" in audit_df:
                for _, row in audit_df.iterrows():
                    flags = any(pd.notna(row.get(key)) and str(row.get(key)).strip() not in {"", "nan", "None", "False", "[]"}
                                for key in ("suspect_flags", "flags"))
                    flags = flags or any(row.get(key, False) is True or str(row.get(key, "")).lower() == "true"
                                         for key in ("suspicious", "suspect", "saturation_suspected"))
                    if flags:
                        flagged.add(row.record_id)
            if not flagged:
                skipped(category, title, "No spectra explicitly flagged by the quality audit")
            else:
                overlay(category, title + " (solid: flagged; dashed: others)", records, "dark_subtracted", flagged=flagged)
        elif category == "batch_comparison":
            if len({r.batch for r in records}) < 2:
                skipped(category, title, "Only one experimental batch is available")
            else:
                overlay(category, title + " (settings may differ)", records)

    with plt.rc_context(_STYLE):
        for category, title in CATEGORIES:
            try:
                renderer(category, title)
            except Exception as error:
                # Batch failures remain explicit and reportable, never counted as success.
                plt.close("all")
                manifest.append({"category": category, "title": title, "status": "error", "reason": f"{type(error).__name__}: {error}",
                                 "path": None, "paths": {}, "data_path": None, "reference": reference})
    return manifest
