"""Conservative mathematical candidate detection on physical Raman-shift axes.

Candidate prominence widths are *not* deconvolved FWHM. No chemical assignments
are performed here. Noise estimates are descriptive single-spectrum estimates.
"""

from __future__ import annotations

from typing import Any

import numpy as np
from scipy.integrate import trapezoid
from scipy.signal import find_peaks, peak_widths


DEFAULTS: dict[str, Any] = {
    "prominence": None,
    "prominence_sigma": 5.0,
    "min_height": None,
    "min_snr": 5.0,
    "noise_sigma_floor": None,
    "distance_cm1": 8.0,
    "min_width_cm1": 1.0,
    "max_width_cm1": None,
    "local_noise_window_cm1": 80.0,
    "broad_band_fwhm_cm1": 100.0,
}


def validate_xy(x: Any, y: Any, minimum_points: int = 3) -> tuple[np.ndarray, np.ndarray]:
    """Reject ambiguity instead of silently sorting, interpolating, or dropping rows."""
    xx, yy = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    if xx.ndim != 1 or yy.ndim != 1 or xx.shape != yy.shape:
        raise ValueError("x and y must be equal-length one-dimensional arrays")
    if len(xx) < minimum_points:
        raise ValueError(f"at least {minimum_points} points are required")
    if not np.isfinite(xx).all() or not np.isfinite(yy).all():
        raise ValueError("x and y must contain only finite numeric values")
    if not np.all(np.diff(xx) > 0):
        raise ValueError("Raman-shift axis must be strictly increasing without duplicates")
    return xx, yy


def noise_sigma(x: np.ndarray, y: np.ndarray) -> float:
    """Robust local-linear interpolation residual noise for irregular sampling.

    Interior points are compared with the line through their immediate neighbors.
    The residual is normalized for the variance of independent equal-variance
    point noise. Smooth curvature can inflate this estimate, so it is not an
    independently calibrated detector noise or experimental standard deviation.
    """
    if len(y) < 5:
        return float("nan")
    a = (x[2:] - x[1:-1]) / (x[2:] - x[:-2])
    b = 1 - a
    residual = (y[1:-1] - a * y[:-2] - b * y[2:]) / np.sqrt(1 + a**2 + b**2)
    median = np.median(residual)
    estimate = float(np.median(np.abs(residual - median)) / 0.6744897501960817)
    floor = np.finfo(float).eps * max(1.0, float(np.max(np.abs(y)))) * 32
    return max(estimate, floor)


def detect_peaks(x: Any, y: Any, config: dict | None = None) -> list[dict]:
    """Find positive local maxima without assuming a number of components.

    All distance/width filtering happens after converting interpolated index
    locations back to cm-1, so irregular sampling does not change their units.
    Area is a descriptive baseline-relative trapezoidal integral between the
    SciPy prominence bases. It can include neighboring contributions.
    """
    xx, yy = validate_xy(x, y)
    cfg = {**DEFAULTS, **(config or {})}
    if np.ptp(yy) <= np.finfo(float).eps * max(1, float(np.max(np.abs(yy)))) * 32:
        return []
    floor = cfg["noise_sigma_floor"]
    if floor is not None and (not np.isfinite(floor) or float(floor) < 0):
        raise ValueError("noise_sigma_floor must be a finite nonnegative number or None")
    noise_floor = 0.0 if floor is None else float(floor)
    # Smoothing creates correlated point noise, invalidating the independent-
    # neighbor variance assumption in noise_sigma. A pre-smoothing estimate,
    # expressed in the current intensity units, conservatively anchors SNR.
    global_noise = max(noise_sigma(xx, yy), noise_floor)
    prominence = cfg["prominence"]
    if prominence is None:
        prominence = float(cfg["prominence_sigma"]) * global_noise
    indices, properties = find_peaks(yy, prominence=float(prominence), height=cfg["min_height"])
    if len(indices) == 0:
        return []
    widths, _, left_ips, right_ips = peak_widths(yy, indices, rel_height=0.5)
    del widths
    grid_indices = np.arange(len(xx))
    left_x = np.interp(left_ips, grid_indices, xx)
    right_x = np.interp(right_ips, grid_indices, xx)
    results = []
    for j, index in enumerate(indices):
        fwhm = float(right_x[j] - left_x[j])
        if fwhm < float(cfg["min_width_cm1"]):
            continue
        if cfg["max_width_cm1"] is not None and fwhm > float(cfg["max_width_cm1"]):
            continue
        local = np.abs(xx - xx[index]) <= max(float(cfg["local_noise_window_cm1"]) / 2, fwhm)
        local_noise = max(noise_sigma(xx[local], yy[local]), noise_floor) if np.count_nonzero(local) >= 7 else global_noise
        snr = float(properties["prominences"][j] / local_noise)
        if snr < float(cfg["min_snr"]):
            continue
        left, right = int(properties["left_bases"][j]), int(properties["right_bases"][j])
        baseline = np.interp(xx[left:right + 1], [xx[left], xx[right]], [yy[left], yy[right]])
        area = float(trapezoid(yy[left:right + 1] - baseline, xx[left:right + 1]))
        broad = fwhm >= float(cfg["broad_band_fwhm_cm1"])
        warnings = ["Half-prominence width is not deconvolved FWHM; area can contain neighboring bands."]
        if floor is not None:
            warnings.append("SNR is conservatively anchored to the supplied pre-smoothing noise floor in current intensity units.")
        if broad:
            warnings.append("A broad band does not establish multiple independent peaks.")
        results.append({
            "position_cm1": float(xx[index]), "height": float(yy[index]),
            "prominence": float(properties["prominences"][j]), "fwhm_cm1": fwhm,
            "area": area, "snr": snr, "noise_sigma": float(local_noise),
            "overlap": False, "broad_band": bool(broad),
            "label": "candidate_broad_band" if broad else "candidate_peak",
            "width_method": "half_prominence", "method": "scipy.signal.find_peaks",
            "parameters": dict(cfg), "warnings": warnings,
            "left_half_prominence_cm1": float(left_x[j]),
            "right_half_prominence_cm1": float(right_x[j]),
            "left_base_cm1": float(xx[left]), "right_base_cm1": float(xx[right]),
        })
    # Find_peaks distance is in samples; physical-coordinate nonmaximum
    # suppression keeps the more prominent maximum with deterministic ties.
    selected: list[dict] = []
    for item in sorted(results, key=lambda v: (-v["prominence"], v["position_cm1"])):
        if all(abs(item["position_cm1"] - other["position_cm1"]) >= float(cfg["distance_cm1"]) for other in selected):
            selected.append(item)
    selected.sort(key=lambda v: v["position_cm1"])
    for j, first in enumerate(selected):
        for second in selected[j + 1:]:
            separation = second["position_cm1"] - first["position_cm1"]
            if separation < 0.75 * (first["fwhm_cm1"] + second["fwhm_cm1"]):
                first["overlap"] = second["overlap"] = True
    for item in selected:
        if item["overlap"]:
            item["warnings"].append("Neighboring maxima may overlap; detection is not component resolution.")
    return selected
