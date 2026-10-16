"""Model-compared mathematical band decomposition with explicit failure states.

Gaussian, Lorentzian and pseudo-Voigt peaks have height/FWHM parameters. The
pseudo-Voigt mixing parameter is the Lorentzian fraction. The fitted baseline
is linear. Analytic areas integrate the peak over the whole real axis. Neither
components nor their areas represent chemical phases or phase fractions.
"""

from __future__ import annotations

from typing import Any

import numpy as np
from scipy.optimize import least_squares
from scipy.signal import find_peaks, peak_widths

from .peaks import noise_sigma, validate_xy


DEFAULTS: dict[str, Any] = {
    "range": None, "shapes": ["gaussian", "lorentzian", "pseudo_voigt"],
    "max_components": 3, "multistart": 6, "seed": 0, "criterion": "bic",
    "min_fwhm_cm1": None, "max_fwhm_cm1": None, "min_component_snr": 3.0,
    "max_parameter_correlation": 0.98, "boundary_tolerance": 0.01,
    "stability_rtol": 0.1, "max_nfev": 3000, "initial_components": None,
}


def _profile(t: np.ndarray, center: float, width: float, height: float,
             shape: str, eta: float = 0.5) -> np.ndarray:
    q = ((t - center) / width) ** 2
    if shape == "gaussian":
        return height * np.exp(-4 * np.log(2) * q)
    lorentzian = 1 / (1 + 4 * q)
    if shape == "lorentzian":
        return height * lorentzian
    return height * (eta * lorentzian + (1 - eta) * np.exp(-4 * np.log(2) * q))


def _curve(t: np.ndarray, p: np.ndarray, shape: str, count: int) -> np.ndarray:
    value = p[0] + p[1] * (t - 0.5)
    stride = 4 if shape == "pseudo_voigt" else 3
    for j in range(count):
        i = 2 + j * stride
        value = value + _profile(t, p[i], p[i + 1], p[i + 2], shape,
                                 p[i + 3] if stride == 4 else 0.5)
    return value


def _criteria(sse: float, n: int, parameters: int, floor: float) -> dict:
    # Include one fitted common residual-variance parameter in k.
    k = parameters + 1
    base = n * np.log(max(sse / n, floor))
    return {"sse": float(sse), "rmse": float(np.sqrt(sse / n)),
            "bic": float(base + k * np.log(n)),
            "aicc": float(base + 2 * k + 2 * k * (k + 1) / (n - k - 1)) if n > k + 1 else float("inf"),
            "n_parameters": parameters, "criterion_parameters_including_variance": k}


def _parameters(p: np.ndarray, shape: str, count: int, low: float,
                span: float, scale: float) -> dict:
    stride = 4 if shape == "pseudo_voigt" else 3
    components = []
    for j in range(count):
        i = 2 + stride * j
        components.append({"center_cm1": float(low + span * p[i]),
                           "fwhm_cm1": float(span * p[i + 1]),
                           "height": float(scale * p[i + 2]),
                           "eta": float(p[i + 3]) if stride == 4 else None})
    return {"baseline_offset_at_window_midpoint": float(p[0] * scale),
            "baseline_slope_per_cm1": float(p[1] * scale / span), "components": components}


def _area(height: float, width: float, shape: str, eta: float) -> float:
    g = np.sqrt(np.pi / (4 * np.log(2)))
    coefficient = g if shape == "gaussian" else np.pi / 2 if shape == "lorentzian" else eta * np.pi / 2 + (1 - eta) * g
    return float(height * width * coefficient)


def _failure(x: np.ndarray, y: np.ndarray, cfg: dict, reason: str) -> dict:
    missing = np.full(y.shape, np.nan)
    return {"status": "failed", "components": [], "models": [], "x": x, "y": y,
            "fitted": missing, "baseline": missing.copy(), "residual": missing.copy(),
            "reasons": [reason], "diagnostics": {}, "config": cfg}


def fit_window(x: Any, y: Any, config: dict | None = None) -> dict:
    """Fit zero through max_components and classify identifiability conservatively.

    Conditional covariance assumes independent homoscedastic residuals and a
    locally linear least-squares approximation. It is withheld whenever selected
    parameters are unreliable. It never estimates sample-to-sample uncertainty.
    """
    cfg = {**DEFAULTS, **(config or {})}
    try:
        xx, yy = validate_xy(x, y)
    except (ValueError, TypeError) as exc:
        return _failure(np.asarray(x), np.asarray(y), cfg, str(exc))
    window = cfg["range"]
    if window is not None:
        if len(window) != 2 or not np.isfinite(window).all() or window[0] >= window[1]:
            return _failure(xx, yy, cfg, "range must contain finite increasing lower/upper bounds")
        inside = (xx >= window[0]) & (xx <= window[1])
        xx, yy = xx[inside], yy[inside]
    if len(xx) < 12:
        return _failure(xx, yy, cfg, "at least 12 points are required in the selected window")
    shapes = cfg["shapes"]
    if isinstance(shapes, str):
        shapes = [shapes]
    if not shapes or any(s not in ("gaussian", "lorentzian", "pseudo_voigt") for s in shapes):
        return _failure(xx, yy, cfg, "shapes must be gaussian, lorentzian or pseudo_voigt")
    criterion = str(cfg["criterion"]).lower()
    if criterion not in ("aicc", "bic"):
        return _failure(xx, yy, cfg, "criterion must be aicc or bic")
    if int(cfg["max_components"]) < 0 or int(cfg["multistart"]) < 1:
        return _failure(xx, yy, cfg, "max_components must be nonnegative and multistart positive")
    low, span = float(xx[0]), float(xx[-1] - xx[0])
    t = (xx - low) / span
    min_width = float(cfg["min_fwhm_cm1"] or 2 * np.median(np.diff(xx))) / span
    max_width = float(cfg["max_fwhm_cm1"] or span * 0.8) / span
    if not 0 < min_width < max_width:
        return _failure(xx, yy, cfg, "FWHM bounds must be positive and increasing")
    scale = max(float(np.ptp(yy)), float(np.std(yy)), np.finfo(float).eps * max(1.0, float(np.max(np.abs(yy)))))
    z = yy / scale
    baseline_p = np.linalg.lstsq(np.column_stack([np.ones_like(t), t - 0.5]), z, rcond=None)[0]
    base_curve = _curve(t, baseline_p, "none", 0)
    numerical_floor = (np.finfo(float).eps * max(1, float(np.max(np.abs(yy)))) * 64)**2
    null = {"shape": "none", "n_components": 0, "status": "converged",
            **_criteria(float(np.sum((yy - base_curve * scale)**2)), len(yy), 2, numerical_floor),
            "parameters": _parameters(baseline_p, "none", 0, low, span, scale),
            "initial_parameters": [], "constraints": {"baseline": "linear unrestricted"},
            "starts": [], "optimizer": "numpy.linalg.lstsq"}
    models = [null]
    internals: dict[tuple[str, int], dict] = {( "none", 0): {"p": baseline_p}}
    # Exact constant/linear signal needs no nonlinear optimization. A numerical
    # floor common to every model prevents exact-fit roundoff selecting peaks.
    flat = null["rmse"] <= np.sqrt(numerical_floor)
    rng = np.random.default_rng(int(cfg["seed"]))
    initial_noise = noise_sigma(xx, yy) / scale
    residual_signal = z - base_curve
    maxima, props = find_peaks(residual_signal, prominence=max(3 * initial_noise, 0.015))
    order = np.argsort(props["prominences"])[::-1] if len(maxima) else np.array([], dtype=int)
    maxima = maxima[order]
    maxima_widths = peak_widths(residual_signal, maxima, rel_height=0.5)[0] * np.median(np.diff(t)) if len(maxima) else np.array([])
    for shape in dict.fromkeys(shapes):
        for count in range(1, int(cfg["max_components"]) + 1):
            stride = 4 if shape == "pseudo_voigt" else 3
            parameter_count = 2 + stride * count
            constraints = {"baseline": "linear unrestricted", "center_cm1": [low, low + span],
                           "fwhm_cm1": [min_width * span, max_width * span],
                           "height": [0.0, 8 * scale], "eta": [0.0, 1.0] if stride == 4 else None}
            model: dict = {"shape": shape, "n_components": count, "n_parameters": parameter_count,
                           "status": "not_run", "constraints": constraints, "initial_parameters": [],
                           "starts": [], "aicc": float("inf"), "bic": float("inf"),
                           "optimizer": "scipy.optimize.least_squares"}
            models.append(model)
            if flat or len(xx) <= parameter_count + 3:
                model["status"] = "skipped_exact_linear_signal" if flat else "insufficient_degrees_of_freedom"
                continue
            lower = np.r_[[-np.inf, -np.inf], np.tile([0.0, min_width, 0.0] + ([0.0] if stride == 4 else []), count)]
            upper = np.r_[[np.inf, np.inf], np.tile([1.0, max_width, 8.0] + ([1.0] if stride == 4 else []), count)]
            solutions = []
            for start in range(int(cfg["multistart"])):
                parts = []
                supplied = cfg["initial_components"] or []
                for j in range(count):
                    if start == 0 and j < len(supplied):
                        c = supplied[j]
                        center = (float(c["center_cm1"]) - low) / span
                        width = float(c["fwhm_cm1"]) / span
                        height = float(c["height"]) / scale
                        eta = float(c.get("eta", 0.5))
                    elif j < len(maxima) and start < 2:
                        center = float(t[maxima[j]])
                        width = float(maxima_widths[j])
                        height = max(0.01, float(residual_signal[maxima[j]]))
                        eta = 0.5
                    else:
                        center = float(t[np.argmax(residual_signal)]) if start == 0 and j == 0 else float(rng.uniform(0.05, 0.95))
                        width = float(rng.uniform(max(min_width, 0.015), max(min_width * 1.01, min(max_width, 0.25))))
                        height = max(0.01, float(np.max(residual_signal)) / count)
                        eta = float(rng.uniform(0.15, 0.85))
                    if start > 0:
                        center += float(rng.normal(0, 0.02))
                        width *= float(rng.uniform(0.65, 1.5))
                    parts.extend([center, width, height] + ([eta] if stride == 4 else []))
                p0 = np.r_[baseline_p, parts]
                p0 = np.maximum(lower + 1e-9, np.minimum(upper - 1e-9, p0))
                model["initial_parameters"].append(_parameters(p0, shape, count, low, span, scale))
                try:
                    result = least_squares(lambda p: _curve(t, p, shape, count) - z, p0,
                                           bounds=(lower, upper), max_nfev=int(cfg["max_nfev"]),
                                           xtol=1e-10, ftol=1e-10, gtol=1e-10)
                    sse = float(np.sum(result.fun**2) * scale**2)
                    model["starts"].append({"success": bool(result.success), "sse": sse,
                                             "nfev": int(result.nfev), "message": str(result.message)})
                    if result.success and np.isfinite(result.x).all():
                        solutions.append((sse, result))
                except (ValueError, FloatingPointError, np.linalg.LinAlgError) as exc:
                    model["starts"].append({"success": False, "message": str(exc)})
            if not solutions:
                model["status"] = "optimization_failed"
                continue
            solutions.sort(key=lambda item: item[0])
            sse, best = solutions[0]
            model.update(_criteria(sse, len(yy), parameter_count, numerical_floor))
            model["status"] = "converged"
            model["parameters"] = _parameters(best.x, shape, count, low, span, scale)
            internals[(shape, count)] = {"p": best.x, "result": best, "solutions": solutions,
                                         "lower": lower, "upper": upper}
    available = [m for m in models if m["status"] == "converged"]
    selected = min(available, key=lambda m: (m[criterion], m["n_components"]))
    shape, count = selected["shape"], selected["n_components"]
    internal = internals[(shape, count)]
    p = internal["p"]
    fitted = _curve(t, p, shape, count) * scale
    baseline = (p[0] + p[1] * (t - 0.5)) * scale
    residual = yy - fitted
    diagnostics: dict = {"selected_shape": shape, "n_components": count, "seed": int(cfg["seed"]),
                         "criterion": criterion, "bic": selected["bic"], "aicc": selected["aicc"],
                         "rmse": selected["rmse"], "baseline_offset": float(p[0] * scale),
                         "baseline_slope_per_cm1": float(p[1] * scale / span),
                         "area_definition": "analytic_component_integral_over_entire_real_axis",
                         "parameter_correlations": [], "boundary_parameters": [],
                         "criterion_assumptions": "ordinary_least_squares_independent_homoscedastic_residuals",
                         "interpretation": "mathematical_peak_decomposition_only_no_phase_or_replicate_inference"}
    reasons: list[str] = []
    components: list[dict] = []
    covariance = None
    failed_models = [m for m in models if m["status"] == "optimization_failed"]
    diagnostics["failed_models"] = [{"shape": m["shape"], "n_components": m["n_components"]} for m in failed_models]
    if failed_models:
        reasons.append("Model comparison is incomplete because one or more requested nonlinear models failed.")
    if count == 0:
        if failed_models and len(available) == 1:
            reasons.append("All requested nonlinear models failed; absence of peaks cannot be assessed.")
            status = "failed"
        elif failed_models:
            status = "unidentifiable"
        else:
            reasons.append("The zero-component linear baseline model is preferred; no reliable peak is established.")
            status = "no_peak"
    else:
        stride = 4 if shape == "pseudo_voigt" else 3
        labels = ["baseline_offset", "baseline_slope"]
        for j in range(count):
            labels.extend([f"component_{j + 1}_{k}" for k in (["center", "fwhm", "height", "eta"] if stride == 4 else ["center", "fwhm", "height"])])
        jac = internal["result"].jac
        singular_values = np.linalg.svd(jac, compute_uv=False)
        rank = int(np.count_nonzero(singular_values > singular_values[0] * 1e-10))
        diagnostics["jacobian_rank"] = rank
        diagnostics["jacobian_condition"] = float(singular_values[0] / max(singular_values[-1], np.finfo(float).tiny))
        if rank < len(p) or not np.isfinite(diagnostics["jacobian_condition"]):
            reasons.append("Singular or near-singular covariance: parameters are not locally identifiable.")
        else:
            try:
                covariance = np.linalg.inv(jac.T @ jac) * float(np.sum((residual / scale)**2)) / (len(yy) - len(p))
                std = np.sqrt(np.maximum(np.diag(covariance), 0))
                corr = covariance / np.maximum(np.outer(std, std), np.finfo(float).tiny)
                for a in range(len(p)):
                    for b in range(a + 1, len(p)):
                        if abs(corr[a, b]) >= float(cfg["max_parameter_correlation"]):
                            diagnostics["parameter_correlations"].append({"parameters": [labels[a], labels[b]], "correlation": float(corr[a, b])})
                if diagnostics["parameter_correlations"]:
                    reasons.append("Strong parameter correlation makes the selected decomposition unreliable.")
            except np.linalg.LinAlgError:
                reasons.append("Singular covariance: uncertainty cannot be estimated reliably.")
                covariance = None
        tol = float(cfg["boundary_tolerance"])
        for i in range(2, len(p)):
            lo, hi = internal["lower"][i], internal["upper"][i]
            if min(p[i] - lo, hi - p[i]) <= tol * (hi - lo):
                diagnostics["boundary_parameters"].append(labels[i])
        if diagnostics["boundary_parameters"]:
            reasons.append("One or more selected parameters are at a fitting boundary.")
        # Compare parameter vectors only among statistically near-optimal starts;
        # sort by center first so label exchange is not treated as instability.
        sorted_parts = lambda q: np.array(sorted(np.asarray(q[2:]).reshape(count, stride), key=lambda v: v[0]))
        best_parts = sorted_parts(p)
        near = [(s, r) for s, r in internal["solutions"] if s <= selected["sse"] + max(selected["sse"] * 0.01, numerical_floor * len(yy))]
        deviations = []
        for _, alternate in near:
            parts = sorted_parts(alternate.x)
            denominator = np.maximum(np.abs(best_parts), np.array([max(min_width, 0.01), min_width, 0.01] + ([0.1] if stride == 4 else [])))
            deviations.append(float(np.max(np.abs(parts - best_parts) / denominator)))
        stability = max(deviations, default=0.0)
        diagnostics["stability"] = {"near_optimal_converged_starts": len(near), "converged_starts": len(internal["solutions"]),
                                    "requested_starts": int(cfg["multistart"]), "max_relative_parameter_deviation": stability,
                                    "assessed": int(cfg["multistart"]) > 1}
        if stability > float(cfg["stability_rtol"]):
            reasons.append("Near-optimal multistart solutions have unstable component parameters.")
        if len(internal["solutions"]) < min(2, int(cfg["multistart"])):
            reasons.append("Too few converged starts to assess parameter stability.")
        alternative = [m for m in available if m["n_components"] != count and m[criterion] - selected[criterion] < 2]
        if alternative:
            reasons.append("Competing component counts have less than 2 criterion units separation.")
        residual_noise = max(noise_sigma(xx, residual), np.sqrt(numerical_floor))
        for j in range(count):
            i = 2 + stride * j
            center, width, height = low + span * p[i], span * p[i + 1], scale * p[i + 2]
            eta = float(p[i + 3]) if stride == 4 else 0.5
            local_uncertainty = None
            if covariance is not None:
                std = np.sqrt(np.maximum(np.diag(covariance), 0))
                gradient = np.zeros_like(p)
                gradient[i + 1] = _area(height, span, shape, eta)
                gradient[i + 2] = _area(scale, width, shape, eta)
                if stride == 4:
                    gradient[i + 3] = height * width * (np.pi / 2 - np.sqrt(np.pi / (4 * np.log(2))))
                local_uncertainty = {"center_cm1": float(std[i] * span), "fwhm_cm1": float(std[i + 1] * span),
                                     "height": float(std[i + 2] * scale), "area": float(np.sqrt(max(0, gradient @ covariance @ gradient))),
                                     "eta": float(std[i + 3]) if stride == 4 else None, "sigma_level": 1,
                                     "scope": "conditional_local_fit_only_not_experimental_replicates",
                                     "assumptions": "independent_homoscedastic_residuals_local_linearization"}
            component = {"shape": shape, "center_cm1": float(center), "fwhm_cm1": float(width),
                         "height": float(height), "eta": eta if stride == 4 else None,
                         "area": _area(height, width, shape, eta), "snr": float(height / residual_noise),
                         "uncertainty": local_uncertainty,
                         "fitted": _profile(t, p[i], p[i + 1], p[i + 2], shape, eta) * scale}
            if center - xx[0] < width / 2 or xx[-1] - center < width / 2:
                reasons.append("A selected peak is truncated by the window boundary; its width and area are unreliable.")
            if component["snr"] < float(cfg["min_component_snr"]):
                reasons.append("A selected component has insufficient residual-based signal-to-noise ratio.")
            local_gaps = np.diff(xx)[np.abs(xx[1:] - center) <= width]
            if not len(local_gaps) or width < 2 * float(np.max(local_gaps)):
                reasons.append("A component is undersampled on the local Raman-shift axis.")
            components.append(component)
        components.sort(key=lambda c: c["center_cm1"])
        for j, first in enumerate(components):
            for second in components[j + 1:]:
                if second["center_cm1"] - first["center_cm1"] < 0.5 * (first["fwhm_cm1"] + second["fwhm_cm1"]):
                    first["overlap"] = second["overlap"] = True
        for component in components:
            component.setdefault("overlap", False)
        status = "unidentifiable" if reasons else "ok"
        if status != "ok":
            for component in components:
                component["uncertainty"] = None
    diagnostics["models_considered"] = len(models)
    diagnostics["max_components_reached"] = bool(count > 0 and count == int(cfg["max_components"]))
    lag1 = None
    if len(residual) >= 5 and np.std(residual[:-1]) > np.sqrt(numerical_floor) and np.std(residual[1:]) > np.sqrt(numerical_floor):
        lag1 = float(np.corrcoef(residual[:-1], residual[1:])[0, 1])
    correlated_residuals = lag1 is not None and abs(lag1) > 0.5
    diagnostics["residual_lag1_correlation"] = lag1
    diagnostics["residual_correlation_threshold"] = 0.5
    adequacy_reasons = []
    if correlated_residuals:
        adequacy_reasons.append("Residuals are strongly correlated; the independent-noise covariance assumption is not met. Smoothing and model mismatch are competing causes.")
    if diagnostics["max_components_reached"]:
        adequacy_reasons.append("The selected model reaches the configured component-count cap; adequacy against higher counts has not been established.")
    diagnostics["adequacy_reasons"] = adequacy_reasons
    diagnostics["model_adequacy"] = "not_established" if adequacy_reasons else "not_falsified_by_current_diagnostics"
    if not count:
        diagnostics["uncertainty_status"] = "not_applicable_no_components"
    elif status != "ok":
        diagnostics["uncertainty_status"] = "withheld_unidentifiable"
    elif correlated_residuals:
        diagnostics["uncertainty_status"] = "withheld_correlated_residuals"
        for component in components:
            component["uncertainty"] = None
    else:
        diagnostics["uncertainty_status"] = "conditional_local_fit_only_not_experimental_replicates"
    return {"status": status, "components": components, "models": models,
            "x": xx, "y": yy, "fitted": fitted, "baseline": baseline, "residual": residual,
            "reasons": list(dict.fromkeys(reasons)), "diagnostics": diagnostics, "config": cfg}
