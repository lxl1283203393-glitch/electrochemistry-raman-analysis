"""Tests use independent simulated truth and never imply experimental repeats."""

import importlib

import numpy as np
import pytest


def fit(x, y, config=None):
    return importlib.import_module("optosky_raman.fitting").fit_window(x, y, config)


def peak(x, center, fwhm, height, shape="gaussian", eta=0.4):
    q = ((x - center) / fwhm) ** 2
    g = np.exp(-4 * np.log(2) * q)
    l = 1 / (1 + 4 * q)
    return height * (g if shape == "gaussian" else l if shape == "lorentzian" else eta * l + (1 - eta) * g)


@pytest.mark.parametrize("shape", ["gaussian", "lorentzian", "pseudo_voigt"])
def test_single_known_peak_recovers_center_width_height_area(shape):
    x = np.linspace(250, 750, 501)
    y = 2 + 0.003 * (x - 500) + peak(x, 493.2, 29.5, 40, shape)
    y += np.random.default_rng(3).normal(0, 0.04, len(x))
    result = fit(x, y, {"range": [300, 700], "shapes": [shape], "max_components": 2, "multistart": 4, "seed": 8})
    assert result["status"] == "ok", result["reasons"]
    assert len(result["components"]) == 1
    component = result["components"][0]
    assert component["center_cm1"] == pytest.approx(493.2, abs=0.15)
    assert component["fwhm_cm1"] == pytest.approx(29.5, abs=0.3)
    assert component["height"] == pytest.approx(40, rel=0.01)
    area_g = 40 * 29.5 * np.sqrt(np.pi / (4 * np.log(2)))
    area_l = 40 * 29.5 * np.pi / 2
    expected = area_g if shape == "gaussian" else area_l if shape == "lorentzian" else 0.4 * area_l + 0.6 * area_g
    assert component["area"] == pytest.approx(expected, rel=0.02)
    assert component["uncertainty"]["scope"] == "conditional_local_fit_only_not_experimental_replicates"
    assert np.allclose(result["y"] - result["fitted"], result["residual"])
    assert result["diagnostics"]["n_components"] == 1
    assert any(model["n_components"] == 0 for model in result["models"])


def test_two_resolvable_peaks_are_selected_over_single_and_three_components():
    x = np.linspace(100, 800, 701)
    y = 3 + peak(x, 360, 24, 35) + peak(x, 415, 32, 25)
    y += np.random.default_rng(8).normal(0, 0.15, len(x))
    result = fit(x, y, {"range": [270, 500], "shapes": ["gaussian"], "max_components": 3, "multistart": 5, "seed": 13})
    assert result["status"] == "ok", result["reasons"]
    assert [c["center_cm1"] for c in result["components"]] == pytest.approx([360, 415], abs=0.3)
    assert result["diagnostics"]["n_components"] == 2
    assert all("aicc" in model and "bic" in model for model in result["models"])
    assert all("constraints" in model and "initial_parameters" in model for model in result["models"])


@pytest.mark.parametrize("noisy", [False, True])
def test_flat_or_pure_noise_selects_null_model(noisy):
    x = np.linspace(100, 500, 401)
    y = -3 + 0.004 * x
    if noisy:
        y = y + np.random.default_rng(35).normal(0, 1, len(x))
    result = fit(x, y, {"shapes": ["gaussian"], "max_components": 2, "multistart": 3})
    assert result["status"] == "no_peak", result["reasons"]
    assert result["components"] == []
    assert result["diagnostics"]["n_components"] == 0


def test_window_truncated_peak_cannot_be_reported_as_reliable():
    x = np.linspace(400, 600, 401)
    y = 2 + peak(x, 402, 30, 40)
    result = fit(x, y, {"shapes": ["gaussian"], "max_components": 1, "multistart": 4})
    assert result["status"] == "unidentifiable"
    assert any("boundary" in reason or "truncated" in reason for reason in result["reasons"])
    assert all(c["uncertainty"] is None for c in result["components"])


def test_seed_reproduces_selected_model_parameters_and_residual():
    x = np.linspace(200, 700, 501)
    y = peak(x, 420, 30, 20) + np.random.default_rng(18).normal(0, 0.2, len(x))
    config = {"shapes": ["gaussian"], "max_components": 2, "multistart": 4, "seed": 17}
    a, b = fit(x, y, config), fit(x, y, config)
    assert a["status"] == b["status"] == "ok"
    assert np.array_equal(a["fitted"], b["fitted"])
    assert a["components"][0]["center_cm1"] == b["components"][0]["center_cm1"]
    assert a["diagnostics"]["seed"] == 17


def test_invalid_window_returns_file_local_failure_not_fabricated_components():
    x = np.linspace(200, 700, 501)
    result = fit(x, np.ones_like(x), {"range": [800, 900]})
    assert result["status"] == "failed"
    assert result["components"] == []
    assert any("points" in reason for reason in result["reasons"])


def test_invalid_duplicate_axis_returns_failed_status():
    result = fit([1, 1, 2, 3], [1, 2, 1, 1])
    assert result["status"] == "failed"
    assert any("strictly increasing" in reason for reason in result["reasons"])


def test_highly_correlated_overlap_is_unidentifiable_even_with_excellent_residual():
    x = np.linspace(350, 450, 401)
    y = 2 + peak(x, 397, 25, 40) + peak(x, 403, 25, 35)
    y += np.random.default_rng(2).normal(0, 0.00001, len(x))
    result = fit(x, y, {"shapes": ["gaussian"], "max_components": 2, "multistart": 4, "seed": 1})
    assert result["status"] == "unidentifiable"
    assert result["diagnostics"]["n_components"] == 2
    assert result["diagnostics"]["rmse"] < 0.001
    assert result["diagnostics"]["parameter_correlations"]
    assert any("correlation" in reason for reason in result["reasons"])
    assert all(c["uncertainty"] is None for c in result["components"])


def test_width_bound_and_negative_baseline_are_not_hidden():
    x = np.linspace(350, 450, 401)
    y = -50 + peak(x, 400, 25, 40)
    result = fit(x, y, {"shapes": ["gaussian"], "max_components": 1, "multistart": 2,
                        "min_fwhm_cm1": 30, "max_fwhm_cm1": 60})
    assert result["status"] == "unidentifiable"
    assert any("fwhm" in item for item in result["diagnostics"]["boundary_parameters"])
    assert np.array_equal(result["y"], y)
    assert np.all(result["y"] < 0)
    assert result["models"][1]["constraints"]["fwhm_cm1"] == [30, 60]


def test_optimizer_failure_is_not_misreported_as_absence_of_peaks():
    x = np.linspace(200, 700, 501)
    y = peak(x, 420, 30, 20) + np.random.default_rng(18).normal(0, 0.2, len(x))
    result = fit(x, y, {"shapes": ["gaussian"], "max_components": 2,
                        "multistart": 2, "max_nfev": 1})
    assert result["status"] == "failed"
    assert any("failed" in reason.lower() for reason in result["reasons"])
    assert any(model["status"] == "optimization_failed" for model in result["models"])


def test_correlated_residuals_withhold_uncertainty_without_inventing_experimental_error():
    x = np.linspace(250, 750, 501)
    rng = np.random.default_rng(44)
    correlated = np.zeros(len(x))
    innovations = rng.normal(0, 0.2, len(x))
    for i in range(1, len(x)):
        correlated[i] = 0.9 * correlated[i - 1] + innovations[i]
    y = 2 + peak(x, 493.2, 29.5, 40) + correlated
    result = fit(x, y, {"shapes": ["gaussian"], "max_components": 1, "multistart": 3, "seed": 8})
    assert result["status"] == "ok", result["reasons"]
    assert result["diagnostics"]["residual_lag1_correlation"] > 0.5
    assert result["diagnostics"]["uncertainty_status"] == "withheld_correlated_residuals"
    assert result["diagnostics"]["model_adequacy"] == "not_established"
    assert any("independent" in reason for reason in result["diagnostics"]["adequacy_reasons"])
    assert all(c["uncertainty"] is None for c in result["components"])
