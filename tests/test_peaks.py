"""Synthetic peak features test numerical behavior, not chemical assignments."""

import importlib

import numpy as np
import pytest
from scipy.signal import savgol_filter


def detect(x, y, config=None):
    return importlib.import_module("optosky_raman.peaks").detect_peaks(x, y, config)


def gaussian(x, center, fwhm, height):
    return height * np.exp(-4 * np.log(2) * ((x - center) / fwhm) ** 2)


def test_known_gaussian_position_and_width_on_irregular_axis():
    rng = np.random.default_rng(4)
    x = np.cumsum(rng.uniform(0.25, 0.8, 1600))
    y = gaussian(x, 390, 24, 50)
    peaks = detect(x, y, {"prominence": 5, "min_snr": 0, "distance_cm1": 5})
    assert len(peaks) == 1
    assert peaks[0]["position_cm1"] == pytest.approx(390, abs=0.7)
    assert peaks[0]["fwhm_cm1"] == pytest.approx(24, abs=0.6)
    assert peaks[0]["prominence"] == pytest.approx(50, abs=0.1)
    assert peaks[0]["width_method"] == "half_prominence"


def test_distance_is_enforced_in_physical_units_on_irregular_axis():
    x = np.r_[np.arange(0, 100, 0.1), np.arange(100, 200, 2.0)]
    y = gaussian(x, 40, 2, 10) + gaussian(x, 46, 2, 8)
    peaks = detect(x, y, {"prominence": 1, "min_snr": 0, "distance_cm1": 10})
    assert len(peaks) == 1
    assert peaks[0]["position_cm1"] == pytest.approx(40)


def test_flat_and_small_noise_spectra_do_not_create_candidates():
    x = np.linspace(100, 1000, 901)
    assert detect(x, np.full(x.shape, -3.0)) == []
    noise = np.random.default_rng(25).normal(0, 1, len(x))
    assert detect(x, noise, {"min_snr": 8, "prominence_sigma": 8}) == []


def test_wide_band_is_labeled_without_claiming_hidden_components():
    x = np.linspace(100, 1000, 1801)
    peaks = detect(x, gaussian(x, 500, 180, 100), {"prominence": 5, "min_snr": 0})
    assert len(peaks) == 1
    assert peaks[0]["broad_band"] is True
    assert peaks[0]["label"] == "candidate_broad_band"
    assert any("independent" in note for note in peaks[0]["warnings"])


def test_neighboring_resolved_maxima_can_still_be_overlapping():
    x = np.linspace(100, 800, 1401)
    y = gaussian(x, 400, 36, 100) + gaussian(x, 438, 36, 100)
    peaks = detect(x, y, {"prominence": 2, "min_snr": 0, "distance_cm1": 5})
    assert len(peaks) == 2
    assert all(peak["overlap"] for peak in peaks)
    assert all(peak["method"] == "scipy.signal.find_peaks" for peak in peaks)


@pytest.mark.parametrize("x,y", [([1, 1, 2], [1, 2, 1]), ([1, 3, 2], [1, 2, 1]), ([1, 2, 3], [1, np.nan, 1])])
def test_ambiguous_axes_and_nonfinite_signal_are_rejected(x, y):
    with pytest.raises(ValueError):
        detect(x, y)


def test_raw_noise_floor_prevents_smoothed_noise_from_becoming_high_snr_candidates():
    x = np.linspace(200, 3800, 2048)
    raw_noise = np.random.default_rng(42).normal(0, 20, len(x))
    smooth_noise = savgol_filter(raw_noise, 7, 2)
    peaks = detect(x, smooth_noise, {"noise_sigma_floor": 20,
                                    "min_snr": 5, "prominence_sigma": 5})
    assert len(peaks) < 5
    assert all(item["noise_sigma"] >= 20 for item in peaks)


def test_noise_floor_is_recorded_and_does_not_erase_strong_known_band():
    x = np.linspace(200, 800, 601)
    y = gaussian(x, 500, 30, 250)
    peaks = detect(x, y, {"noise_sigma_floor": 10, "min_snr": 5, "prominence_sigma": 5})
    assert len(peaks) == 1
    assert peaks[0]["noise_sigma"] == pytest.approx(10)
    assert peaks[0]["snr"] == pytest.approx(25, rel=0.01)
    assert peaks[0]["parameters"]["noise_sigma_floor"] == 10
