import pytest
import yaml

from optosky_raman.config import load_config, ConfigError


def test_deep_merge_preserves_defaults_and_null_experiment_metadata(tmp_path):
    p = tmp_path / 'experiment.yaml'
    p.write_text('preprocessing:\n  smooth:\n    enabled: false\nexperiment:\n  ph: 14\n', encoding='utf-8')
    config = load_config(p)
    assert config['preprocessing']['smooth']['window_length'] == 7
    assert config['preprocessing']['smooth']['enabled'] is False
    assert config['experiment']['ph'] == 14
    assert config['experiment']['sample_name'] is None
    assert load_config()['experiment']['ph'] is None


@pytest.mark.parametrize('content', [
    {'preprocessing': {'smoth': {}}},
    {'experiment': {'ph': 20}},
    {'preprocessing': {'smooth': {'window_length': 4}}},
    {'preprocessing': {'smooth': {'enabled': 'false'}}},
    {'fitting': {'windows': [[800, 200]]}},
    {'filename_overrides': {'x.csv': {'potential_v': 1, 'invented': True}}},
    {'filename_overrides': {'x.csv': {'potential_v': -1}}},
    [],
    {'preprocessing': {'smooth': {'window_length': None}}},
    {'preprocessing': {'despike': {'kernel_size': 4}}},
    {'experiment': {'sample_name': 123}},
    {'multivariate': {'method': 'random'}},
    {'qc': {'axis_min_cm1': 100, 'axis_max_cm1': 10}},
    {'plotting': {'formats': []}},
    {'fitting': {'criterion': 'aic'}},
    {'fitting': {'seed': -1}},
    {'fitting': {'min_fwhm_cm1': 10, 'max_fwhm_cm1': 2}},
    {'fitting': {'initial_components': {'center_cm1': 2}}},
    {'fitting': {'initial_components': [{'center_cm1': 2, 'fwhm_cm1': -1, 'height': 3}]}},
    {'fitting': {'initial_components': [{'center_cm1': 2, 'fwhm_cm1': 1, 'height': 3, 'unknown': 1}]}},
])
def test_invalid_or_unknown_configuration_fails_loudly(tmp_path, content):
    p = tmp_path / 'invalid.yaml'
    p.write_text(yaml.safe_dump(content), encoding='utf-8')
    with pytest.raises(ConfigError):
        load_config(p)
