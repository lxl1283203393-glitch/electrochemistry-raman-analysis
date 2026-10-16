import pytest

from optosky_raman.potential import parse_potential, convert_to_rhe


@pytest.mark.parametrize('name,state,value', [
    ('-1.353ing cv _2026_10_12_15_48_55.csv', 'applied_potential', -1.353),
    ('-1.023cv.csv', 'applied_potential', -1.023),
    ('0cv.csv', 'applied_potential', 0.0),
    ('no sample.csv', 'blank', None),
    ('OPEN CV _2026_09_06.csv', 'open_circuit', None),
    ('open circuit potential.csv', 'open_circuit', None),
    ('unbiased.csv', 'unbiased', None),
    ('sample_2026_10_12.csv', 'unknown', None),
    ('sample1.2.csv', 'unknown', None),
    ('20s 1.273_2026_09_06.csv', 'applied_potential', 1.273),
])
def test_signed_potential_and_state_are_separate(name, state, value):
    result = parse_potential(name)
    assert result['state'] == state
    assert result['potential_v'] == value


def test_unknown_suffix_has_no_invented_experiment_meaning():
    result = parse_potential('-1.353ing cv.csv')
    assert result['suffix'] == 'ing'
    assert 'unknown_suffix_meaning' in result['warnings']


def test_documented_override_preserves_positive_filename_label():
    result = parse_potential('1.023CV.csv', {'potential_v': -1.023, 'reason': 'operator confirmed missing sign'})
    assert result['potential_v'] == -1.023
    assert result['filename_potential_v'] == 1.023
    assert result['provenance']['override']['reason'] == 'operator confirmed missing sign'


def test_configured_mv_is_applied_and_explicit_unit_wins():
    assert parse_potential('-1023cv.csv', potential_unit='mV')['potential_v'] == -1.023
    assert parse_potential('-1.023V.csv', potential_unit='mV')['potential_v'] == -1.023


def test_explicit_potential_override_can_resolve_unknown_filename():
    result = parse_potential('mystery.csv', {'potential_v': -1.2, 'reason': 'operator log'})
    assert result['state'] == 'applied_potential'
    assert result['potential_v'] == -1.2
    assert result['filename_state'] == 'unknown'


def test_rhe_missing_conditions_is_null_and_explicit_offset_is_usable():
    missing = convert_to_rhe(-1.023, {})
    assert missing['potential_rhe_v'] is None
    assert missing['missing']
    offset = convert_to_rhe(-1.023, {'rhe_offset_v': 1.023})
    assert offset['potential_rhe_v'] == pytest.approx(0)


def test_known_reference_at_25c_requires_ph_and_approximate_conversion_optin():
    valid = convert_to_rhe(-1, {'reference_electrode': 'Ag/AgCl', 'reference_fill': 'saturated KCl', 'ph': 14, 'temperature_c': 25})
    assert valid['potential_rhe_v'] == pytest.approx(0.02524, abs=0.0001)
    unknown = convert_to_rhe(-1, {'reference_electrode': 'Ag/AgCl', 'ph': 14, 'temperature_c': 25})
    assert unknown['potential_rhe_v'] is None
    approx = convert_to_rhe(-1, {'reference_electrode': 'Ag/AgCl', 'ph': 14, 'temperature_c': 25, 'allow_approximate_rhe': True})
    assert approx['potential_rhe_v'] is not None
    assert approx['approximate'] is True and approx['assumptions']
