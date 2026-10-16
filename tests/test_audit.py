import numpy as np

from optosky_raman.audit import audit_records
from optosky_raman.importer import parse_file
from test_import import write_csv


def test_audit_rebuilds_internal_time_order_and_does_not_remove_suspects(tmp_path):
    p1 = write_csv(tmp_path, name='-1.2cv.csv')
    p2 = write_csv(tmp_path, name='-1.1cv.csv')
    p2.write_text(p2.read_text(encoding='utf-8').replace('13:51:42', '13:50:42'), encoding='utf-8')
    records = [parse_file(p1, 'batch'), parse_file(p2, 'batch')]
    df = audit_records(records, None)
    assert len(df) == 2
    ordered = df.sort_values('acquisition_order')
    assert list(ordered['filename']) == ['-1.1cv.csv', '-1.2cv.csv']
    assert list(ordered['interval_s'])[1] == 60
    assert list(ordered['negative_count'])[0] == 0
    assert all(df['axis_consistent'])
    assert all(df['raw_dark_equal'])


def test_mismatched_axes_and_multiplicative_drop_are_reported_without_deletion(tmp_path):
    a = write_csv(tmp_path, name='-1.1cv.csv', rows=['1,1,200,11,10', '2,1,201,21,20', '3,1,202,31,30'])
    b = write_csv(tmp_path, name='-1.2cv.csv', rows=['1,1,200,2,1', '2,1,201,3,2', '3,1,203,4,3'])
    records = [parse_file(a, 'same'), parse_file(b, 'same')]
    df = audit_records(records, None)
    assert len(df) == 2
    assert not all(df['axis_consistent'])
    assert 'signal_median' in df
    assert 'suspect_flags' in df
    assert df['similarity_max'].isna().all()


def test_audit_marks_known_setting_change_before_unknown_optical_scaling(tmp_path):
    a = write_csv(tmp_path, name='-1.1cv.csv', rows=['1,1,200,11,10', '2,1,201,21,20', '3,1,202,31,30'])
    b = write_csv(tmp_path, name='-1.2cv.csv', rows=['1,1,200,101,100', '2,1,201,201,200', '3,1,202,301,300'])
    b.write_text(b.read_text(encoding='utf-8').replace('10000','20000').replace('13:51:42','13:52:42'), encoding='utf-8')
    df = audit_records([parse_file(a, 'same'), parse_file(b, 'same')])
    flags = df.iloc[1]['suspect_flags']
    assert 'acquisition_settings_changed' in flags
    assert 'check_optics_bubbles_focus' not in flags
