from types import SimpleNamespace
import numpy as np
import pytest
from optosky_raman.preprocessing import process_record, subtract_background


def make_record():
    x = np.linspace(200, 900, 351)
    y = -4 + 0.015*x + 70*np.exp(-4*np.log(2)*((x-495)/45)**2)
    return SimpleNamespace(x=x,channels={'dark_subtracted':y.copy(),'baseline_subtracted':y-5},
        stages={},history=[],warnings=[],metadata={},path='synthetic.csv')


@pytest.mark.parametrize('method',['asls','arpls','airpls'])
def test_actual_ramanspy_stages_preserve_original_and_peak(method):
    rec=make_record(); before=rec.channels['dark_subtracted'].copy()
    process_record(rec, {'baseline':{'method':method,'lam':1e6},'despike':{'enabled':False},
                        'smooth':{'enabled':True,'window_length':7,'polyorder':2},
                        'sensitivity':{'enabled':False},'normalization':'none'})
    np.testing.assert_array_equal(rec.channels['dark_subtracted'],before)
    np.testing.assert_array_equal(rec.stages['dark_subtracted'],before)
    assert abs(rec.x[np.argmax(rec.stages['unified'])]-495)<5
    assert rec.history[-1]['framework']=='RamanSPy'
    assert np.isfinite(rec.stages['unified']).all()


def test_no_negative_clipping():
    rec=make_record(); rec.channels['dark_subtracted']-=100
    process_record(rec,{'baseline':{'method':'none'},'despike':{'enabled':False},
                       'smooth':{'enabled':False},'sensitivity':{'enabled':False}})
    assert np.min(rec.stages['unified'])<0


def test_blank_subtraction_requires_matching_conditions_and_axis():
    a,b=make_record(),make_record()
    with pytest.raises(ValueError,match='conditions'):
        subtract_background(a,b,{'enabled':True,'conditions_verified':False})
    b.x=b.x+0.1
    with pytest.raises(ValueError,match='axis'):
        subtract_background(a,b,{'enabled':True,'conditions_verified':True})


def test_background_subtraction_retains_negative_artifact_diagnostics():
    a,b=make_record(),make_record()
    out,details=subtract_background(a,b,{'enabled':True,'conditions_verified':True,'scale':2})
    assert (out<0).any()
    assert details['scale']==2
    assert details['negative_fraction']>0


def test_descending_axis_processing_restores_all_stages_to_original_order():
    ascending = make_record()
    descending = make_record()
    descending.x = descending.x[::-1].copy()
    descending.channels = {name: values[::-1].copy() for name, values in descending.channels.items()}
    original_axis = descending.x.copy()
    original_signal = descending.channels['dark_subtracted'].copy()
    options = {'baseline': {'method': 'asls', 'lam': 1e6}, 'despike': {'enabled': False},
               'smooth': {'enabled': True, 'window_length': 7, 'polyorder': 2},
               'sensitivity': {'enabled': True, 'methods': ['asls', 'arpls', 'airpls'],
                               'lam_factors': [1], 'smooth_windows': [5, 7]}, 'normalization': 'none'}
    process_record(ascending, options)
    process_record(descending, options)
    np.testing.assert_array_equal(descending.x, original_axis)
    np.testing.assert_array_equal(descending.channels['dark_subtracted'], original_signal)
    assert abs(descending.x[np.argmax(descending.stages['unified'])] - 495) < 5
    assert ascending.stages.keys() == descending.stages.keys()
    for name in ascending.stages:
        np.testing.assert_allclose(descending.stages[name][::-1], ascending.stages[name], rtol=1e-10, atol=1e-10,
                                   err_msg=f'Axis alignment lost in stage {name}')


@pytest.mark.parametrize('field,first,second',[
    ('measurement_position','front','back'),('objective_magnification',50,10),
    ('laser_wavelength_nm',532,785),('averages',1,3),
])
def test_blank_subtraction_rejects_known_optical_and_average_conflicts(field,first,second):
    sample,blank=make_record(),make_record()
    sample.metadata[field]=first; blank.metadata[field]=second
    with pytest.raises(ValueError,match=field):
        subtract_background(sample,blank,{'conditions_verified':True})


def test_preprocessing_failure_does_not_publish_partial_or_misaligned_stages():
    rec=make_record()
    rec.x=np.arange(6.,-1.,-1.)
    rec.channels={'dark_subtracted':np.arange(10.,80.,10.),'baseline_subtracted':np.arange(7.)}
    rec.stages={'previous_valid':np.array([1.,2.])}
    rec.history=[{'operation':'import'}]
    rec.warnings=['previous_warning']
    options={'baseline':{'method':'none'},'despike':{'enabled':False},
             'smooth':{'enabled':True,'window_length':7,'polyorder':2},
             'sensitivity':{'enabled':True,'methods':[],'smooth_windows':[5,7,9]}}
    with pytest.raises(ValueError,match='sensitivity SavGol window: 9'):
        process_record(rec,options)
    assert set(rec.stages)=={'previous_valid'}
    np.testing.assert_array_equal(rec.stages['previous_valid'],[1.,2.])
    assert rec.history==[{'operation':'import'}]
    assert rec.warnings==['previous_warning']


def test_successful_reprocessing_does_not_retain_stale_sensitivity_arrays():
    rec=make_record()
    rec.x=rec.x[::-1].copy()
    rec.channels={name:y[::-1].copy() for name,y in rec.channels.items()}
    options={'baseline':{'method':'none'},'despike':{'enabled':False},'smooth':{'enabled':False},
             'sensitivity':{'enabled':True,'methods':['asls'],'lam_factors':[1],'smooth_windows':[]}}
    process_record(rec,options)
    assert 'asls_lam_1' in rec.stages
    options['sensitivity']['enabled']=False
    process_record(rec,options)
    assert 'asls_lam_1' not in rec.stages
    np.testing.assert_array_equal(rec.stages['unified'],rec.channels['dark_subtracted'])
    assert rec.history[-2]['parameters']['sensitivity']['enabled'] is True
