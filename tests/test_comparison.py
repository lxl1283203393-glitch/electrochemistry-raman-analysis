from types import SimpleNamespace
from datetime import datetime, timedelta
import numpy as np
from optosky_raman.comparison import compare_records


def test_multiplicative_drop_and_time_potential_confounding():
    x=np.linspace(200,800,301); y=10+np.exp(-((x-500)/40)**2)*100
    records=[]
    for i,s in enumerate([1,0.15,0.8]):
        records.append(SimpleNamespace(id=str(i),batch='b',path=SimpleNamespace(name=f'{i}.csv'),
            x=x,channels={'dark_subtracted':y*s},stages={'unified':y*s},
            state='applied_potential',potential_v=-1.1-i*.01,
            acquired_at=datetime(2026,10,12)+timedelta(minutes=i),
            metadata={'integration_time_s':10,'laser_power_mw':30},warnings=[]))
    out=compare_records(records,[],{})
    assert out['adjacent'].iloc[0]['scale_ratio']<0.2
    assert bool(out['adjacent'].iloc[0]['common_scale_suspected'])
    assert out['confounding']['b']['time_potential_confounded']
    assert len(out['difference_spectra'])==602


def test_insufficient_axis_overlap_returns_status_without_aborting_comparisons():
    records=[]
    for i,x in enumerate([np.arange(5.),np.arange(3.,8.)]):
        records.append(SimpleNamespace(id=str(i),batch='b',path=SimpleNamespace(name=f'{i}.csv'),
            x=x,channels={'dark_subtracted':np.arange(1.,6.)},stages={'unified':np.arange(1.,6.)},
            state='applied_potential',potential_v=-1.1-i*.01,
            acquired_at=datetime(2026,10,12)+timedelta(minutes=i),
            metadata={'integration_time_s':10,'laser_power_mw':30},warnings=[]))
    out=compare_records(records,[],{})
    row=out['adjacent'].iloc[0]
    assert row['status']=='insufficient_overlap'
    assert not row['common_scale_suspected']
    assert not row['candidate_change']
    assert out['difference_spectra'].empty
