from types import SimpleNamespace
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from optosky_raman.config import merge_config
from optosky_raman.pipeline import _multivariate


def test_ramanspy_pca_scores_follow_record_order(tmp_path):
    matrix=np.random.default_rng(4).uniform(size=(4,20)); x=np.arange(20.)
    records=[SimpleNamespace(id=str(i),x=x,stages={'unified':y},state='applied_potential') for i,y in enumerate(matrix)]
    (tmp_path/'tables').mkdir()
    cfg=merge_config({'multivariate':{'enabled':True,'method':'pca','normalization':'none','n_components':2}})
    out=_multivariate(records,cfg,tmp_path)
    scores=pd.read_csv(tmp_path/'tables/multivariate_scores.csv')
    truth=PCA(n_components=2).fit_transform(matrix)
    assert out['status']=='computed'
    np.testing.assert_allclose(scores[['component_1','component_2']].to_numpy(),truth,rtol=1e-9,atol=1e-9)


def test_nmf_does_not_clip_negative_data(tmp_path):
    records=[SimpleNamespace(id=str(i),x=np.arange(20.),stages={'unified':np.arange(20.)-5},state='applied_potential') for i in range(4)]
    cfg=merge_config({'multivariate':{'enabled':True,'method':'nmf','normalization':'none'}})
    out=_multivariate(records,cfg,tmp_path)
    assert out['status']=='skipped'
    assert 'negative' in out['reason']


def test_ramanspy_nmf_valid_nonnegative_input_runs(tmp_path):
    x=np.arange(30.)
    basis=np.stack([np.exp(-((x-8)/3)**2),np.exp(-((x-21)/4)**2)])
    weights=np.random.default_rng(2).uniform(.2,2,(6,2))
    records=[SimpleNamespace(id=str(i),x=x,stages={'unified':y},state='applied_potential') for i,y in enumerate(weights@basis)]
    (tmp_path/'tables').mkdir()
    cfg=merge_config({'multivariate':{'enabled':True,'method':'nmf','normalization':'none'}})
    out=_multivariate(records,cfg,tmp_path)
    assert out['status']=='computed'
    scores=pd.read_csv(tmp_path/'tables/multivariate_scores.csv').iloc[:,1:].to_numpy()
    components=pd.read_csv(tmp_path/'tables/multivariate_components.csv').iloc[:,:2].to_numpy().T
    assert (scores>=0).all()
    np.testing.assert_allclose(scores@components,weights@basis,atol=.02)


def test_randomized_pca_uses_configured_seed(tmp_path):
    matrix=np.random.default_rng(120).normal(size=(40,600)); x=np.arange(600.)
    records=[SimpleNamespace(id=str(i),x=x,stages={'unified':y},state='applied_potential') for i,y in enumerate(matrix)]
    cfg=merge_config({'multivariate':{'enabled':True,'method':'pca','normalization':'none','seed':42}})
    (tmp_path/'tables').mkdir()
    _multivariate(records,cfg,tmp_path)
    before=(tmp_path/'tables/multivariate_scores.csv').read_bytes()
    _multivariate(records,cfg,tmp_path)
    assert (tmp_path/'tables/multivariate_scores.csv').read_bytes()==before
