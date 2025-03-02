'''
Module with tests for MCParPdf class
'''
import os
import copy

import ROOT
import zfit
import pytest
from ROOT                  import RDataFrame
from dmu.logging.log_store import LogStore
from rx_data.rdf_getter    import RDFGetter
from rx_fitter.mc_par_pdf  import MCParPdf
from rx_fitter             import components as cmp

log = LogStore.add_logger('rx_fitter:test_mc_par_pdf')
# ------------------------------------------
class Data:
    '''
    data class
    '''
    os.environ['FITDIR'] = '/tmp/tests/rx_fitter/mc_par_pdf/fits'
    cache_dir            = '/tmp/tests/rx_fitter/mc_par_pdf'

    cfg = {
            'fitting':
            {
                'error_method'  : 'minuit_hesse',
                'weights_column': 'weights',
                'ntries'        : 3,
                'pvalue'        : 0.02,
                },
            'plotting' :
            {
                'nbins'   : 50,
                'stacked' : True,
                },
            }

    obs = zfit.Space('B_M', limits=[4500, 6000])

    RDFGetter.samples = {
        'main'       : '/home/acampove/external_ssd/Data/samples/main.yaml',
        'mva'        : '/home/acampove/external_ssd/Data/samples/mva.yaml',
        'hop'        : '/home/acampove/external_ssd/Data/samples/hop.yaml',
        'cascade'    : '/home/acampove/external_ssd/Data/samples/cascade.yaml',
        'jpsi_misid' : '/home/acampove/external_ssd/Data/samples/jpsi_misid.yaml'}
# ------------------------------------------
@pytest.fixture(scope='session', autouse=True)
def _initialize():
    LogStore.set_level('rx_fitter:mc_par_pdf', 10)
# ------------------------------------------
def _get_rdf(cfg : dict) -> RDataFrame:
    nbrem   = cfg['nbrem']
    cuts    = {'nbrem' : f'nbrem == {nbrem}'}
    trigger = cfg['trigger']
    q2bin   = cfg['q2bin']
    sample  = cfg['name']

    return cmp.get_rdf(sample, q2bin, trigger, cuts)
# ------------------------------------------
def test_read():
    '''
    Used to read inputs
    '''
    cfg            = copy.deepcopy(Data.cfg)
    cfg['name'   ] = 'Bu_JpsiK_ee_eq_DPC'
    cfg['q2bin'  ] = 'jpsi'
    cfg['trigger'] = 'Hlt2RD_BuToKpEE_MVA'
    cfg['nbrem'  ] = 1
    cfg['fvers'  ] = None
    cfg['create' ] = False
    cfg['shared' ] = ['mu']
    cfg['model'  ] = ['suj']
    cfg['pfloat' ] = ['mu', 'sg']

    rdf = _get_rdf(cfg=cfg)
    obj = MCParPdf(rdf=rdf, obs=Data.obs, cfg=cfg)
    fcm = obj.get_fcomp()

    fcm.run()
# ------------------------------------------
def test_create():
    '''
    Used to create a new version
    '''
    cfg            = copy.deepcopy(Data.cfg)
    cfg['name'   ] = 'Bu_JpsiK_ee_eq_DPC'
    cfg['q2bin'  ] = 'jpsi'
    cfg['trigger'] = 'Hlt2RD_BuToKpEE_MVA'
    cfg['nbrem'  ] = 1
    cfg['fvers'  ] = None
    cfg['create' ] = True
    cfg['shared' ] = ['mu']
    cfg['model'  ] = ['suj']
    cfg['pfloat' ] = ['mu', 'sg']

    rdf = _get_rdf(cfg=cfg)
    obj = MCParPdf(rdf=rdf, obs=Data.obs, cfg=cfg)
    fcm = obj.get_fcomp()

    fcm.run()
# ------------------------------------------
