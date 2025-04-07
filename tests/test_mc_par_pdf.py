'''
Module with tests for MCParPdf class
'''
import copy

import ROOT
import zfit
import numpy
import pytest
from ROOT                  import RDataFrame, RDF
from dmu.logging.log_store import LogStore
from rx_fitter.mc_par_pdf  import MCParPdf

log = LogStore.add_logger('rx_fitter:test_mc_par_pdf')
# ------------------------------------------
class Data:
    '''
    data class
    '''
    cfg = {
            'component_name'   : 'Bu_JpsiK_ee_eq_DPC',
            'fvers'  : None,
            'create' : False,
            'q2bin'  : 'jpsi',
            'nbrem'  : 1,
            'trigger': 'Hlt2RD_BuToKpEE_MVA',
            'shared' : ['mu'],
            'model'  : ['cbl'],
            'pfloat' : ['mu', 'sg'],
            'output' :
            {
                'out_dir' : '/tmp/tests/rx_fitter/mc_par_pdf/fits',
                },
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
# ------------------------------------------
@pytest.fixture(scope='session', autouse=True)
def _initialize():
    LogStore.set_level('rx_fitter:mc_par_pdf', 10)
# ------------------------------------------
def _get_rdf() -> RDataFrame:
    arr_mass = numpy.random.normal(loc=5280, scale=100, size=1_000)
    arr_brem = numpy.random.randint(0, 3, size=1_000)

    rdf = RDF.FromNumpy({'B_M' : arr_mass, 'nbrem' : arr_brem})

    return rdf
# ------------------------------------------
def test_read():
    '''
    Used to read inputs
    '''
    cfg = copy.deepcopy(Data.cfg)
    rdf = _get_rdf()
    obj = MCParPdf(rdf=rdf, obs=Data.obs, cfg=cfg)
    fcm = obj.get_fcomp()

    fcm.run()
# ------------------------------------------
def test_create():
    '''
    Used to create a new version
    '''
    cfg            = copy.deepcopy(Data.cfg)
    cfg['create' ] = True

    rdf = _get_rdf()
    obj = MCParPdf(rdf=rdf, obs=Data.obs, cfg=cfg)
    fcm = obj.get_fcomp()

    fcm.run()
# ------------------------------------------
def test_fix_pars():
    '''
    Used to create a new version with parameters fixed from old version
    '''
    cfg            = copy.deepcopy(Data.cfg)
    cfg['fvers'  ] = 'v2'
    cfg['create' ] = True

    rdf = _get_rdf()
    obj = MCParPdf(rdf=rdf, obs=Data.obs, cfg=cfg)
    fcm = obj.get_fcomp()

    fcm.run()
# ------------------------------------------
