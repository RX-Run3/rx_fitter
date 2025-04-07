'''
Module with tests for MCParPdf class
'''
import copy
from importlib.resources import files

import yaml
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
    cfg : dict
    obs = zfit.Space('B_M', limits=[4500, 6000])
# ------------------------------------------
def _load_config(name : str) -> None:
    cfg_path = files('rx_fitter_data').joinpath(f'tests/mc_par_pdf/{name}.yaml')

    with open(cfg_path, encoding='utf-8') as ifile:
        data = yaml.safe_load(ifile)

    return data
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
    cfg = _load_config('read')
    rdf = _get_rdf()
    obj = MCParPdf(rdf=rdf, obs=Data.obs, cfg=cfg)
    fcm = obj.get_fcomp()

    fcm.run()
# ------------------------------------------
def test_create():
    '''
    Used to create a new version
    '''
    cfg            = _load_config('read')
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
    cfg            = _load_config('read')
    cfg['fvers'  ] = 'v2'
    cfg['create' ] = True

    rdf = _get_rdf()
    obj = MCParPdf(rdf=rdf, obs=Data.obs, cfg=cfg)
    fcm = obj.get_fcomp()

    fcm.run()
# ------------------------------------------
