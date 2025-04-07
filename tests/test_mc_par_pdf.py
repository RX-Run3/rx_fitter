'''
Module with tests for MCParPdf class
'''
import os
import shutil
from importlib.resources import files

import yaml
import ROOT
import zfit
import numpy
import pytest
from ROOT                  import RDataFrame, RDF
from dmu.logging.log_store import LogStore
from dmu.stats.utilities   import print_pdf
from rx_calibration.hltcalibration.fit_component import NoFitDataFound
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
    arr_mass = numpy.random.normal(loc=5280, scale=200, size=1_000)
    arr_brem = numpy.random.randint(0, 3, size=1_000)

    rdf = RDF.FromNumpy({'B_M' : arr_mass, 'nbrem' : arr_brem})

    return rdf
# ------------------------------------------
def test_read_fail():
    '''
    Will require failure due to missing RDF and fit parameters
    '''
    cfg = _load_config('read')
    out_dir = cfg['output']['out_dir']

    if os.path.isdir(out_dir):
        shutil.rmtree(out_dir) # otherwise old fit will allow loading of parameters, and test won't fail

    obj = MCParPdf(rdf=None, obs=Data.obs, cfg=cfg)

    with pytest.raises(NoFitDataFound):
        obj.get_pdf(must_load_pars=True)
# ------------------------------------------
def test_create():
    '''
    Used to create a new version
    '''
    cfg = _load_config('read')
    rdf = _get_rdf()
    obj = MCParPdf(rdf=rdf, obs=Data.obs, cfg=cfg)
    pdf = obj.get_pdf()

    print_pdf(pdf)
# ------------------------------------------
def test_read():
    '''
    Used to read input parameters
    '''
    cfg = _load_config('read')
    obj = MCParPdf(rdf=None, obs=Data.obs, cfg=cfg)
    pdf = obj.get_pdf(must_load_pars=True)

    print_pdf(pdf)
# ------------------------------------------
def test_read_reparametrize():
    '''
    Used to read inputs
    '''
    cfg = _load_config('read_reparametrize')
    obj = MCParPdf(rdf=None, obs=Data.obs, cfg=cfg)
    pdf = obj.get_pdf()

    print_pdf(pdf)
# ------------------------------------------
