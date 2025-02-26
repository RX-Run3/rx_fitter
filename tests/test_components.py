'''
Module with functions to test functions in components.py
'''
import ROOT
import zfit
import pytest

from dmu.logging.log_store  import LogStore
from rx_data.rdf_getter     import RDFGetter
from rx_fitter              import components as cmp

# --------------------------------------------------------------
@pytest.fixture(scope='session', autouse=True)
def _intiailize():
    LogStore.set_level('rx_fitter:prec', 10)

    RDFGetter.samples = {
        'main'       : '/home/acampove/external_ssd/Data/samples/main.yaml',
        'mva'        : '/home/acampove/external_ssd/Data/samples/mva.yaml',
        'hop'        : '/home/acampove/external_ssd/Data/samples/hop.yaml',
        'cascade'    : '/home/acampove/external_ssd/Data/samples/cascade.yaml',
        'jpsi_misid' : '/home/acampove/external_ssd/Data/samples/jpsi_misid.yaml'}
# --------------------------------------------------------------
@pytest.mark.parametrize('nbrem', [0, 1, 2])
def test_signal(nbrem : int):
    '''
    Testing creation of PDF from MC sample
    '''
    #obs=zfit.Space('B_const_mass_M', limits=(5100, 5500))
    obs=zfit.Space('B_M', limits=(4500, 6000))
    trigger = 'Hlt2RD_BuToKpEE_MVA'

    cmp_sig = cmp.get_mc(obs    = obs,
                         sample = 'Bu_JpsiK_ee_eq_DPC',
                         trigger= trigger,
                         q2bin  = 'jpsi',
                         nbrem  = nbrem)
    cmp_sig.run()
# --------------------------------------------------------------
@pytest.mark.parametrize('nbrem', [-1])
def test_prec(nbrem : int):
    '''
    Testing creation of PDF from MC sample
    '''
    cmp.Data.cfg['out_dir'] = '/tmp/tests/rx_fitter/components/prec'

    obs     = zfit.Space('B_M', limits=(4500, 6000))
    trigger = 'Hlt2RD_BuToKpEE_MVA'
    cmp_prc = cmp.get_prc(obs    = obs,
                         q2bin  = 'jpsi',
                         trigger= trigger,
                         nbrem  = nbrem)

    cmp_prc.run()
# --------------------------------------------------------------
