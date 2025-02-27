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
@pytest.mark.parametrize('nbrem' , [0, 1, 2])
@pytest.mark.parametrize('mass'  , ['B_const_mass_M', 'B_M'])
@pytest.mark.parametrize('sample', ['Bu_JpsiK_ee_eq_DPC', 'Bu_JpsiPi_ee_eq_DPC'])
def test_signal(nbrem : int, mass : str, sample : str):
    '''
    Testing creation of PDF from MC sample
    '''
    limits = [5100, 5500] if mass == 'B_const_mass_M' else [4500, 6000]
    obs=zfit.Space(mass, limits=limits)
    trigger = 'Hlt2RD_BuToKpEE_MVA'

    cmp_sig = cmp.get_mc(obs    = obs,
                         sample = sample,
                         trigger= trigger,
                         q2bin  = 'jpsi',
                         nbrem  = nbrem)
    cmp_sig.run()
# --------------------------------------------------------------
@pytest.mark.parametrize('nbrem', [0, 1, 2])
@pytest.mark.parametrize('mass' , ['B_const_mass_M'])
def test_prec(nbrem : int, mass : str):
    '''
    Testing creation of PDF from MC sample
    '''
    cmp.Data.cfg['out_dir'] = f'/tmp/tests/rx_fitter/components/prec/{mass}'

    obs     = zfit.Space(mass, limits=(4500, 6000))
    trigger = 'Hlt2RD_BuToKpEE_MVA'
    cmp_prc = cmp.get_prc(obs    = obs,
                         q2bin  = 'jpsi',
                         trigger= trigger,
                         nbrem  = nbrem)
    cmp_prc.run()
# --------------------------------------------------------------
