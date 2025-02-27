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
def _get_mass_range(mass : str, sample : str) -> list[int]:
    if mass == 'B_const_mass_M' and sample == 'Bu_JpsiPi_ee_eq_DPC':
        return [5100, 6000]

    if mass == 'B_const_mass_M' and sample == 'Bu_JpsiK_ee_eq_DPC':
        return [5100, 5500]

    if mass == 'B_M':
        return [4500, 6000]

    raise ValueError(f'Invalid mass and sample: {mass}/{sample}')
# --------------------------------------------------------------
@pytest.mark.parametrize('nbrem' , [0, 1, 2])
@pytest.mark.parametrize('mass'  , ['B_const_mass_M', 'B_M'])
@pytest.mark.parametrize('sample', ['Bu_JpsiPi_ee_eq_DPC'])
def test_signal(nbrem : int, mass : str, sample : str):
    '''
    Testing creation of PDF from MC sample
    '''
    limits = _get_mass_range(mass, sample)

    obs=zfit.Space(mass, limits=limits)
    trigger = 'Hlt2RD_BuToKpEE_MVA'

    cmp.Data.cfg['fitting']['ntries'] = 15

    cmp_sig = cmp.get_mc(obs    = obs,
                         sample = sample,
                         trigger= trigger,
                         q2bin  = 'jpsi',
                         nbrem  = nbrem)
    cmp_sig.run()
# --------------------------------------------------------------
@pytest.mark.parametrize('mass'     , ['B_const_mass_M', 'B_M'])
@pytest.mark.parametrize('cut, name', [
    ('nbrem == 0', 'bz'),
    ('nbrem == 1', 'bo'),
    ('nbrem >= 2', 'bt')])
def test_prec_brem(mass : str, cut : str, name : str):
    '''
    Testing creation of PDF from MC sample with brem cut
    '''
    cmp.Data.cfg['out_dir'] = f'/tmp/tests/rx_fitter/components/prec/{mass}_{name}'

    obs     = zfit.Space(mass, limits=(4500, 6000))
    trigger = 'Hlt2RD_BuToKpEE_MVA'
    cmp_prc = cmp.get_prc(
            name   = name,
            obs    = obs,
            q2bin  = 'jpsi',
            trigger= trigger,
            cuts   = {'nbrem' : cut, 'core' : 'B_const_mass_M > 5150'})

    cmp_prc.run()
# --------------------------------------------------------------
def test_combinatorial():
    '''
    Testing creation of PDF used for combinatorial
    '''

    obs=zfit.Space('B_M', limits=[4500, 6000])
    cmp_sig = cmp.get_cb(obs=obs, kind='exp')
    cmp_sig.run()
# --------------------------------------------------------------
