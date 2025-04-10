'''
Module with functions testing SignalScales class
'''

import pytest
from dmu.logging.log_store   import LogStore
from rx_fitter.signal_scales import FitParameters

log=LogStore.add_logger('rx_fitter:test_signal_scales')
# --------------------------------------------------------------
class Data:
    '''
    Data class
    '''
    l_sig_par = [
            'ar_dscb_Signal_002_1_reso_flt',
            'mu_Signal_000_scale_flt',
            'mu_Signal_001_scale_flt',
            'mu_Signal_002_scale_flt',
            'nl_dscb_Signal_001_1_reso_flt',
            'nr_dscb_Signal_002_1_reso_flt',
            'sg_Signal_000_reso_flt',
            'sg_Signal_001_reso_flt',
            'sg_Signal_002_reso_flt',
            ]

    l_brem_frac = [
            'frac_brem_000',
            'frac_brem_001',
            'frac_brem_002',
            ]

    l_invalid = [
            'sBd_Kstee_eq_btosllball05_DPC',
            'sBu_Kstee_Kpi0_eq_btosllball05_DPC',
            'ap_hypexp',
            'bt_hypexp',
            'mu_hypexp',
            'ncmb',
            'nsig',
            ]
# --------------------------------------------------------------
@pytest.fixture(scope='session', autouse=True)
def _intialize():
    LogStore.set_level('rx_fitter:signal_scales', 10)
# ------------------------------------
def test_get_data():
    '''
    Tests getting dataframe with parameters
    '''
    obj = FitParameters()
    df  = obj.get_data()

    log.info(df)

    assert len(df) > 0
# ------------------------------------
@pytest.mark.parametrize('name', Data.l_sig_par)
def test_get_mass_scales(name : str):
    '''
    Tests mass scales and resolutions
    '''
    obj      = FitParameters()
    val, err = obj.get_parameter_scale(name=name)

    log.info(f'Value: {val:.3f}')
    log.info(f'Error: {err:.3f}')
# ------------------------------------
@pytest.mark.parametrize('name', Data.l_brem_frac)
def test_get_brem_values(name : str):
    '''
    Tests retrieval of data brem fractions
    '''
    obj      = FitParameters()
    val, err = obj.get_brem_fraction(name=name)

    log.info(f'Value: {val:.3f}')
    log.info(f'Error: {err:.3f}')
