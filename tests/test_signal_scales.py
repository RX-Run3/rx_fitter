'''
Module with functions testing SignalScales class
'''

import pytest
from dmu.logging.log_store   import LogStore
from rx_fitter.signal_scales import SignalScales

log=LogStore.add_logger('rx_fitter:test_signal_scales')
# --------------------------------------------------------------
@pytest.fixture(scope='session', autouse=True)
def _intialize():
    LogStore.set_level('rx_fitter:signal_scales', 10)
# ------------------------------------
def test_get_data():
    '''
    Tests getting dataframe with parameters
    '''
    obj = SignalScales()
    df  = obj.get_data()

    log.info(df)

    assert len(df) > 0
# ------------------------------------
