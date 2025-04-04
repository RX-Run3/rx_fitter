'''
Module used to test PrecScales class
'''

import os
import pytest
import matplotlib.pyplot as plt

from conftest                    import ScalesData
from dmu.logging.log_store       import LogStore
from rx_efficiencies.decay_names import DecayNames as dn
from rx_fitter.prec_scales       import PrecScales

log=LogStore.add_logger('rx_fitter:test_prec_scales')
# -----------------------------------
@pytest.fixture(scope='session', autouse=True)
def _initialize():
    LogStore.set_level('rx_fitter:prec_scales', 10)
#-------------------------------
def _plot_df(df, trig):
    df = df[df.trig == trig]
    ax = None
    for proc, df_p in df.groupby('kind'):
        label = dn.tex_from_decay(proc)
        ax=df_p.plot(x='year', y='val', yerr='err', ax=ax, label=label, linestyle='none', marker='o')

    os.makedirs('tests/scales/all_datasets', exist_ok=True)

    plt_path = f'tests/scales/all_datasets/{trig}.png'
    log.info(f'Saving to: {plt_path}')
    plt.grid()
    plt.ylabel('Scale')
    plt.ylim(0.0, 0.5)
    plt.xlabel('')
    plt.tight_layout()
    plt.savefig(plt_path)
    plt.close('all')
#-------------------------------
@pytest.mark.parametrize('q2bin'  , ['low', 'central', 'high'])
@pytest.mark.parametrize('process', dn.get_decays())
def test_all_datasets(q2bin : str, process : str):
    '''
    Tests retrieval of scales between signal and PRec yields
    '''
    signal   = 'bpkpee'

    obj      = PrecScales(proc=process, q2bin=q2bin)
    val, err = obj.get_scale(signal=signal)

    if process != signal:
        assert val  < 1   # Prec should be smaller than signal
    else:
        assert val == 1   # If this runs on signal, scale is 1

    if val > 0:          # If no events pass selection, error will be zero, as well as value
        assert err < val # Error smaller than value

    ScalesData.collect(process, q2bin, val, err)
#-------------------------------
