'''
Module with functions to test functions in components.py
'''
import copy

import ROOT
import zfit
import pytest

from dmu.logging.log_store  import LogStore
from rx_data.rdf_getter     import RDFGetter
from rx_fitter              import components as cmp

# --------------------------------------------------------------
class Data:
    '''
    Data class
    '''
    cfg = {
            'input': {
                'q2bin': 'jpsi',
                'trigger': 'Hlt2RD_BuToKpEE_MVA',
                'samples': {
                    'main': '/home/acampove/external_ssd/Data/samples/main.yaml',
                    'mva': '/home/acampove/external_ssd/Data/samples/mva.yaml',
                    'hop': '/home/acampove/external_ssd/Data/samples/hop.yaml',
                    'cascade': '/home/acampove/external_ssd/Data/samples/cascade.yaml',
                    'jpsi_misid': '/home/acampove/external_ssd/Data/samples/jpsi_misid.yaml'
                    },
                'selection': {
                    'mass': 'B_const_mass_M > 5160'
                    }
            },
            'output': {
                'fit_dir': '/tmp/tests/rx_fitter/components',
            },
            'fitting': {
                'range': {
                    'B_M': [
                        4500,
                        6000
                        ],
                    'B_const_mass_M': [
                        5160,
                        5500
                        ]
                    },
                'components': {
                    'Signal': true,
                    'Cabibbo': false,
                    'PRec': false,
                    'combinatorial': false,
                    'data': false
                    },
                'config': {
                    'data': {
                        'fitting': {
                            'error_method': 'minuit_hesse'
                            },
                        'plotting': {
                            'nbins': 30,
                            'stacked': true,
                            'd_leg': {
                                'Bu_JpsiK_ee_eq_DPC': '$B^+\\to K^+J/\\psi(\\to e^+e^-)$',
                                'Bu_JpsiPi_ee_eq_DPC': '$B^+\\to \\pi^+J/\\psi(\\to e^+e^-)$',
                                'combinatorial': 'Combinatorial'
                                }
                            }
                        },
                    'Signal': {
                        'sample': 'Bu_JpsiK_ee_eq_DPC',
                        'fitting': {
                            'error_method': 'minuit_hesse',
                            'weights_column': 'weights',
                            'ntries': 20,
                            'pvalue': 0.02
                            },
                        'plotting': {
                            'nbins': 30,
                            'stacked': true
                            }
                        },
                    'Cabibbo': {
                        'sample': 'Bu_JpsiPi_ee_eq_DPC',
                        'fitting': {
                            'error_method': 'minuit_hesse',
                            'weights_column': 'weights',
                            'ntries': 20,
                            'pvalue': 0.02
                            },
                        'plotting': {
                            'nbins': 30,
                            'stacked': true
                            }
                        },
                    'PRec': {
                        'bw': 20
                        },
                    'combinatorial': {
                        'kind': 'exp'
                        }
                    }
            },
            'brem': {
                    0 : 'nbrem == 0',
                    1 : 'nbrem == 1',
                    2 : 'nbrem >= 2'
            },
            'components': {
                    'Signal': {
                        '0': {
                            'model': [
                                'suj',
                                'suj'
                                ],
                            'pfloat': [
                                'mu',
                                'sg'
                                ],
                            'shared': [
                                'mu'
                                ],
                            'fvers' : 'v2',
                            'create': True,
                            },
                        '1': {
                            'model': [
                                'suj',
                                'dscb'
                                ],
                            'pfloat': [
                                'mu',
                                'sg'
                                ],
                            'shared': [
                                'mu'
                                ],
                            'fvers' : 'v2',
                            'create': True,
                            },
                        '2': {
                            'model': [
                                'suj',
                                'dscb'
                                ],
                            'pfloat': [
                                'mu',
                                'sg'
                                ],
                            'shared': [
                                'mu'
                                ],
                            'fvers' : 'v2',
                            'create': True,
                            }
                        },
                    'Cabibbo': {
                        '0': {
                            'model': [
                                'suj'
                                ],
                            'pfloat': [],
                            'shared': []
                            },
                        '1': {
                            'model': [
                                'suj'
                                ],
                            'pfloat': [],
                            'shared': []
                            },
                        '2': {
                            'model': [
                                'suj'
                                ],
                            'pfloat': [],
                            'shared': []
                            }
                        }
            }
  }
# --------------------------------------------------------------
@pytest.fixture(scope='session', autouse=True)
def _intiailize():
    LogStore.set_level('rx_fitter:prec'              , 10)
    LogStore.set_level('rx_fitter:components'        , 10)
    LogStore.set_level('rx_calibration:fit_component', 10)

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
@pytest.mark.parametrize('mass'  , ['B_M'])
@pytest.mark.parametrize('sample', ['Bu_JpsiK_ee_eq_DPC'])
def test_signal(nbrem : int, mass : str, sample : str):
    '''
    Testing creation of PDF from MC sample
    '''
    limits = _get_mass_range(mass, sample)

    obs=zfit.Space(mass, limits=limits)
    trigger = 'Hlt2RD_BuToKpEE_MVA'

    cmp.Data.cfg['fitting']['ntries'] = 15

    cmp.Data.cfg['out_dir'] = '/tmp/tests/rx_fitter/components/signal'
    cmp_sig = cmp.get_mc(obs    = obs,
                         sample = sample,
                         trigger= trigger,
                         q2bin  = 'jpsi',
                         model  = ['cbl', 'cbr'],
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
    cmp.Data.cfg['out_dir'] = '/tmp/tests/rx_fitter/components/prec'

    obs     = zfit.Space(mass, limits=(4500, 6000))
    trigger = 'Hlt2RD_BuToKpEE_MVA'
    cmp_prc = cmp.get_prc(
            name   = name,
            obs    = obs,
            q2bin  = 'jpsi',
            trigger= trigger,
            bw     = 20,
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
