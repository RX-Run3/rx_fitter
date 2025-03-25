'''
Module with functions to test functions in components.py
'''
import copy

import ROOT
import zfit
import pytest

from zfit.core.interfaces   import ZfitSpace as zobs
from dmu.logging.log_store  import LogStore
from rx_fitter              import components as cmp

log=LogStore.add_logger('rx_fitter:test_components')
# --------------------------------------------------------------
class Data:
    '''
    Data class
    '''
    cfg = {
            'input': {
                'q2bin'   : 'jpsi',
                'trigger' : 'Hlt2RD_BuToKpEE_MVA',
                'samples': {
                    'main'        : '/home/acampove/external_ssd/Data/samples/main.yaml',
                    'mva'         : '/home/acampove/external_ssd/Data/samples/mva.yaml',
                    'hop'         : '/home/acampove/external_ssd/Data/samples/hop.yaml',
                    'cascade'     : '/home/acampove/external_ssd/Data/samples/cascade.yaml',
                    'jpsi_misid'  : '/home/acampove/external_ssd/Data/samples/jpsi_misid.yaml',
                    'brem_track_2': '/home/acampove/external_ssd/Data/samples/brem_track_2.yaml'
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
                    'Signal': True,
                    'Cabibbo': False,
                    'PRec': False,
                    'combinatorial': False,
                    'data': False
                    },
                'config': {
                    'data': {
                        'fitting': {
                            'error_method': 'minuit_hesse'
                            },
                        'plotting': {
                            'nbins': 30,
                            'stacked': True,
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
                            'ntries': 2,
                            'pvalue': 0.02
                            },
                        'plotting': {
                            'nbins': 30,
                            'stacked': True
                            }
                        },
                    'Cabibbo': {
                        'sample': 'Bu_JpsiPi_ee_eq_DPC',
                        'fitting': {
                            'error_method': 'minuit_hesse',
                            'ntries': 2,
                            'pvalue': 0.02
                            },
                        'plotting': {
                            'nbins': 30,
                            'stacked': True
                            }
                        },
                    'PRec': {
                        'cfg_kde':
                        {
                        'bandwidth': 20,
                        'padding'  : {'lowermirror': 0.5, 'uppermirror': 0.5},
                            },
                        'sample' : [
                            'Bu_JpsiX_ee_eq_JpsiInAcc',
                            'Bd_JpsiX_ee_eq_JpsiInAcc',
                            'Bs_JpsiX_ee_eq_JpsiInAcc',
                            ],
                        'weights' : {
                            'dec' : 1,
                            'sam' : 1,
                            },
                        'plotting' : {
                                'nbins'   : 30,
                                'stacked' : True,
                                },
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
                        0: {
                            'model' : ['cbl'],
                            'pfloat': [
                                'mu',
                                'sg'
                                ],
                            'shared': [
                                'mu'
                                ],
                            'fvers' : 'v2',
                            'create': True,
                            'weights' : 'weights',
                            },
                        1: {
                            'model': ['dscb'],
                            'pfloat': [
                                'mu',
                                'sg'
                                ],
                            'shared': [
                                'mu'
                                ],
                            'fvers' : 'v2',
                            'create': True,
                            'weights' : 'weights',
                            },
                        2: {
                            'model' : ['dscb'],
                            'pfloat': [
                                'mu',
                                'sg'
                                ],
                            'shared': [
                                'mu'
                                ],
                            'fvers' : 'v2',
                            'create': True,
                            'weights' : 'weights',
                            }
                        },
                    'Cabibbo': {
                        0: {
                            'model': [
                                'suj'
                                ],
                            'pfloat': [],
                            'shared': []
                            },
                        1: {
                            'model': [
                                'suj'
                                ],
                            'pfloat': [],
                            'shared': []
                            },
                        2: {
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
# --------------------------------------------------------------
def _get_obs(mass : str, cfg : dict) -> zobs:
    [min_mass, max_mass] = cfg['fitting']['range'][mass]

    obs = zfit.Space(mass, limits=(min_mass, max_mass))

    return obs
# --------------------------------------------------------------
@pytest.mark.parametrize('nbrem', [0, 1, 2])
@pytest.mark.parametrize('mass' , ['B_M'])
@pytest.mark.parametrize('name' , ['Signal'])
def test_mc_reuse(nbrem : int, mass : str, name : str):
    '''
    Testing reuse of old fit
    '''
    cfg            = copy.deepcopy(Data.cfg)
    cfg['out_dir'] = f'/tmp/tests/rx_fitter/components/test_mc/{name}_{mass}_{nbrem:03}'
    d_cmp_set      = cfg['components'][name][nbrem]
    d_cmp_set['create'] = False
    d_cmp_set['fvers' ] = None

    obs            = _get_obs(mass, cfg)

    cmp_sig = cmp.get_mc(obs=obs, component_name=name, nbrem=nbrem, cfg=cfg)
    cmp_sig.run()
# --------------------------------------------------------------
@pytest.mark.parametrize('nbrem', [0, 1, 2])
@pytest.mark.parametrize('mass' , ['B_M'])
@pytest.mark.parametrize('name' , ['Signal'])
def test_mc_create(nbrem : int, mass : str, name : str):
    '''
    Testing creation of PDF from MC sample
    '''
    cfg            = copy.deepcopy(Data.cfg)
    cfg['out_dir'] = f'/tmp/tests/rx_fitter/components/test_mc/{name}_{mass}_{nbrem:03}'
    d_cmp_set      = cfg['components'][name][nbrem]
    d_cmp_set['fvers'] = None

    obs            = _get_obs(mass, cfg)

    cmp_sig = cmp.get_mc(obs=obs, component_name=name, nbrem=nbrem, cfg=cfg)
    cmp_sig.run()
# --------------------------------------------------------------
@pytest.mark.parametrize('nbrem', [0, 1, 2])
@pytest.mark.parametrize('mass' , ['B_M'])
@pytest.mark.parametrize('name' , ['Signal'])
def test_mc_fix(nbrem : int, mass : str, name : str):
    '''
    Testing creation of PDF from MC sample with tails fixed from other version
    '''
    cfg            = copy.deepcopy(Data.cfg)
    cfg['out_dir'] = f'/tmp/tests/rx_fitter/components/test_mc/{name}_{mass}_{nbrem:03}'
    d_cmp_set      = cfg['components'][name][nbrem]
    d_cmp_set['fvers'] = 'v1'

    obs     = _get_obs(mass, cfg)
    cmp_sig = cmp.get_mc(obs=obs, component_name=name, nbrem=nbrem, cfg=cfg)
    cmp_sig.run()
# --------------------------------------------------------------
@pytest.mark.parametrize('nbrem',                 [0, 1, 2])
@pytest.mark.parametrize('mass' , ['B_const_mass_M', 'B_M'])
def test_prec_brem(mass : str, nbrem : int):
    '''
    Testing creation of PDF from MC sample with brem cut
    '''
    log.info('')

    cfg            = copy.deepcopy(Data.cfg)
    cfg['out_dir'] = f'/tmp/tests/rx_fitter/components/test_prec_brem/{mass}_{nbrem:03}/v1'

    obs            = _get_obs(mass, cfg)
    cmp_prc        = cmp.get_prc(obs, nbrem, cfg)
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
