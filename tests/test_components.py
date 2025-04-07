'''
Module with functions to test functions in components.py
'''
from importlib.resources import files

import ROOT
import zfit
import yaml
import pytest

from zfit.core.interfaces   import ZfitSpace as zobs
from dmu.logging.log_store  import LogStore
from dmu.stats.utilities    import print_pdf
from rx_fitter              import components as cmp

log=LogStore.add_logger('rx_fitter:test_components')
# --------------------------------------------------------------
class Data:
    '''
    Data class
    '''
    mass = 'B_M_brem_track_2'
# --------------------------------------------------------------
@pytest.fixture(scope='session', autouse=True)
def _intiailize():
    LogStore.set_level('rx_fitter:prec'              , 10)
    LogStore.set_level('rx_fitter:components'        , 10)
    LogStore.set_level('rx_fitter:mc_par_pdf'        , 10)
    LogStore.set_level('rx_calibration:fit_component', 10)
# --------------------------------------------------------------
def _load_config(test : str) -> dict:
    cfg_path = files('rx_fitter_data').joinpath(f'tests/components/{test}.yaml')
    with open(cfg_path, encoding='utf-8') as ifile:
        cfg = yaml.safe_load(ifile)

    return cfg
# --------------------------------------------------------------
def _get_obs(mass : str, cfg : dict) -> zobs:
    [min_mass, max_mass] = cfg['fitting']['range'][mass]

    obs = zfit.Space(mass, limits=(min_mass, max_mass))

    return obs
# --------------------------------------------------------------
def _get_brem_cut(nbrem : int, kind : str) -> str:
    l1 = 'L1_HASBREMADDED'
    l2 = 'L2_HASBREMADDED'
    if kind is not None:
        l1 = f'{l1}_{kind}'
        l2 = f'{l2}_{kind}'

    expr = f'int({l1}) + int({l2})'

    if nbrem in [0, 1]:
        cut = f'{expr} == {nbrem}'
    else:
        cut = f'{expr} >= {nbrem}'

    log.info(f'Using brem requirement: {cut}')

    return cut
# --------------------------------------------------------------
def _get_fitting_range(kind : str) -> dict[str:list[int]]:
    if kind is None:
        return {Data.mass : [4500, 6000]}

    return {f'B_M_{kind}' : [4500, 6000]}
# --------------------------------------------------------------
@pytest.mark.parametrize('nbrem', [0, 1, 2])
@pytest.mark.parametrize('kind' , ['brem_track_2'])
@pytest.mark.parametrize('mass' , ['B_M_brem_track_2'])
def test_brem_definitions(nbrem : int, kind : str, mass : str):
    '''
    Will test old and new brem definition
    '''
    log.info()

    cfg                      = _load_config('mc')
    out_dir                  = cfg['output']['out_dir']
    cfg['output']['out_dir'] = f'{out_dir}/test_brem_definitions'

    d_cmp_set                = cfg['components']['Signal'][nbrem]
    d_cmp_set['create']      = False
    d_cmp_set['fvers' ]      = None
    cfg['brem'][nbrem]       = _get_brem_cut(nbrem, kind)
    cfg['fitting']['range']  = _get_fitting_range(kind)

    obs     = _get_obs(mass, cfg)
    cmp_sig = cmp.get_mc(obs=obs, component_name='Signal', nbrem=nbrem, cfg=cfg)
    cmp_sig.run()
# --------------------------------------------------------------
@pytest.mark.parametrize('nbrem', [0, 1, 2])
@pytest.mark.parametrize('mass' , ['B_M_brem_track_2'])
@pytest.mark.parametrize('name' , ['Signal'])
def test_mc_create(nbrem : int, mass : str, name : str):
    '''
    Testing creation of PDF from MC sample
    '''
    log.info('')

    cfg                      = _load_config('mc_create')
    out_dir                  = cfg['output']['out_dir']
    cfg['output']['out_dir'] = f'{out_dir}/test_mc_create'

    d_cmp_set                = cfg['components'][name][nbrem]
    d_cmp_set['fvers']       = None

    obs     = _get_obs(mass, cfg)
    cmp_sig = cmp.get_mc(obs=obs, component_name=name, nbrem=nbrem, cfg=cfg)
    pdf     = cmp_sig.pdf

    print_pdf(pdf)
# --------------------------------------------------------------
@pytest.mark.parametrize('nbrem', [0, 1, 2])
@pytest.mark.parametrize('mass' , ['B_M_brem_track_2'])
@pytest.mark.parametrize('name' , ['Signal'])
def test_mc_reuse(nbrem : int, mass : str, name : str):
    '''
    Testing reuse of old fit
    '''
    log.info('')

    cfg                      = _load_config('mc_reuse')
    out_dir                  = cfg['output']['out_dir']
    cfg['output']['out_dir'] = f'{out_dir}/test_mc_create'

    d_cmp_set                = cfg['components'][name][nbrem]
    d_cmp_set['create']      = False
    d_cmp_set['fvers' ]      = None

    obs     = _get_obs(mass, cfg)
    cmp_sig = cmp.get_mc(obs=obs, component_name=name, nbrem=nbrem, cfg=cfg)
    pdf     = cmp_sig.pdf

    print_pdf(pdf)
# --------------------------------------------------------------
@pytest.mark.parametrize('nbrem', [0, 1, 2])
@pytest.mark.parametrize('mass' , ['B_M_brem_track_2'])
@pytest.mark.parametrize('name' , ['Signal'])
def test_mc_fix(nbrem : int, mass : str, name : str):
    '''
    Testing creation of PDF from MC sample with tails fixed from other version
    '''
    log.info('')

    cfg                      = _load_config('mc_reuse')
    out_dir                  = cfg['output']['out_dir']
    cfg['output']['out_dir'] = f'{out_dir}/test_mc_create'

    d_cmp_set                = cfg['components'][name][nbrem]
    d_cmp_set['fvers']       = 'v1'

    obs     = _get_obs(mass, cfg)
    cmp_sig = cmp.get_mc(obs=obs, component_name=name, nbrem=nbrem, cfg=cfg)
    pdf     = cmp_sig.pdf

    print_pdf(pdf)
# --------------------------------------------------------------
@pytest.mark.parametrize('nbrem', [0, 1, 2])
@pytest.mark.parametrize('mass' , ['B_M_brem_track_2'])
@pytest.mark.parametrize('name' , ['Signal'])
def test_mc_reparametrized(nbrem : int, mass : str, name : str):
    '''
    Testing creation of PDF from MC sample with tails fixed from other version
    '''
    log.info('')

    cfg                      = _load_config('mc_reparametrize')
    out_dir                  = cfg['output']['out_dir']
    cfg['output']['out_dir'] = f'{out_dir}/test_mc_create'

    d_cmp_set                = cfg['components'][name][nbrem]
    d_cmp_set['fvers']       = 'v1'

    obs  = _get_obs(mass, cfg)
    pdf  = cmp.get_mc_reparametrized(obs=obs, component_name=name, nbrem=nbrem, cfg=cfg)

    print_pdf(pdf)
# --------------------------------------------------------------
@pytest.mark.parametrize('mass' , ['B_M_brem_track_2'])
@pytest.mark.parametrize('name' , ['Signal'])
def test_mc_brem_reparametrized(mass : str, name : str):
    '''
    Test building full signal PDF for electron channel with Brem reparametrization
    '''
    log.info('')

    cfg                      = _load_config('mc_reparametrize')
    out_dir                  = cfg['output']['out_dir']
    cfg['output']['out_dir'] = f'{out_dir}/test_mc_create'

    for nbrem in [0, 1, 2]:
        d_cmp_set                = cfg['components'][name][nbrem]
        d_cmp_set['fvers']       = 'v1'

    obs  = _get_obs(mass, cfg)
    pdf  = cmp.get_mc_reparametrized(obs=obs, component_name=name, cfg=cfg, nbrem=None)

    print_pdf(pdf)
# --------------------------------------------------------------
@pytest.mark.parametrize('nbrem',                 [0, 1, 2])
@pytest.mark.parametrize('mass' , ['B_const_mass_M', 'B_M_brem_track_2'])
def test_prec_brem(mass : str, nbrem : int):
    '''
    Testing creation of PDF from MC sample with brem cut
    '''
    log.info('')

    cfg                      = _load_config('mc')
    out_dir                  = cfg['output']['out_dir']
    cfg['output']['out_dir'] = f'{out_dir}/test_prec_brem/{mass}_{nbrem:03}/v1'

    obs            = _get_obs(mass, cfg)
    cmp_prc        = cmp.get_prc(obs, nbrem, cfg)
    cmp_prc.run()
# --------------------------------------------------------------
@pytest.mark.parametrize('kind', ['exp', 'hypexp', 'modexp'])
def test_combinatorial(kind : str):
    '''
    Testing creation of PDF used for combinatorial
    '''
    log.info('')

    cfg            = _load_config(test='combinatorial')
    cfg['name']    = 'Combinatorial'
    out_dir        = cfg['out_dir']
    cfg['out_dir'] = f'{out_dir}/{kind}'

    obs= zfit.Space('B_M', limits=[4500, 6000])
    pdf= cmp.get_cb(obs=obs, kind=kind, cfg=cfg)
    print_pdf(pdf)
# --------------------------------------------------------------
@pytest.mark.parametrize('nbrem', [0, 1, 2, None])
@pytest.mark.parametrize('q2bin' , ['low', 'central', 'high'])
@pytest.mark.parametrize('sample', ['Bu_Kstee_Kpi0_eq_btosllball05_DPC', 'Bd_Kstee_eq_btosllball05_DPC', 'Bs_phiee_eq_Ball_DPC'])
def test_bxhsee(nbrem : int, q2bin : str, sample : str):
    '''
    Test B(u,d,s) -> K*ee
    '''
    log.info('')

    cfg                   = _load_config(test='bxhsee')
    cfg['input']['q2bin'] = q2bin

    obs = zfit.Space('B_M_brem_track_2', limits=(4500, 6000))
    pdf = cmp.get_kde(obs=obs, sample=sample, nbrem=nbrem, cfg=cfg)

    if pdf is None:
        return

    print_pdf(pdf)
# --------------------------------------------------------------
def test_bxhsee_first():
    '''
    Test Bu -> K*ee
    '''
    log.info('')

    nbrem = 1
    sample= 'Bu_Kstee_Kpi0_eq_btosllball05_DPC'

    cfg = _load_config(test='bxhsee_first')
    obs = zfit.Space('B_M_brem_track_2', limits=(4500, 6000))
    pdf = cmp.get_kde(obs=obs, sample=sample, nbrem=nbrem, cfg=cfg)

    if pdf is None:
        return

    print_pdf(pdf)
# --------------------------------------------------------------
