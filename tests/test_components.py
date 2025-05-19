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
from rx_selection           import selection  as sel
from rx_fitter              import components as cmp

log=LogStore.add_logger('rx_fitter:test_components')
# --------------------------------------------------------------
class Data:
    '''
    Data class
    '''
    mass    = 'B_M_brem_track_2'
    out_dir = '/tmp/tests/rx_fitter/components'
# --------------------------------------------------------------
@pytest.fixture(scope='session', autouse=True)
def _intialize():
    LogStore.set_level('rx_fitter:prec'              , 10)
    LogStore.set_level('rx_fitter:components'        , 10)
    LogStore.set_level('rx_fitter:mc_par_pdf'        , 10)
    LogStore.set_level('rx_calibration:fit_component', 10)
# --------------------------------------------------------------
def _load_config(test : str) -> dict:
    cfg_path = files('rx_fitter_data').joinpath(f'rare_fit/v1/rk_ee/{test}.yaml')
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
    log.info('')

    cfg                      = _load_config('signal')
    out_dir                  = cfg['output']['out_dir']
    cfg['output']['out_dir'] = f'{out_dir}/test_brem_definitions'

    d_cmp_set                = cfg['components']['Signal'][nbrem]
    d_cmp_set['create']      = False
    cfg['brem'][nbrem]       = _get_brem_cut(nbrem, kind)
    cfg['fitting']['range']  = _get_fitting_range(kind)
    cfg['input']['q2bin']    = 'central'

    obs = _get_obs(mass, cfg)
    pdf = cmp.get_mc(obs=obs, component_name='Signal', nbrem=nbrem, cfg=cfg)

    print_pdf(pdf)
# --------------------------------------------------------------
@pytest.mark.parametrize('nbrem', [0, 1, 2])
@pytest.mark.parametrize('mass' , ['B_M_brem_track_2'])
@pytest.mark.parametrize('name' , ['Signal'])
def test_mc_create(nbrem : int, mass : str, name : str):
    '''
    Testing creation of PDF from MC sample
    '''
    log.info('')

    cfg                      = _load_config('signal')
    out_dir                  = cfg['output']['out_dir']
    cfg['output']['out_dir'] = f'{out_dir}/test_mc_create'
    cfg['input']['q2bin']    = 'central'

    obs = _get_obs(mass, cfg)
    pdf = cmp.get_mc(obs=obs, component_name=name, nbrem=nbrem, cfg=cfg)

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

    cfg                      = _load_config('signal')
    out_dir                  = cfg['output']['out_dir']
    cfg['output']['out_dir'] = f'{out_dir}/test_mc_create'
    d_cmp_set                = cfg['components'][name][nbrem]
    d_cmp_set['create']      = False
    cfg['input']['q2bin']    = 'central'

    obs = _get_obs(mass, cfg)
    pdf = cmp.get_mc(obs=obs, component_name=name, nbrem=nbrem, cfg=cfg)

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

    cfg                      = _load_config('signal')
    out_dir                  = cfg['output']['out_dir']
    cfg['input']['q2bin']    = 'central'
    cfg['output']['out_dir'] = f'{out_dir}/test_mc_create'

    obs = _get_obs(mass, cfg)
    pdf = cmp.get_mc(obs=obs, component_name=name, nbrem=nbrem, cfg=cfg)

    print_pdf(pdf)
# --------------------------------------------------------------
@pytest.mark.parametrize('mass' , ['B_M_brem_track_2'])
@pytest.mark.parametrize('name' , ['Signal'])
def test_mc_reparametrized(mass : str, name : str):
    '''
    Testing creation of PDF from MC sample with tails fixed from other version
    '''
    log.info('')

    cfg                      = _load_config('signal')
    out_dir                  = cfg['output']['out_dir']
    cfg['input']['q2bin']    = 'central'
    cfg['output']['out_dir'] = f'{out_dir}/test_mc_create'

    obs  = _get_obs(mass, cfg)
    pdf  = cmp.get_mc_reparametrized(obs=obs, component_name=name, l_nbrem=[0, 1, 2], cfg=cfg)

    print_pdf(pdf)
# --------------------------------------------------------------
@pytest.mark.parametrize('mass'   , ['B_M_brem_track_2'])
@pytest.mark.parametrize('name'   , ['Signal'])
@pytest.mark.parametrize('kind'   , ['signal'])
@pytest.mark.parametrize('l_nbrem', [[0], [1], [2], [1,2], [0,1,2]])
def test_mc_brem_reparametrized(mass : str, name : str, kind : str, l_nbrem : list[str]):
    '''
    Test building full signal PDF for electron channel with Brem reparametrization
    '''
    log.info('')

    cfg                   = _load_config(kind)
    cfg['input']['q2bin'] = 'central'

    obs  = _get_obs(mass, cfg)
    pdf  = cmp.get_mc_reparametrized(obs=obs, component_name=name, cfg=cfg, l_nbrem=l_nbrem)

    brem = '_'.join(map(str, l_nbrem))
    print_pdf(pdf, txt_path=f'{Data.out_dir}/mc_brem_reparametrized/pdf_{brem}.txt')
# --------------------------------------------------------------
@pytest.mark.parametrize('nbrem', [[0], [1], [2], [0,1,2], [1,2]])
@pytest.mark.parametrize('mass' , ['B_const_mass_M', 'B_const_mass_psi2S_M', 'B_M_brem_track_2'])
@pytest.mark.parametrize('q2bin', ['jpsi', 'psi2'])
def test_prec_brem(mass : str, nbrem : list[int], q2bin : str):
    '''
    Testing creation of PDF from MC sample with brem cut
    '''
    log.info('')

    cfg                      = _load_config('prec')
    _set_brem_category(l_brem=nbrem, cfg=cfg)
    nbrem     = [ str(elm) for elm in nbrem ]
    brem_name = '_'.join(nbrem)

    cfg['input']['q2bin']    = q2bin
    cfg['output']['out_dir'] = f'{Data.out_dir}/prec/{mass}/{q2bin}/{brem_name}'

    obs     = _get_obs(mass, cfg)
    cmp_prc = cmp.get_prc(obs, cfg)

    if cmp_prc is not None:
        log.info(f'Component was built for mass/nbrem: {mass}/{nbrem}')
        cmp_prc.run()
    else:
        log.warning(f'No component was built for mass/nbrem: {mass}/{nbrem}')
# --------------------------------------------------------------
@pytest.mark.parametrize('q2bin', ['low', 'central', 'high'])
def test_combinatorial(q2bin : str):
    '''
    Testing creation of PDF used for combinatorial
    '''
    log.info('')

    cfg            = _load_config(test='combinatorial')
    cfg['out_dir'] = f'{Data.out_dir}/combinatorial/{q2bin}'

    obs= zfit.Space('B_M', limits=[4500, 6000])
    pdf= cmp.get_cb(obs=obs, q2bin=q2bin, cfg=cfg)
    print_pdf(pdf)
# --------------------------------------------------------------
@pytest.mark.parametrize('nbrem', [[0], [1], [2], [0,1,2], [1,2]])
@pytest.mark.parametrize('q2bin' , ['low', 'central', 'high'])
@pytest.mark.parametrize('sample', ['Bu_Kstee_Kpi0_eq_btosllball05_DPC', 'Bd_Kstee_eq_btosllball05_DPC', 'Bs_phiee_eq_Ball_DPC'])
def test_bxhsee(nbrem : list[int], q2bin : str, sample : str):
    '''
    Test B(u,d,s) -> K*ee
    '''
    log.info('')
    cfg                     = _load_config(test='bxhsee')
    _set_brem_category(l_brem=nbrem, cfg=cfg)
    nbrem     = [ str(elm) for elm in nbrem ]
    brem_name = '_'.join(nbrem)

    cfg['input']['q2bin']   = q2bin
    cfg['output']['out_dir']= f'{Data.out_dir}/bxhsee_{brem_name}'

    obs = zfit.Space('B_M_brem_track_2', limits=(4500, 6000))
    pdf = cmp.get_kde(obs=obs, sample=sample, cfg=cfg)

    if pdf is None:
        log.warning(f'No KDE can be built for {nbrem}/{q2bin}/{sample}')
        return

    print_pdf(pdf)
# --------------------------------------------------------------
@pytest.mark.parametrize('mass'  , ['B_M_smr_brem_track_2', 'B_M_brem_track_2'])
@pytest.mark.parametrize('nbrem' , [[0], [1], [2], [0,1,2], [1,2]])
@pytest.mark.parametrize('q2bin' , ['low', 'central', 'high'])
@pytest.mark.parametrize('sample', ['Bu_JpsiK_ee_eq_DPC', 'Bu_psi2SK_ee_eq_DPC'])
def test_cc_leakage(
        nbrem  : list[int],
        mass   : str,
        q2bin  : str,
        sample : str):
    '''
    Builds KDE for leaked ccbar component
    '''
    log.info('')
    cfg                      = _load_config(test='ccbar_leak')
    _set_brem_category(l_brem=nbrem, cfg=cfg)

    nbrem     = [ str(elm) for elm in nbrem ]
    brem_name = '_'.join(nbrem)

    cfg['input']['q2bin']    = q2bin
    cfg['output']['out_dir'] = f'{Data.out_dir}/leakage_{brem_name}'

    obs = zfit.Space(mass, limits=(4500, 6000))
    pdf = cmp.get_kde(obs=obs, sample=sample, cfg=cfg)

    if pdf is None:
        log.warning(f'No KDE can be built for {nbrem}/{q2bin}/{sample}')
        return

    print_pdf(pdf)
# --------------------------------------------------------------
def test_jpsi_leakage():
    '''
    Builds KDE for jpsi leakage in central bin
    '''
    q2bin = 'central'
    mass  = 'B_M_smr_brem_track_2'
    sample= 'Bu_JpsiK_ee_eq_DPC'

    log.info('')
    cfg = _load_config(test='ccbar_leak')

    cfg['input']['q2bin']    = q2bin
    cfg['output']['out_dir'] = f'{Data.out_dir}/leakage_jpsi'

    obs = zfit.Space(mass, limits=(4500, 6000))
    cmp.get_kde(obs=obs, sample=sample, cfg=cfg)
# --------------------------------------------------------------
def _set_brem_category(l_brem : list[int], cfg : dict) -> None:
    sel.reset_custom_selection()

    l_cut    = [ cfg['brem'][nbrem] for nbrem in l_brem ]
    l_cut    = [ f'({cut})'         for cut   in l_cut  ]
    brem_cut = ' || '.join(l_cut)

    log.info(f'Overriding selection with: {brem_cut}')

    sel.set_custom_selection(d_cut = {'brem' : brem_cut})
# --------------------------------------------------------------
