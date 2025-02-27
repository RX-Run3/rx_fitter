'''
Module with functions needed to provide fit components
'''
# pylint: disable=too-many-positional-arguments, too-many-function-args, too-many-arguments

import os
import copy

from zfit.core.interfaces                        import ZfitSpace as zobs
from ROOT                                        import RDataFrame
from dmu.stats.model_factory                     import ModelFactory
from dmu.logging.log_store                       import LogStore
from rx_selection                                import selection as sel
from rx_data.rdf_getter                          import RDFGetter
from rx_calibration.hltcalibration.fit_component import FitComponent
from rx_fitter.prec                              import PRec

log = LogStore.add_logger('rx_fitter:components')
# ------------------------------------
class Data:
    '''
    Data class
    '''
    cache_dir = '/tmp/cache/rx_fits'
    cfg       = {
            'out_dir': 'plots/fit',
            'fitting':
            {
                'error_method'  : 'minuit_hesse',
                'weights_column': 'weights',
                'ntries'        : 20,
                'pvalue'        : 0.02,
                },
            'plotting' :
            {
                'nbins'   : 50,
                'stacked' : True,
                },
            }
# ---------------------------------
def _update_selection(d_sel : dict[str,str]) -> dict[str,str]:
    if 'selection' not in Data.cfg:
        log.info('Not updating selection')
        return d_sel

    log.info('Updating selection')
    d_cut = Data.cfg['selection']
    d_sel.update(d_cut)

    return d_sel
# ---------------------------------
def _get_rdf(sample : str, q2bin : str, trigger : str, nbrem : int) -> RDataFrame:
    out_path = f'{Data.cache_dir}/{sample}_{q2bin}.root'
    if os.path.isfile(out_path):
        log.info('DataFrame already cached, reloading')
        rdf = RDataFrame('tree', out_path)
        return rdf

    log.info('DataFrame not cached')

    gtr = RDFGetter(sample=sample, trigger=trigger)
    rdf = gtr.get_rdf()

    analysis = 'MM' if 'MuMu' in trigger else 'EE'

    d_sel = sel.selection(project='RK', analysis=analysis, q2bin=q2bin, process=sample)
    d_sel = _update_selection(d_sel)
    for cut_name, cut_value in d_sel.items():
        log.info(f'{cut_name:<20}{cut_value}')
        if cut_name == 'mass':
            cut_value = '(1)'

        rdf = rdf.Filter(cut_value, cut_name)

    rdf = rdf.Define('nbrem', 'L1_BremMultiplicity + L2_BremMultiplicity')

    if   nbrem == -1:
        pass
    elif nbrem in [0, 1]:
        rdf = rdf.Filter(f'nbrem == {nbrem}', 'nbrem')
    elif nbrem == 2:
        rdf = rdf.Filter(f'nbrem >= {nbrem}', 'nbrem')
    else:
        raise ValueError(f'Invalid brem argument: {nbrem}')

    rep = rdf.Report()
    rep.Print()

    return rdf
# ------------------------------------
def _get_model(sample : str, q2bin : str, trigger : str, nbrem : int, model : list[str]) -> list[str]:
    if model is not None:
        return model

    log.info('Model not passed, will pick default')

    is_sig = sample  == 'Bu_JpsiK_ee_eq_DPC'
    is_trg = trigger == 'Hlt2RD_BuToKpEE_MVA'
    is_jps = q2bin   == 'jpsi'
    is_brm = nbrem   in [0, 1, 2]

    if is_sig and is_jps and is_brm and is_trg:
        return {
                0 : ['suj', 'suj'],
                1 : ['suj', 'cbr'],
                2 : ['suj', 'suj']}[nbrem]


    if sample == 'Bu_JpsiPi_ee_eq_DPC':
        return {
                0 : ['suj'],
                1 : ['suj'],
                2 : ['suj']}[nbrem]

    raise ValueError(f'Cannot assign default model for: {sample}/{q2bin}/{trigger}/{nbrem}')
# ------------------------------------
def get_mc(obs : zobs, sample : str, q2bin : str, trigger : str, nbrem : int, model : list[str] = None) -> FitComponent:
    '''
    Will return FitComponent object for given MC sample
    '''
    model          = _get_model(sample, q2bin, trigger, nbrem, model)
    brem_name      = 'all' if nbrem not in [0, 1, 2] else nbrem
    model_name     = '_'.join(model)
    mass           = obs.obs[0]
    cfg            = copy.deepcopy(Data.cfg)
    cfg['name']    = sample
    out_dir        = cfg['out_dir']
    cfg['out_dir'] = f'{out_dir}/{q2bin}/{sample}_{trigger}/{mass}_{brem_name}/{model_name}'

    rdf   = _get_rdf(sample, q2bin, trigger, nbrem)
    rdf   = rdf.Define('weights', '1')

    mod   = ModelFactory(sample, obs, model, ['mu', 'sg'])
    pdf   = mod.get_pdf()

    obj   = FitComponent(cfg=cfg, rdf=rdf, pdf=pdf, obs=obs)

    return obj
# ------------------------------------
def get_prc(obs : zobs, q2bin : str, trigger : str, nbrem : int) -> FitComponent:
    '''
    Function returning FitComponent object for Partially reconstructed background
    '''
    mass        = obs.obs[0]
    cfg         = copy.deepcopy(Data.cfg)
    cfg['name'] = 'PRec'
    out_dir     = cfg['out_dir']
    cfg['out_dir'] = f'{out_dir}/{trigger}_{q2bin}_{nbrem:03}'

    bw     = {'jpsi' :  5, 'psi2' : 10}[q2bin]
    l_samp = [
            'Bu_JpsiX_ee_eq_JpsiInAcc',
            'Bd_JpsiX_ee_eq_JpsiInAcc',
            'Bs_JpsiX_ee_eq_JpsiInAcc']

    d_wgt= {'dec' : 1, 'sam' : 1}
    obj  = PRec(samples=l_samp, trig=trigger, q2bin=q2bin, d_weight=d_wgt)

    if   nbrem == -1:
        pass
    elif nbrem in [0, 1]:
        obj.cuts = {'nbrem' : f'nbrem == {nbrem}'}
    elif nbrem == 2:
        obj.cuts = {'nbrem' : f'nbrem >= {nbrem}'}
    else:
        raise ValueError(f'Invalid brem argument: {nbrem}')

    pdf=obj.get_sum(mass=mass, name='PRec', obs=obs, bandwidth=bw)
    fcm= FitComponent(cfg=cfg, rdf=None, pdf=pdf)

    return fcm
# ------------------------------------
def get_cb(obs : zobs, kind : str) -> FitComponent:
    '''
    Returns fit component for combinatorial fit
    '''
    cfg            = copy.deepcopy(Data.cfg)
    cfg['name']    = 'combinatorial'
    out_dir        = cfg['out_dir']
    cfg['out_dir'] = f'{out_dir}/{kind}'

    mod   = ModelFactory(preffix='', obs=obs, l_pdf = [kind], l_shared = [])
    pdf   = mod.get_pdf()

    obj   = FitComponent(cfg=cfg, rdf=None, pdf=pdf, obs=obs)

    return obj
# ------------------------------------
