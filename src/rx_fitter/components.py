'''
Module with functions needed to provide fit components
'''
# pylint: disable=too-many-positional-arguments, too-many-function-args, too-many-arguments, too-many-locals

import os
import copy

from zfit.core.interfaces                        import ZfitSpace as zobs
from ROOT                                        import RDataFrame
from dmu.generic                                 import version_management as vman
from dmu.stats.model_factory                     import ModelFactory
from dmu.logging.log_store                       import LogStore
from rx_selection                                import selection as sel
from rx_data.rdf_getter                          import RDFGetter
from rx_calibration.hltcalibration.fit_component import FitComponent
from rx_fitter.mc_par_pdf                        import MCParPdf
from rx_fitter.prec                              import PRec

log = LogStore.add_logger('rx_fitter:components')
# ------------------------------------
class Data:
    '''
    Data class
    '''
    cache_dir = '/tmp/cache/rx_fits'
# ---------------------------------
def get_rdf(sample : str, q2bin : str, trigger : str, cuts : dict[str,str] = None) -> RDataFrame:
    '''
    Function that returns a ROOT dataframe for a given dataset, MC or real data
    '''
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

    if cuts is not None:
        log.warning('Overriding default selection')
        for name, expr in cuts.items():
            log.info(f'   {name:<20}{expr}')
            rdf = rdf.Filter(expr, name)

    rep = rdf.Report()
    rep.Print()

    return rdf
# ------------------------------------
def _get_last_version(path : str) -> str:
    [init, fnal] = path.split('VERS')
    init         = vman.get_last_version(dir_path=init, version_only=False)

    return f'{init}{fnal}'
# ------------------------------------
def get_mc(obs : zobs, **kwargs) -> FitComponent:
    '''
    Will return FitComponent object for given MC sample
    '''
    cfg     = copy.deepcopy(Data.cfg)
    cfg.update(kwargs)

    name    = kwargs['name'   ]
    nbrem   = kwargs['nbrem'  ]
    q2bin   = kwargs['q2bin'  ]
    trigger = kwargs['trigger']

    bcut  = f'nbrem == {nbrem}' if nbrem in [0, 1] else f'nbrem >= {nbrem}'
    d_cut = {'nbrem' : bcut}
    rdf   = get_rdf(name, q2bin, trigger, d_cut)
    rdf   = rdf.Define('weights', '1')

    obj   = MCParPdf(rdf=rdf, obs=obs, cfg=cfg)

    return obj.get_fcomp()
# ------------------------------------
def get_prc(name : str, obs : zobs, q2bin : str, trigger : str, cuts : dict[str,str] = None, bw : int = None) -> FitComponent:
    '''
    Function returning FitComponent object for Partially reconstructed background
    '''
    mass        = obs.obs[0]
    cfg         = copy.deepcopy(Data.cfg)
    cfg['name'] = 'PRec'
    out_dir        = f'{Data.fit_dir}/mc/{q2bin}/VERS/binclusive_{trigger}/{mass}_{name}/kde'
    cfg['out_dir'] = _get_last_version(out_dir)

    bw     = {'jpsi' :  5, 'psi2' : 10}[q2bin] if bw is None else bw

    l_samp = [
            'Bu_JpsiX_ee_eq_JpsiInAcc',
            'Bd_JpsiX_ee_eq_JpsiInAcc',
            'Bs_JpsiX_ee_eq_JpsiInAcc']

    d_wgt= {'dec' : 1, 'sam' : 1}
    obj  = PRec(samples=l_samp, trig=trigger, q2bin=q2bin, d_weight=d_wgt)
    if cuts is not None:
        obj.cuts = cuts

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
    cfg['out_dir'] = f'/tmp/components/{kind}'

    mod   = ModelFactory(preffix='cmb', obs=obs, l_pdf = [kind], l_shared = [], l_float= [])
    pdf   = mod.get_pdf()

    obj   = FitComponent(cfg=cfg, rdf=None, pdf=pdf, obs=obs)
    obj.run()

    return obj
# ------------------------------------
