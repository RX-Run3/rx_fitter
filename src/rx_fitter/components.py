'''
Module with functions needed to provide fit components
'''
# pylint: disable=too-many-positional-arguments, too-many-function-args

import os
import copy

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
    cache_dir = '/home/acampove/Data/RX_run3/cache/rx_fits'
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
def _get_rdf(sample : str, q2bin : str, trigger : str) -> RDataFrame:
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
    for cut_name, cut_value in d_sel.items():
        log.info(f'{cut_name:<20}{cut_value}')
        rdf = rdf.Filter(cut_value, cut_name)

    return rdf
# ------------------------------------
def get_mc(obs, sample : str, q2bin : str, trigger : str, model : list[str]) -> FitComponent:
    '''
    Will return FitComponent object for given MC sample
    '''
    mass           = obs.obs[0]
    cfg            = copy.deepcopy(Data.cfg)
    cfg['name']    = sample
    out_dir        = cfg['out_dir']
    cfg['out_dir'] = f'{out_dir}/{q2bin}/{sample}_{trigger}/{mass}'

    rdf   = _get_rdf(sample, q2bin, trigger)
    rdf   = rdf.Define('weights', '1')

    mod   = ModelFactory(sample, obs, model, ['mu', 'sg'])
    pdf   = mod.get_pdf()

    obj   = FitComponent(cfg=cfg, rdf=rdf, pdf=pdf, obs=obs)

    return obj
# ------------------------------------
def get_prc(obs, q2bin : str, trigger : str) -> FitComponent:
    '''
    Function returning FitComponent object for Partially reconstructed background
    '''
    mass   = obs.obs[0]
    cfg    = copy.deepcopy(Data.cfg)
    bw     = {'jpsi' :  5, 'psi2' : 10}[q2bin]
    l_samp = [
            'Bu_JpsiX_ee_eq_JpsiInAcc',
            'Bd_JpsiX_ee_eq_JpsiInAcc',
            'Bs_JpsiX_ee_eq_JpsiInAcc']

    d_wgt= {'dec' : 1, 'sam' : 1}
    obj=PRec(samples=l_samp, trig=trigger, q2bin=q2bin, d_weight=d_wgt)
    pdf=obj.get_sum(mass=mass, name='PRec', obs=obs, bandwidth=bw)
    fcm= FitComponent(cfg=cfg, rdf=None, pdf=pdf)

    return fcm
# ------------------------------------
