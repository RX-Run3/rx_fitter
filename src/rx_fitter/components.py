'''
Module with functions needed to provide fit components
'''
# pylint: disable=too-many-positional-arguments, too-many-function-args, too-many-arguments, too-many-locals

import os
import copy
from typing import Union

import zfit
import pandas as pnd
from zfit.core.interfaces                        import ZfitSpace as zobs
from zfit.core.basepdf                           import BasePDF   as zpdf
from ROOT                                        import RDataFrame, RDF
from dmu.stats.model_factory                     import ModelFactory
from dmu.logging.log_store                       import LogStore
from dmu.generic                                 import hashing
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
def get_rdf(sample : str, q2bin : str, trigger : str) -> RDataFrame:
    '''
    Function that returns a ROOT dataframe for a given dataset, MC or real data
    '''
    gtr  = RDFGetter(sample=sample, trigger=trigger)
    rdf  = gtr.get_rdf()
    d_sel=sel.selection(trigger=trigger, q2bin=q2bin, process=sample)

    for cut_name, cut_value in d_sel.items():
        log.debug(f'{cut_name:<20}{cut_value}')
        rdf = rdf.Filter(cut_value, cut_name)

    if log.getEffectiveLevel() < 20:
        rep = rdf.Report()
        rep.Print()

    return rdf
# ------------------------------------
def _get_mc_rdf(cfg : dict, component_name : str, nbrem : int) -> RDataFrame:
    d_cmp_set = cfg['components'][component_name][nbrem]
    if not d_cmp_set['create']:
        log.info('Will not redo fit, not recalculating dataframe')
        return None

    log.info('Making ROOT dataframe with input data to run fit')
    d_inp   = cfg['input']
    trigger = d_inp['trigger']
    q2bin   = d_inp['q2bin'  ]
    brem_cut= cfg['brem'][nbrem]

    d_cmp   = cfg['fitting']['config'][component_name]
    sample  = d_cmp['sample']
    rdf     = get_rdf(sample, q2bin, trigger)
    # This is needed to get fit per-brem category
    # It does not override analysis selection
    rdf     = rdf.Filter(brem_cut, f'Brem {nbrem}')
    rdf     = rdf.Define('weights', '1')

    if 'max_entries' in cfg['input']:
        max_entries = cfg['input']['max_entries']
        log.warning(f'Limitting dataframe to {max_entries} entries')
        rdf = rdf.Range(max_entries)

    return rdf
# ------------------------------------
def get_mc(obs : zobs, component_name : str, nbrem : int, cfg : dict) -> FitComponent:
    '''
    Will return FitComponent object for given MC sample
    '''
    cfg     = copy.deepcopy(cfg)
    rdf     = _get_mc_rdf(cfg, component_name, nbrem)

    d_inp   = cfg['input']
    d_cmp   = cfg['fitting']['config'][component_name]
    d_fit   = d_cmp['fitting']
    d_plt   = d_cmp['plotting']
    trigger = d_inp['trigger']
    q2bin   = d_inp['q2bin'  ]

    cfg['component_name'] = f'{component_name}_{nbrem:03}'
    cfg['q2bin'  ]        = q2bin
    cfg['trigger']        = trigger
    cfg['nbrem'  ]        = nbrem

    cmp_cfg         = cfg['components'][component_name][nbrem]
    if 'fvers' in cmp_cfg:
        cfg['fvers'] = cmp_cfg['fvers']

    cfg['create'  ] = cmp_cfg['create' ]
    cfg['shared'  ] = cmp_cfg['shared' ]
    cfg['model'   ] = cmp_cfg['model'  ]
    cfg['pfloat'  ] = cmp_cfg['pfloat' ]
    cfg['fitting' ] = d_fit
    cfg['plotting'] = d_plt
    cfg['fitting']['weights_column'] = cmp_cfg['weights']

    obj   = MCParPdf(rdf=rdf, obs=obs, cfg=cfg)

    return obj.get_pdf()
# ------------------------------------
def _get_mc_reparametrized_brem(obs : zobs, component_name : str, cfg : dict, nbrem : int) -> zpdf:
    cfg     = copy.deepcopy(cfg)
    d_inp   = cfg['input']
    trigger = d_inp['trigger']
    q2bin   = d_inp['q2bin'  ]

    d_cmp   = cfg['fitting']['config'][component_name]
    d_fit   = d_cmp['fitting']

    cfg['component_name'] = component_name
    cfg['q2bin'  ]        = q2bin
    cfg['trigger']        = trigger
    cfg['nbrem'  ]        = nbrem
    cmp_cfg               = cfg['components'][component_name][nbrem]

    if 'fvers' in cmp_cfg:
        cfg['fvers'] = cmp_cfg['fvers']

    cfg['reparametrize'] = cmp_cfg['reparametrize']

    cfg['create'  ] = cmp_cfg['create' ]
    cfg['shared'  ] = cmp_cfg['shared' ]
    cfg['model'   ] = cmp_cfg['model'  ]
    cfg['pfloat'  ] = cmp_cfg['pfloat' ]
    cfg['fitting' ] = d_fit
    cfg['fitting']['weights_column'] = cmp_cfg['weights']

    obj = MCParPdf(rdf=None, obs=obs, cfg=cfg)
    pdf = obj.get_pdf()

    return pdf
# ------------------------------------
def get_mc_reparametrized(obs : zobs, component_name : str, cfg : dict, l_nbrem : list[int]) -> zpdf:
    '''
    Will return reparametrized fit component. The MC fit is expected to have been done already and this
    will only load those parameters:

    - No RDF needed
    - No plotting needed
    '''

    l_pdf = [ _get_mc_reparametrized_brem(obs, component_name, cfg, nbrem) for nbrem in l_nbrem ]
    if len(l_pdf) == 1:
        return l_pdf[0]

    l_frc = [       zfit.Parameter(f'frac_brem_{nbrem:03}', 0.3, 0, 1)     for nbrem in l_nbrem ]
    pdf   = zfit.pdf.SumPDF(pdfs=l_pdf, fracs=l_frc)

    return pdf
# ------------------------------------
def get_prc(obs : zobs, nbrem : int, cfg : dict) -> Union[FitComponent,None]:
    '''
    Function returning FitComponent object for Partially reconstructed background
    build from cocktail charmonium MC, NOT the exclusive rare MC
    '''
    mass     = obs.obs[0]
    q2bin    = cfg['input']['q2bin']
    trigger  = cfg['input']['trigger']
    l_samp   = cfg['fitting']['config']['PRec']['sample'  ]
    d_wgt    = cfg['fitting']['config']['PRec']['weights' ]
    d_plt    = cfg['fitting']['config']['PRec']['plotting']
    cfg_kde  = cfg['fitting']['config']['PRec']['cfg_kde' ]
    out_dir  = cfg['output']['out_dir']

    obj      = PRec(samples=l_samp, trig=trigger, q2bin=q2bin, d_weight=d_wgt)
    obj.cuts = { 'nbrem' : cfg['brem'][nbrem] }
    pdf      = obj.get_sum(mass=mass, name='PRec', obs=obs, **cfg_kde)

    if pdf is None:
        log.warning('No PRec PDF can be built')
        return None

    cfg['name']    = 'PRec'
    cfg['plotting']= d_plt
    cfg['out_dir'] = out_dir

    fcm= FitComponent(cfg=cfg, rdf=None, pdf=pdf)

    return fcm
# ------------------------------------
def get_cb(obs : zobs, q2bin : str, cfg : dict) -> FitComponent:
    '''
    Returns fit component for combinatorial fit
    '''
    kind = cfg['q2'][q2bin]['model']
    cfg['name'] = 'Combinatorial'

    d_fix= None
    if 'fix' in cfg['q2'][q2bin]:
        d_fix= cfg['q2'][q2bin]['fix']

    mod  = ModelFactory(preffix='cmb', obs=obs, l_pdf = [kind], l_shared = [], l_float= [], d_fix=d_fix)
    pdf  = mod.get_pdf()

    obj  = FitComponent(cfg=cfg, rdf=None, pdf=pdf, obs=obs)

    return obj.get_pdf()
# ------------------------------------
def get_kde(obs : zobs, sample : str, l_nbrem : list[int], cfg : dict) -> zpdf:
    '''
    Function returning zfit PDF object for Samples that need to be modelled with a KDE

    obs    : zfit observable
    sample : Sample name, e.g.
    l_nbrem: Brem category list e.g. [0, 1, 2]
    cfg    : Dictionary with configuration
    '''

    hsh = hashing.hash_object(obj=[obs.to_json(), sample, l_nbrem, cfg])

    nbrem    = '_'.join(map(str, l_nbrem))
    mass     = obs.obs[0]
    q2bin    = cfg['input']['q2bin']
    trigger  = cfg['input']['trigger']
    d_plt    = cfg['fitting']['config'][sample]['plotting']
    out_dir  = cfg['output']['out_dir']
    out_dir  = f'{out_dir}/{sample}/{q2bin}/{mass}_{nbrem}/{hsh}'

    d_plt['title'] = f'{sample}; {l_nbrem}'
    cfg['name']    = sample
    cfg['plotting']= d_plt
    cfg['out_dir'] = out_dir

    if os.path.isfile(f'{out_dir}/data.json'):
        data_path = f'{out_dir}/data.json'
        log.debug(f'JSON file with data found, loading: {data_path}')

        df        = pnd.read_json(data_path)
        rdf       = RDF.FromPandas(df)
        fcm       = FitComponent(cfg=cfg, rdf=rdf, pdf=None, obs=obs)

        return fcm.get_pdf()

    rdf = get_rdf(sample=sample, q2bin=q2bin, trigger=trigger)
    fcm = FitComponent(cfg=cfg, rdf=rdf, pdf=None, obs=obs)
    pdf = fcm.get_pdf()

    return pdf
# ------------------------------------
