'''
Module with functions needed to provide fit components
'''
# pylint: disable=too-many-positional-arguments, too-many-function-args, too-many-arguments, too-many-locals

import os
import copy

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
def get_rdf(sample : str, q2bin : str, trigger : str, cuts : dict[str,str] = None) -> RDataFrame:
    '''
    Function that returns a ROOT dataframe for a given dataset, MC or real data
    '''
    gtr = RDFGetter(sample=sample, trigger=trigger)
    rdf = gtr.get_rdf()

    analysis = 'MM' if 'MuMu' in trigger else 'EE'

    d_sel=sel.selection(project='RK', analysis=analysis, q2bin=q2bin, process=sample)
    d_sel.update(cuts)
    for cut_name, cut_value in d_sel.items():
        log.debug(f'{cut_name:<20}{cut_value}')
        rdf = rdf.Filter(cut_value, cut_name)

    if log.getEffectiveLevel() < 20:
        rep = rdf.Report()
        rep.Print()

    return rdf
# ------------------------------------
def _get_cuts(nbrem : int, cfg : dict) -> dict[str,str]:
    d_cut          = {}
    d_cut['nbrem'] = cfg['brem'][nbrem]

    if 'input'     not in cfg:
        return d_cut

    if 'selection' not in cfg['input']:
        return d_cut

    log.warning('Overriding default selection')
    cuts= cfg['input']['selection']
    for name, expr in cuts.items():
        log.debug(f'{name:<20}{expr}')

    cuts.update(d_cut)

    return cuts
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

    d_cmp   = cfg['fitting']['config'][component_name]
    sample  = d_cmp['sample']
    cuts    = _get_cuts(nbrem, cfg)
    rdf     = get_rdf(sample, q2bin, trigger, cuts)
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

    d_inp                 = cfg['input']
    d_cmp                 = cfg['fitting']['config'][component_name]
    d_fit                 = d_cmp['fitting']
    d_plt                 = d_cmp['plotting']
    trigger               = d_inp['trigger']
    q2bin                 = d_inp['q2bin'  ]

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
def _get_brem_reparametrization(pdf_z : zpdf, pdf_o : zpdf, pdf_t : zpdf) -> zpdf:
    frac_z = zfit.Parameter('frac_brem_000', 0.3, 0, 1)
    frac_o = zfit.Parameter('frac_brem_001', 0.4, 0, 1)
    frac_t = zfit.Parameter('frac_brem_002', 0.3, 0, 1)

    pdf = zfit.pdf.SumPDF(pdfs=[pdf_z, pdf_o, pdf_t], fracs=[frac_z, frac_o, frac_t])

    return pdf
# ------------------------------------
def get_mc_reparametrized(obs : zobs, component_name : str, cfg : dict, nbrem : int) -> zpdf:
    '''
    Will return reparametrized fit component. The MC fit is expected to have been done already and this
    will only load those parameters:

    - No RDF needed
    - No plotting needed
    '''
    if nbrem is None:
        pdf_z = get_mc_reparametrized(obs, component_name, cfg, 0)
        pdf_o = get_mc_reparametrized(obs, component_name, cfg, 1)
        pdf_t = get_mc_reparametrized(obs, component_name, cfg, 2)
        pdf   = _get_brem_reparametrization(pdf_z=pdf_z, pdf_o=pdf_o, pdf_t=pdf_t)

        return pdf

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
def get_prc(obs : zobs, nbrem : int, cfg : dict) -> FitComponent:
    '''
    Function returning FitComponent object for Partially reconstructed background
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
    obj.cuts = _get_cuts(nbrem, cfg)

    pdf=obj.get_sum(mass=mass, name='PRec', obs=obs, **cfg_kde)

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
def _get_brem_cut(cfg : dict, l_nbrem : list[int]) -> str:
    l_cut = [ cfg['brem'][nbrem] for nbrem in l_nbrem ]
    cut   = '||'.join(l_cut)

    log.info(f'Using brem cut: {cut}')

    return cut
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

    d_cut = {}
    if 'selection' in cfg:
        d_sel = cfg['selection']
        d_cut.update(d_sel)

    d_cut['nbrem'] = _get_brem_cut(cfg=cfg, l_nbrem=l_nbrem)

    rdf = get_rdf(sample=sample, q2bin=q2bin, trigger=trigger, cuts=d_cut)
    fcm = FitComponent(cfg=cfg, rdf=rdf, pdf=None, obs=obs)
    pdf = fcm.get_pdf()

    return pdf
# ------------------------------------
