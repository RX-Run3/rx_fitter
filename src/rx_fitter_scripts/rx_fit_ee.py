'''
Script used to do fits to data
'''
import os
import copy
import argparse
from dataclasses         import dataclass
from importlib.resources import files

import ROOT
import zfit
import yaml
from ROOT                    import RDataFrame, RDF
from rx_data.rdf_getter      import RDFGetter
from rx_selection            import selection    as sel
from rx_calibration.hltcalibration.fit_component import FitComponent
from rx_calibration.hltcalibration.dt_fitter     import DTFitter

from dmu.logging.log_store   import LogStore
from dmu.rdataframe.atr_mgr  import AtrMgr
from dmu.stats.model_factory import ModelFactory

log = LogStore.add_logger('rx_fitter:rx_fit_ee')
# ---------------------------------
@dataclass
class Data:
    '''
    Class used to share attributes
    '''
    q2_bin  : str
    fit_cfg : dict
    min_mass = 4600
    max_mass = 6200

    l_col    = [
            'hop_mass',
            'B_M',
            'Jpsi_M',
            'mva_cmb',
            'mva_prc',
            'swp_cascade_mass_swp',
            'swp_jpsi_misid_mass_swp']

    trigger  = 'Hlt2RD_BuToKpEE_MVA'
    cmb_wp   = 0.00
    prc_wp   = 0.00
    obs      = zfit.Space('mass', limits=(min_mass, max_mass))

    RDFGetter.samples= {
            'main' : '/home/acampove/external_ssd/Data/samples/main.yaml',
            'mva'  : '/home/acampove/external_ssd/Data/samples/mva.yaml',
            'jps'  : '/home/acampove/external_ssd/Data/samples/jpsi_misid.yaml',
            'cas'  : '/home/acampove/external_ssd/Data/samples/cascade.yaml',
            'hop'  : '/home/acampove/external_ssd/Data/samples/hop.yaml'} 

    cache_dir = '/home/acampove/external_ssd/Data/cache/rx_fit/ee'

    l_use_cut = [
            'jpsi_misid',
            'cascade',
            'bdt',
            'hop',
            'q2',
            'mass']

    mc_cfg = {
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

    dt_cfg = {
            'error_method' : 'minuit_hesse',
            'out_dir'      : 'plots/fit/data',
            'plotting'     :
            {
                'nbins'   : 50,
                'stacked' : True,
                'd_leg'   : {
                    'SumPDF_ext'     : r'$B^+\to K^+\mu^+\mu^-$',
                    'Exponential_ext': 'Combinatorial',
                    }
                },
            }
# ---------------------------------
def _parse_args() -> None:
    parser = argparse.ArgumentParser(description='Script used to fit mass distributions')
    parser.add_argument('-q', '--q2bin' , type=str, help='q2 bin' , choices=['low', 'central', 'jpsi', 'psi2', 'high'], required=True)
    args = parser.parse_args()

    Data.q2_bin = args.q2bin
# ---------------------------------
def _load_rdf(sample : str) -> RDataFrame:
    sample_name = sample.replace('*', 'p')
    out_path    = f'{Data.cache_dir}/{sample_name}_{Data.q2_bin}.root'
    if os.path.isfile(out_path):
        log.info('DataFrame already cached, reloading')
        rdf = RDataFrame('tree', out_path)
        return rdf

    log.info('DataFrame not cached')

    gtr = RDFGetter(sample=sample, trigger=Data.trigger)
    rdf = gtr.get_rdf(columns=Data.l_col)
    _   = AtrMgr(rdf)

    rdf    = _apply_selection(rdf, process = sample)
    d_data = rdf.AsNumpy(['B_M', 'mva_cmb', 'mva_prc'])

    rdf=RDF.FromNumpy(d_data)
    rdf.Snapshot('tree', out_path)

    return rdf
# ---------------------------------
def _get_rdf(sample : str) -> RDataFrame:
    rdf = _load_rdf(sample)
    rdf = rdf.Filter(f'mva_cmb > {Data.cmb_wp}', 'CMB')
    rdf = rdf.Filter(f'mva_prc > {Data.prc_wp}', 'PRC')
    rdf = rdf.Define('mass', 'B_M')

    rep = rdf.Report()
    rep.Print()

    return rdf
# ---------------------------------
def _apply_selection(rdf : RDataFrame, process : str) -> RDataFrame:
    d_sel = sel.selection(project='RK', analysis='EE', q2bin=Data.q2_bin, process=process)
    d_use = { name : cut for name, cut in d_sel.items() if name in Data.l_use_cut }

    for name, cut in d_use.items():
        if name == 'bdt':
            continue
        rdf = rdf.Filter(cut, name)

    log.info(f'Using dataframe for {Data.q2_bin} q2 bin')
    rep = rdf.Report()
    rep.Print()

    return rdf
# ---------------------------------
def _initialize() -> None:
    os.makedirs(Data.cache_dir, exist_ok=True)
    cfg_path = files('rx_fitter_data').joinpath('config/v1.yaml')
    with open(cfg_path, encoding='utf-8') as ifile:
        Data.fit_cfg = yaml.safe_load(ifile)
# ---------------------------------
def _get_mc(sample : str) -> FitComponent:
    cfg            = copy.deepcopy(Data.mc_cfg)
    cfg['name']    = sample
    out_dir        = cfg['out_dir']
    cfg['out_dir'] = f'{out_dir}/{Data.q2_bin}/{sample}'

    rdf   = _get_rdf(sample)
    rdf   = rdf.Define('weights', '1')

    l_pdf, l_shr = _get_fitting_model(sample)
    if l_pdf == ['kde']:
        pdf = None
    else:
        mod   = ModelFactory(obs = Data.obs, l_pdf = l_pdf, l_shared=l_shr)
        pdf   = mod.get_pdf()

    obj   = FitComponent(cfg=cfg, rdf=rdf, pdf=pdf, obs=Data.obs)

    return obj
# ---------------------------------
def _get_fitting_model(sample : str) -> tuple[list[str],list[str]]:
    cfg = Data.fit_cfg[Data.q2_bin][sample]

    return cfg['model'], cfg['shared']
# ---------------------------------
def _get_combinatorial() -> FitComponent:
    cfg            = copy.deepcopy(Data.mc_cfg)
    out_dir        = cfg['out_dir']
    cfg['out_dir'] = f'{out_dir}/{Data.q2_bin}'

    del cfg['fitting']
    cfg['name'] = 'comb'

    mod   = ModelFactory(obs = Data.obs, l_pdf = ['exp'], l_shared=[])
    pdf   = mod.get_pdf()
    obj   = FitComponent(cfg=cfg, rdf=None, pdf=pdf)

    return obj
# ---------------------------------
def _get_cfg() -> dict:
    cfg            = copy.deepcopy(Data.dt_cfg)
    out_dir        = cfg['out_dir']
    cfg['out_dir'] = f'{out_dir}/{Data.q2_bin}'

    return cfg
# ---------------------------------
def main():
    '''
    Process starts here
    '''
    _parse_args()
    _initialize()

    cmp_sig = _get_signal()
    cmp_cmb = _get_combinatorial()
    rdf     = _get_rdf(is_mc=False)
    cfg     = _get_cfg()

    obj = DTFitter(rdf = rdf, components = [cmp_cmb, cmp_sig], cfg=cfg)
    obj.fit()
# ---------------------------------
if __name__ == '__main__':
    main()
