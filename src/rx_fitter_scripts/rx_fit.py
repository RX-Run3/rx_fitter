'''
Script used to do fits to data
'''
import os
import copy
import argparse
from dataclasses         import dataclass

import ROOT
import zfit
from ROOT                    import RDataFrame, RDF
from rx_data.rdf_getter      import RDFGetter
from rx_selection.selection  import load_selection_config
from rx_calibration.hltcalibration.fit_component import FitComponent
from rx_calibration.hltcalibration.dt_fitter     import DTFitter

from dmu.logging.log_store   import LogStore
from dmu.rdataframe.atr_mgr  import AtrMgr
from dmu.stats.model_factory import ModelFactory
from dmu.stats.utilities     import print_pdf

log = LogStore.add_logger('rx_fitter:rx_fit')
# ---------------------------------
@dataclass
class Data:
    '''
    Class used to share attributes
    '''
    q2_bin  : str

    trigger  = 'Hlt2RD_BuToKpMuMu_MVA'
    mss_mm   = '(B_M > 5180) && (B_M < 5600)'
    mva_cmb  = 'mva.mva_cmb > 0.8'
    mva_prc  = 'mva.mva_prc > 0.6'
    obs      = zfit.Space('mass', limits=(5180, 5600))

    RDFGetter.samples_dir = '/home/acampove/Data/RX_run3/NO_q2_bdt_mass_Q2_central_VR_v1'
    cache_dir             = '/home/acampove/Data/RX_run3/cache/rx_fits'

    mc_cfg = {
            'name'   : 'signal',
            'out_dir': 'plots/fit/signal',
            'fitting':
            {
                'error_method'  : 'minuit_hesse',
                'weights_column': 'weights',
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
def _get_rdf(is_mc : bool) -> RDataFrame:
    kind     = 'mc' if is_mc else 'data'
    out_path = f'{Data.cache_dir}/{kind}_{Data.q2_bin}.root'
    if os.path.isfile(out_path):
        log.info('DataFrame already cached, reloading')
        rdf = RDataFrame('tree', out_path)
        return rdf

    log.info('DataFrame not cached')
    if is_mc:
        sample = 'Bu_Kmumu_eq_btosllball05_DPC'
    else:
        sample = 'DATA_24_Mag*_24c*'

    gtr = RDFGetter(sample=sample, trigger=Data.trigger)
    rdf = gtr.get_rdf()
    _   = AtrMgr(rdf)

    rdf      = _apply_selection(rdf)
    arr_mass = rdf.AsNumpy(['B_M'])['B_M']

    rdf=RDF.FromNumpy({'mass' : arr_mass})
    rdf.Snapshot('tree', out_path)

    return rdf
# ---------------------------------
def _apply_selection(rdf : RDataFrame) -> RDataFrame:
    cfg = load_selection_config()
    qsq = cfg['q2_common'][Data.q2_bin]

    rdf = rdf.Filter(qsq, 'q2')
    rdf = rdf.Filter(Data.mva_cmb, 'MVA cmb')
    rdf = rdf.Filter(Data.mva_prc, 'MVA prc')
    rdf = rdf.Filter(Data.mss_mm, 'mass')

    log.info(f'Using dataframe for {Data.q2_bin} q2 bin')
    log.info(f'{"Mass":<20}{Data.mss_mm:<50}')
    rep = rdf.Report()
    rep.Print()

    return rdf
# ---------------------------------
def _initialize() -> None:
    os.makedirs(Data.cache_dir, exist_ok=True)
# ---------------------------------
def _get_signal() -> FitComponent:
    rdf   = _get_rdf(is_mc= True)
    rdf   = rdf.Define('weights', '1')

    l_pdf = ['cbr'] + 2 * ['cbl']
    l_shr = ['mu', 'sg']
    mod   = ModelFactory(obs = Data.obs, l_pdf = l_pdf, l_shared=l_shr)
    pdf   = mod.get_pdf()
    obj   = FitComponent(cfg=Data.mc_cfg, rdf=rdf, pdf=pdf)

    return obj
# ---------------------------------
def _get_combinatorial() -> FitComponent:
    mc_cfg= copy.deepcopy(Data.mc_cfg)
    del mc_cfg['fitting']

    mc_cfg['name'] = 'comb'

    mod   = ModelFactory(obs = Data.obs, l_pdf = ['exp'], l_shared=[])
    pdf   = mod.get_pdf()
    obj   = FitComponent(cfg=mc_cfg, rdf=None, pdf=pdf)

    return obj
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

    obj = DTFitter(rdf = rdf, components = [cmp_cmb, cmp_sig], cfg = Data.dt_cfg)
    obj.fit()
# ---------------------------------
if __name__ == '__main__':
    main()
