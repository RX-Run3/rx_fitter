'''
Script used to validate PDFs needed to fit combinatorial
'''
import argparse
from importlib.resources    import files

import yaml
import zfit
import matplotlib.pyplot as plt

from ROOT                   import RDataFrame
from zfit.core.data         import Data      as zdata
from zfit.core.basepdf      import BasePDF   as zpdf
from dmu.logging.log_store  import LogStore
from dmu.stats.fitter       import Fitter
from dmu.stats.zfit_plotter import ZFitPlotter
from rx_data.rdf_getter     import RDFGetter
from rx_selection           import selection as sel
from rx_fitter              import models

log=LogStore.add_logger('rx_fitter:validate_cmb')
# --------------------------------
class Data:
    '''
    Dataclass
    '''
    obs = zfit.Space('mass', limits=(4500, 6000))

    cfg    : dict
    q2bin  : str
    config : str
    sample : str
    trigger: str
# --------------------------------
def _parse_args() -> None:
    parser = argparse.ArgumentParser(description='Used to perform fits to validate choice of PDF for combinatorial')
    parser.add_argument('-q', '--q2bin'  , type=str, help='Q2bin', choices=['low', 'central', 'high'], required=True)
    parser.add_argument('-c', '--config' , type=str, help='Name of config file', required=True)
    parser.add_argument('-s', '--sample' , type=str, help='Name of sample'     , required=True)
    parser.add_argument('-t', '--trigger', type=str, help='Name of trigger'    , required=True)
    args = parser.parse_args()

    Data.q2bin  = args.q2bin
    Data.config = args.config
    Data.sample = args.sample
    Data.trigger= args.trigger
# --------------------------------
def _apply_selection(rdf : RDataFrame) -> RDataFrame:
    d_sel = sel.selection(project='RK', trigger=Data.trigger, q2bin=Data.q2bin, process=Data.sample)
    if 'selection' in Data.cfg:
        log.info('Updating selection')
        d_cut = Data.cfg['selection']
        d_sel.update(d_cut)

    for cut_name, cut_expr in d_sel.items():
        rdf = rdf.Filter(cut_expr, cut_name)

    return rdf
# --------------------------------
def _get_data() ->  zdata:
    gtr = RDFGetter(sample=Data.sample, trigger=Data.trigger)
    gtr.initialize()
    rdf = gtr.get_rdf()
    rdf = _apply_selection(rdf)

    var = Data.cfg['input']['observable']
    arr_mass = rdf.AsNumpy([var])[var]
    data     = zfit.Data.from_numpy(obs=Data.obs, array=arr_mass)

    return data
# --------------------------------
def _fit(pdf : zpdf, data : zdata) -> None:
    obj = Fitter(pdf, data)
    res = obj.fit()

    return res
# --------------------------------
def _plot(pdf : zpdf, data : zdata) -> None:
    obj= ZFitPlotter(data=data, model=pdf)
    obj.plot(nbins=50)
    plt.savefig('fit.png')
# --------------------------------
def _initialize() -> None:
    conf_path = files('rx_fitter_data').joinpath(f'combinatorial/{Data.config}.yaml')
    with open(conf_path, encoding='utf-8') as ifile:
        Data.cfg = yaml.safe_load(ifile)
# --------------------------------
def main():
    '''
    Start here
    '''
    _parse_args()
    _initialize()

    pdf  = models.get_pdf(obs=Data.obs, name='HypExp')
    data = _get_data()
    res  = _fit(pdf, data)

    print(res)
    _plot(pdf, data)
# --------------------------------
if __name__ == '__main__':
    main()
