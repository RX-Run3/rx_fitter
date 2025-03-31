'''
Script used to validate PDFs needed to fit combinatorial
'''
import os
import re
import argparse
from importlib.resources    import files

import yaml
import ROOT
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
    minx = 4500
    maxx = 6000
    obs  = zfit.Space('mass', limits=(minx, maxx))

    cfg    : dict
    q2bin  : str
    model  : str
    config : str
    sample : str
    trigger: str
    initial: int
    final  : int
    ntries : int
# --------------------------------
def _parse_args() -> None:
    parser = argparse.ArgumentParser(description='Used to perform fits to validate choice of PDF for combinatorial')
    parser.add_argument('-q', '--q2bin'  , type=str, help='Q2bin', choices=['low', 'central', 'high'], required=True)
    parser.add_argument('-m', '--model'  , type=str, help='Fitting model', choices=['HypExp', 'ModExp', 'Exp', 'Pol2', 'Pol3'], required=True)
    parser.add_argument('-c', '--config' , type=str, help='Name of config file', required=True)
    parser.add_argument('-s', '--sample' , type=str, help='Name of sample'     , required=True)
    parser.add_argument('-t', '--trigger', type=str, help='Name of trigger'    , required=True)
    parser.add_argument('-i', '--initial', type=int, help='Index of initial fit', default=0)
    parser.add_argument('-f', '--final'  , type=int, help='Index of final fit, if not passed, will do all', default=1000)
    parser.add_argument('-n', '--ntries' , type=int, help='Maximum number of tries, default 1', default=1)
    args = parser.parse_args()

    Data.q2bin  = args.q2bin
    Data.model  = args.model
    Data.config = args.config
    Data.sample = args.sample
    Data.trigger= args.trigger
    Data.initial= args.initial
    Data.final  = args.final
    Data.ntries = args.ntries
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
def _suffix_from_name(name : str) -> str:
    name = name.replace(' ',  '_')
    name = name.replace('<', 'lt')
    name = name.replace('>', 'gt')
    name = name.replace('=', 'eq')
    name = name.replace('{', '_')
    name = name.replace('}', '_')
    name = name.replace('.', 'p')
    name = name.replace('$', '_')
    name = name.replace('&&', 'and')
    name = re.sub(r'_+', '_', name)

    return name
# --------------------------------
def _get_rdf() ->  zdata:
    gtr = RDFGetter(sample=Data.sample, trigger=Data.trigger)
    gtr.initialize()
    rdf = gtr.get_rdf()
    rdf = _apply_selection(rdf)

    return rdf
# --------------------------------
def _data_from_rdf(rdf : RDataFrame, cut : str) ->  zdata:
    rdf      = rdf.Filter(cut)
    var      = Data.cfg['input']['observable']
    arr_mass = rdf.AsNumpy([var])[var]
    data     = zfit.Data.from_numpy(obs=Data.obs, array=arr_mass)

    return data
# --------------------------------
def _fit(pdf : zpdf, data : zdata) -> None:
    d_retry = {'ntries' : Data.ntries, 'pvalue_thresh' : 0.05, 'ignore_status' : False}

    obj = Fitter(pdf, data)
    res = obj.fit(cfg={'strategy' : {'retry' : d_retry}})

    return res
# --------------------------------
def _get_out_dir() -> str:
    fit_dir = os.environ['FITDIR']
    plt_dir = f'{fit_dir}/{Data.sample}/{Data.trigger}/{Data.q2bin}/{Data.model}'
    plt_dir = plt_dir.replace('*', 'p')

    os.makedirs(plt_dir, exist_ok=True)

    return plt_dir
# --------------------------------
def _plot(pdf : zpdf, data : zdata, name : str) -> None:
    suffix   = _suffix_from_name(name)
    out_dir  = _get_out_dir()
    nentries = data.value().shape[0]
    ext_text = f'Entries={nentries}\n{Data.sample}\n{Data.trigger}'

    obj= ZFitPlotter(data=data, model=pdf)
    obj.plot(nbins=50, title=name, ext_text=ext_text, d_leg={'ZPDF' : Data.model})

    obj.axs[1].set_ylim([-5, +5])
    obj.axs[1].plot([Data.minx, Data.maxx], [+3, +3], linestyle='--', color='red')
    obj.axs[1].plot([Data.minx, Data.maxx], [-3, -3], linestyle='--', color='red')

    plt.savefig(f'{out_dir}/fit_{suffix}.png')
# --------------------------------
def _initialize() -> None:
    conf_path = files('rx_fitter_data').joinpath(f'combinatorial/{Data.config}.yaml')
    with open(conf_path, encoding='utf-8') as ifile:
        Data.cfg = yaml.safe_load(ifile)
# --------------------------------
def _skip_fit(index : int) -> bool:
    if Data.initial <= index <= Data.final:
        return False

    return True
# --------------------------------
def main():
    '''
    Start here
    '''
    _parse_args()
    _initialize()

    pdf  = models.get_pdf(obs=Data.obs, name=Data.model)
    rdf  = _get_rdf()

    d_cutflow = Data.cfg['cutflow']

    index = 0
    for name, cut in d_cutflow.items():
        if _skip_fit(index):
            log.info(f'Skipping {name}/{index}')
            index += 1
            continue

        log.info(f'Fitting {name}/{index}')
        data = _data_from_rdf(rdf, cut)
        _fit(pdf, data)

        _plot(pdf, data, name)

        index += 1
# --------------------------------
if __name__ == '__main__':
    main()
