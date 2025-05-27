'''
Script used to validate PDFs needed to fit combinatorial
'''
import os
import re
import argparse
from importlib.resources    import files

import yaml
import matplotlib.pyplot as plt

from dmu.logging.log_store  import LogStore
from dmu.stats.zfit         import zfit
from dmu.stats.fitter       import Fitter
from dmu.stats.zfit_plotter import ZFitPlotter

from ROOT                   import RDataFrame
from zfit.core.data         import Data      as zdata
from zfit.core.basepdf      import BasePDF   as zpdf
from zfit.core.interfaces   import ZfitSpace as zobs
from rx_data.rdf_getter     import RDFGetter
from rx_selection           import selection as sel
from rx_fitter              import models

log=LogStore.add_logger('rx_fitter:validate_cmb')
# --------------------------------
class Data:
    '''
    Dataclass
    '''
    minx   : float
    maxx   : float
    mass   : str

    wp_cmb : float
    wp_prc : float

    obs    : zobs
    cfg    : dict
    out_dir: str
    q2bin  : str
    q2_kind: str
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
    parser.add_argument('-q', '--q2bin'  , type=str, help='Q2bin'         , choices=['low', 'central', 'high'], required=True)
    parser.add_argument('-m', '--model'  , type=str, help='Fitting model' , choices=['HypExp', 'ModExp', 'Exp', 'Pol2', 'Pol3', 'SUJohnson'], default='SUJohnson')
    parser.add_argument('-k', '--q2_kind', type=str, help='Kind of q2 cut', choices=['def', 'dtf', 'trk'])
    parser.add_argument('-c', '--config' , type=str, help='Name of config file'                           , default='validation')
    parser.add_argument('-s', '--sample' , type=str, help='Name of sample'                                , default='DATA*')
    parser.add_argument('-t', '--trigger', type=str, help='Name of trigger'                               , default='Hlt2RD_BuToKpEE_SameSign_MVA')
    parser.add_argument('-i', '--initial', type=int, help='Index of initial fit'                          , default=0)
    parser.add_argument('-f', '--final'  , type=int, help='Index of final fit, if not passed, will do all', default=1000)
    parser.add_argument('-n', '--ntries' , type=int, help='Maximum number of tries, default 1'            , default=1)
    parser.add_argument('-w', '--wpoint' , nargs=2 , help='Array with two working points, combinatorial and prec')
    args = parser.parse_args()

    if args.wpoint is not None:
        [ cmb, prc] = args.wpoint
        Data.wp_cmb = float(cmb)
        Data.wp_prc = float(prc)

    Data.q2bin  = args.q2bin
    Data.model  = args.model
    Data.q2_kind= args.q2_kind
    Data.config = args.config
    Data.sample = args.sample
    Data.trigger= args.trigger
    Data.initial= args.initial
    Data.final  = args.final
    Data.ntries = args.ntries
# --------------------------------
def _apply_selection(rdf : RDataFrame) -> RDataFrame:
    d_sel = sel.selection(trigger=Data.trigger, q2bin=Data.q2bin, process=Data.sample)
    for cut_name, cut_expr in d_sel.items():
        rdf = rdf.Filter(cut_expr, cut_name)

    return rdf
# --------------------------------
def _suffix_from_name(name : str) -> str:
    name = name.replace('||', 'OR')
    name = name.replace(' ' ,  '_')
    name = name.replace('<' , 'lt')
    name = name.replace('>' , 'gt')
    name = name.replace('=' , 'eq')
    name = name.replace('{' ,  '_')
    name = name.replace('}' ,  '_')
    name = name.replace('.' ,  'p')
    name = name.replace('$' ,  '_')
    name = name.replace('&&','and')
    name = re.sub(r'_+', '_', name)

    return name
# --------------------------------
def _get_rdf() ->  zdata:
    gtr = RDFGetter(sample=Data.sample, trigger=Data.trigger)
    rdf = gtr.get_rdf()
    rdf = _apply_selection(rdf)

    return rdf
# --------------------------------
def _data_from_rdf(rdf : RDataFrame, cut : str) ->  zdata:
    rdf      = rdf.Filter(cut)
    arr_mass = rdf.AsNumpy([Data.mass])[Data.mass]
    data     = zfit.Data.from_numpy(obs=Data.obs, array=arr_mass)

    return data
# --------------------------------
def _fit(pdf : zpdf, data : zdata) -> None:
    fit_cfg = Data.cfg['fitting']

    obj = Fitter(pdf, data)
    res = obj.fit(cfg=fit_cfg)

    return res
# --------------------------------
def _get_out_dir() -> str:
    ana_dir = os.environ['ANADIR']
    out_dir = Data.cfg['output']['path']

    if hasattr(Data, 'q2_kind'):
        out_dir = f'{ana_dir}/{Data.q2_kind}/{out_dir}'
    else:
        out_dir = f'{ana_dir}/{out_dir}'

    os.makedirs(out_dir, exist_ok=True)

    return out_dir
# --------------------------------
def _plot(pdf : zpdf, data : zdata, name : str) -> None:
    suffix   = _suffix_from_name(name)
    nentries = data.value().shape[0]
    ext_text = f'Entries={nentries}\n{Data.sample}\n{Data.trigger}'

    obj= ZFitPlotter(data=data, model=pdf)
    rng= Data.cfg['fitting']['ranges']
    obj.plot(nbins=50, ranges=rng, title=name, ext_text=ext_text, d_leg={'ZPDF' : Data.model})

    obj.axs[0].axvline(x=5280, linestyle='--', color='gray', label='$B^+$')

    obj.axs[1].set_ylim([-5, +5])
    obj.axs[1].plot([Data.minx, Data.maxx], [+3, +3], linestyle='--', color='red')
    obj.axs[1].plot([Data.minx, Data.maxx], [-3, -3], linestyle='--', color='red')

    plot_path = f'{Data.out_dir}/fit_{suffix}.png'
    log.info(f'Saving to: {plot_path}')
    plt.savefig(plot_path)
# --------------------------------
def _override_q2(cuts : dict[str,str]) -> dict[str,str]:
    if Data.q2_kind is None:
        log.debug('Not overriding q2 cut')
        return cuts

    cut = Data.cfg['q2_kind'][Data.q2_kind]
    log.info(f'Using q2 cut: {cut}')
    cuts['q2'] = cut

    return cuts
# --------------------------------
def _initialize() -> None:
    conf_path = files('rx_fitter_data').joinpath(f'combinatorial/{Data.config}.yaml')
    with open(conf_path, encoding='utf-8') as ifile:
        Data.cfg = yaml.safe_load(ifile)

    if 'selection' in Data.cfg:
        log.info('Updating selection')
        d_cut = Data.cfg['selection']
        d_cut = _override_q2(cuts=d_cut)

        sel.set_custom_selection(d_cut = d_cut)

    Data.minx= Data.cfg['model']['observable']['minx']
    Data.maxx= Data.cfg['model']['observable']['maxx']
    Data.mass= Data.cfg['model']['observable']['name']


    Data.out_dir = _get_out_dir()
    Data.obs     = zfit.Space(Data.mass, limits=(Data.minx, Data.maxx))
# --------------------------------
def _skip_fit(index : int) -> bool:
    if Data.initial <= index <= Data.final:
        return False

    return True
# --------------------------------
def _get_cutflow() -> dict[str,str]:
    if hasattr(Data, 'wp_cmb') and hasattr(Data, 'wp_prc'):
        key = f'$BDT_{{cmb}} > {Data.wp_cmb:.2f}$ && $BDT_{{prc}} > {Data.wp_prc:.2f}$'
        cut = f'     mva_cmb > {Data.wp_cmb:.2f}  &&      mva_prc > {Data.wp_prc:.2f}'
        log.warning(f'Overriding WP with: {cut}')

        return {key : cut}

    log.debug('Picking up cutflow from YAML')
    return Data.cfg['cutflow']
# --------------------------------
def main():
    '''
    Start here
    '''
    _parse_args()
    _initialize()

    pdf  = models.get_pdf(obs=Data.obs, name=Data.model)
    rdf  = _get_rdf()

    d_cutflow = _get_cutflow()

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
