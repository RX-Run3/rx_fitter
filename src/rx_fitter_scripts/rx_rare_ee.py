'''
Script used to fit the rare mode
'''
import os
import argparse
from importlib.resources import files

import yaml
import ROOT
import numpy
import zfit

from zfit.core.data              import Data       as zdata
from zfit.core.basepdf           import BasePDF    as zpdf
from zfit.core.interfaces        import ZfitSpace  as zobs
from zfit.core.parameter         import Parameter  as zpar
from zfit.result                 import FitResult  as zres

from dmu.generic                 import hashing
from dmu.generic                 import utilities  as gut
from dmu.logging.log_store       import LogStore
from dmu.stats.zfit_plotter      import ZFitPlotter
from dmu.stats.fitter            import Fitter
from dmu.stats                   import utilities  as stat_utilities

from rx_data.rdf_getter          import RDFGetter
from rx_selection                import selection  as sel
from rx_fitter                   import components as cmp
from rx_fitter.constraint_reader import ConstraintReader

log=LogStore.add_logger('rx_fitter:rx_rare_ee')
# --------------------------
class Data:
    '''
    Data class
    '''
    dry_run      : bool
    q2bin        : str
    trigger      : str
    sample       : str
    mass         : str
    l_nbrem      : list[int]
    minx         : int
    maxx         : int
    obs          : zobs

    brem_def     : str        = 'int(L1_HASBREMADDED_brem_track_2) + int(L2_HASBREMADDED_brem_track_2)'
    cache_dir    : str        = '/tmp/rx_fitter/cache'
    gut.TIMER_ON : bool       = True
    log_level    : int        = 20
    version      : str        = 'v1'
    nsig         : zpar       = zfit.Parameter('nsig', 0, 0, 10_000)
    l_pdf        : list[zpdf] = []
# --------------------------
def _initialize_settings(cfg : dict) -> None:
    Data.l_nbrem = cfg['nbrem'][Data.q2bin]
    Data.mass    = cfg['input']['observable']
    Data.minx    = cfg['input']['minx']
    Data.maxx    = cfg['input']['maxx']
    Data.trigger = cfg['input']['trigger']
    Data.sample  = cfg['input']['sample']

    Data.obs     = zfit.Space(Data.mass, limits=(Data.minx, Data.maxx))
# --------------------------------------------------------------
def _parse_args():
    parser = argparse.ArgumentParser(description='Script used to fit rare mode electron channel data for RK')
    parser.add_argument('-q', '--q2bin'  , type=str, help='q2 bin', required=True, choices=['low', 'central', 'high'])
    parser.add_argument('-l', '--loglv'  , type=int, help='Logging level', default=Data.log_level, choices=[10, 20, 30])
    parser.add_argument('-d', '--dry_run', action='store_true', help='If used, will skip fit')
    args = parser.parse_args()

    Data.q2bin     = args.q2bin
    Data.dry_run   = args.dry_run
    Data.log_level = args.loglv
# --------------------------------------------------------------
def _load_config(component : str) -> dict:
    cfg_path = files('rx_fitter_data').joinpath(f'rare_fit/{Data.version}/rk_ee/{component}.yaml')
    with open(cfg_path, encoding='utf-8') as ifile:
        cfg = yaml.safe_load(ifile)

    return cfg
# --------------------------
def _add_pdf_cmb() -> None:
    cfg  = _load_config(component = 'combinatorial')
    pdf  = cmp.get_cb(obs=Data.obs, q2bin=Data.q2bin, cfg=cfg)
    ncmb = zfit.Parameter('ncmb', 0, 0, 20_000)
    pdf  = pdf.create_extended(ncmb)

    Data.l_pdf.append(pdf)
# --------------------------
def _add_pdf_prc(sample : str) -> None:
    cfg                   = _load_config(component='bxhsee')
    cfg['input']['q2bin'] = Data.q2bin
    pdf                   = cmp.get_kde(obs=Data.obs, sample=sample, l_nbrem=Data.l_nbrem, cfg=cfg)

    if pdf is None:
        log.warning(f'No data found for leakage sample {sample}, skipping')
        return

    scale= zfit.Parameter(f's{sample}', 0, 0, 10)
    nprc = zfit.ComposedParameter(f'n{sample}', lambda x : x['nsig'] * x['scale'], params={'nsig' : Data.nsig, 'scale' : scale})
    pdf.set_yield(nprc)

    Data.l_pdf.append(pdf)
# --------------------------
def _add_pdf_leak(sample : str) -> None:
    cfg                   = _load_config(component='ccbar_leak')
    cfg['input']['q2bin'] = Data.q2bin
    pdf                   = cmp.get_kde(obs=Data.obs, sample=sample, l_nbrem=Data.l_nbrem, cfg=cfg)

    if pdf is None:
        log.warning(f'No data found for leakage sample {sample}, skipping')
        return

    nleak = zfit.Parameter(f'n{sample}', 0, 0, 10_000)
    pdf.set_yield(nleak)

    Data.l_pdf.append(pdf)
# --------------------------
def _add_pdf_sig() -> None:
    cfg  = _load_config(component='signal')
    cfg['input']['q2bin'] = Data.q2bin

    pdf  = cmp.get_mc_reparametrized(obs=Data.obs, component_name='Signal', cfg=cfg, l_nbrem=Data.l_nbrem)
    pdf  = pdf.create_extended(Data.nsig)

    Data.l_pdf.append(pdf)
# --------------------------
@gut.timeit
def _get_pdf() -> zpdf:
    _add_pdf_cmb()
    _add_pdf_prc(sample='Bu_Kstee_Kpi0_eq_btosllball05_DPC')
    _add_pdf_prc(sample='Bd_Kstee_eq_btosllball05_DPC')
    _add_pdf_prc(sample='Bs_phiee_eq_Ball_DPC')

    _add_pdf_leak(sample='Bu_JpsiK_ee_eq_DPC')
    _add_pdf_leak(sample='Bu_psi2SK_ee_eq_DPC')

    _add_pdf_sig()

    pdf = zfit.pdf.SumPDF(Data.l_pdf)

    return pdf
# --------------------------
def _get_brem_cut() -> str:
    l_cut = [ f'(nbrem == {nbrem})' for nbrem in Data.l_nbrem ]
    cut   = '||'.join(l_cut)

    return cut
# --------------------------
def _get_q2_cut(d_sel : dict[str,str]) -> str:
    if Data.q2bin != 'high':
        cut = d_sel['q2']
        return cut

    return Data.high_q2_cut
# --------------------------
@gut.timeit
def _get_data() -> zdata:
    gtr = RDFGetter(sample=Data.sample, trigger=Data.trigger)
    rdf = gtr.get_rdf()
    rdf = rdf.Define('nbrem', Data.brem_def)

    d_sel        = sel.selection(project='RK', trigger=Data.trigger, q2bin=Data.q2bin, process=Data.sample)
    d_sel['q2']  = _get_q2_cut(d_sel)
    d_sel['mva'] = Data.mva_cut
    d_sel['mass']= '(1)'
    d_sel['brem']= _get_brem_cut()

    hsh             = hashing.hash_object([d_sel, Data.sample, Data.trigger, Data.mass])
    data_cache_path = f'{Data.cache_dir}/{hsh}.json'
    if os.path.isfile(data_cache_path):
        log.warning(f'Caching data from: {data_cache_path}')
        l_mass   = gut.load_json(data_cache_path)
        arr_mass = numpy.array(l_mass)
        data = zfit.Data.from_numpy(obs=Data.obs, array=arr_mass)

        return data

    for cut_name, cut_expr in d_sel.items():
        log.info(f'{cut_name:<20}{cut_expr}')
        rdf = rdf.Filter(cut_expr, cut_name)

    rep = rdf.Report()
    rep.Print()

    arr_mass = rdf.AsNumpy([Data.mass])[Data.mass]
    log.info(f'Caching data to: {data_cache_path}')
    gut.dump_json(arr_mass.tolist(), data_cache_path)

    data = zfit.Data.from_numpy(obs=Data.obs, array=arr_mass)

    return data
# --------------------------
def _get_constraints(pdf : zpdf) -> dict[str,tuple[float,float]]:
    s_par  = pdf.get_params()
    l_name = [par.name for par in s_par]

    obj    = ConstraintReader(parameters = l_name, q2bin=Data.q2bin, mva_cut = Data.mva_cut)
    d_cns  = obj.get_constraints()

    return d_cns
# --------------------------
def _get_title() -> str:
    title = f'{Data.mva_cut}; Brem:{Data.l_nbrem}'

    return title
# --------------------------
def _get_extra_text(data : zdata) -> str:
    arr_mass = data.to_numpy()
    mask     = (arr_mass > Data.minx) & (arr_mass < Data.maxx)
    arr_mass = arr_mass[mask]
    nentries = len(arr_mass)

    if Data.q2bin == 'high':
        return f'Entries={nentries:.0f}\n{Data.high_q2_trk}'

    return f'Entries={nentries:.0f}'
# --------------------------
def _initialize() -> None:
    LogStore.set_level('rx_fitter:constraint_reader', Data.log_level)

    fit_dir      = os.environ['FITDIR']
    sample       = Data.sample.replace('*', 'p')
    Data.fit_dir = f'{fit_dir}/{sample}/{Data.trigger}/{Data.version}/{Data.q2bin}'
# --------------------------
@gut.timeit
def _fit(pdf : zpdf, data : zdata, constraints : dict[str,tuple[float,float]]) -> zres:
    cfg = {
            'constraints' : constraints,
            }

    if Data.dry_run:
        log.warning('Running dry run')
    else:
        obj = Fitter(pdf, data)
        res = obj.fit(cfg=cfg)

    title    = f'Brem:{Data.l_nbrem}'
    ext_text = _get_extra_text(data)

    d_leg = {
            'SumPDF_ext'                        : 'Signal',
            'exp_1_ext'                         : 'Combinatorial',
            'Bu_JpsiK_ee_eq_DPC'                : r'$B^+\to J/\psi(\to ee) K^+$',
            'Bu_Kstee_Kpi0_eq_btosllball05_DPC' : r'$B^+\to K^{*+}ee$',
            'Bd_Kstee_eq_btosllball05_DPC'      : r'$B^0\to K^{*0}ee$',
            'Bs_phiee_eq_Ball_DPC'              : r'$B_s\to \phi(\to KK)ee$',
            'suj_1_ext'                         : 'Combinatorial',
            }

    obj   = ZFitPlotter(data=data, model=pdf)
    obj.plot(nbins=50, d_leg=d_leg, stacked=True, title=title, ext_text=ext_text)
    obj.axs[1].set_xlabel(Data.mass)
    obj.axs[0].axvline(x=5280, linestyle='--', color='gray', label='$B^+$')

    obj.axs[1].set_ylim([-5, +5])
    obj.axs[1].plot([Data.minx, Data.maxx], [+3, +3], linestyle='--', color='red')
    obj.axs[1].plot([Data.minx, Data.maxx], [-3, -3], linestyle='--', color='red')

    return res
# --------------------------
def main():
    '''
    Start here
    '''
    _parse_args()
    _initialize()

    data = _get_data()
    pdf  = _get_pdf()
    d_cns= _get_constraints(pdf)

    stat_utilities.print_pdf(pdf=pdf, d_const=d_cns, txt_path=f'{Data.fit_dir}/pre_fit.txt')
    fit_result = _fit(pdf=pdf, data=data, constraints=d_cns)

    stat_utilities.save_fit(
            data   =data,
            model  =pdf,
            res    =fit_result,
            fit_dir=Data.fit_dir,
            d_const=d_cns)
# --------------------------
if __name__ == '__main__':
    main()
