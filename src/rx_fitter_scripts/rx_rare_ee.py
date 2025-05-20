'''
Script used to fit the rare mode
'''
import os
import inspect
import argparse
from typing              import Union
from importlib.resources import files

import yaml
import pandas            as pnd
import matplotlib.pyplot as plt

from dmu.stats.zfit              import zfit
from dmu.generic                 import hashing
from dmu.generic                 import utilities  as gut
from dmu.logging.log_store       import LogStore
from dmu.rdataframe              import utilities  as rut
from dmu.stats.zfit_plotter      import ZFitPlotter
from dmu.stats.fitter            import Fitter
from dmu.stats.fit_stats         import FitStats
from dmu.stats                   import utilities  as stat_utilities

from zfit.core.data              import Data       as zdata
from zfit.core.basepdf           import BasePDF    as zpdf
from zfit.core.interfaces        import ZfitSpace  as zobs
from zfit.core.parameter         import Parameter  as zpar
from zfit.result                 import FitResult  as zres

from rx_misid.misid_pdf          import MisIdPdf
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
    hsh          : str
    sample       : str
    mass         : str
    l_nbrem      : list[int]
    d_custom_sel : dict[str,str]
    d_total_sel  : dict[str,str]
    minx         : int
    maxx         : int
    nbins        : int
    mid_vers     : str
    obs          : zobs
    l_pdf        : list[zpdf]

    cache_dir    : str        = '/tmp/rx_fitter/cache'
    gut.TIMER_ON : bool       = True
    log_level    : int        = 20
    version      : str        = 'v1'
    nsig         : zpar       = zfit.Parameter('nsig', 0, 0, 10_000)
    # --------------------------------
    @staticmethod
    def is_hashable(obj, name : str) -> bool:
        '''
        Will take an object and its name, the the object are attributes of the Data class.
        It will return a bool indicating if the object is hashable
        '''
        if name.startswith('__'):
            return False

        if inspect.isroutine(obj) or inspect.ismethoddescriptor(obj):
            return False

        tp_obj=type(obj)
        tp_str=str(tp_obj)

        if 'zfit' in tp_str:
            return False

        log.debug(f'Will use {name} for hashing')

        return True
    # --------------------------------
    @classmethod
    def get_hash(cls) -> str:
        '''
        Creates hash from class attributes and returns it as a string
        '''
        data = {
                k : v
                for k, v in vars(cls).items()
                if cls.is_hashable(obj=v, name=k)}

        val = hashing.hash_object(data)

        return val
# --------------------------
def _update_selection_with_brem() -> None:
    if set(Data.l_nbrem) == {0, 1, 2}:
        return

    l_brem_cut = [ f'(nbrem == {brem})' for brem in Data.l_nbrem ]
    brem_cut   = ' || '.join(l_brem_cut)

    Data.d_sel['nbrem'] = brem_cut
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
    log.info(30 * '-')
    log.info('Adding combinatorial')
    log.info(30 * '-')

    cfg  = _load_config(component = 'combinatorial')
    pdf  = cmp.get_cb(obs=Data.obs, q2bin=Data.q2bin, cfg=cfg)
    ncmb = zfit.Parameter('ncmb', 0, 0, 20_000)
    pdf  = pdf.create_extended(ncmb)

    Data.l_pdf.append(pdf)
# --------------------------
def _add_pdf_prc(sample : str) -> None:
    log.info(30 * '-')
    log.info(f'Adding: {sample}')
    log.info(30 * '-')
    cfg                   = _load_config(component='bxhsee')
    cfg['input']['q2bin'] = Data.q2bin
    pdf                   = cmp.get_kde(obs=Data.obs, sample=sample, cfg=cfg)

    if pdf is None:
        log.warning(f'No PDF found for PRec sample {sample}, skipping')
        return

    scale= zfit.Parameter(f's{sample}', 0, 0, 10)
    nprc = zfit.ComposedParameter(f'n{sample}', lambda x : x['nsig'] * x['scale'], params={'nsig' : Data.nsig, 'scale' : scale})
    pdf.set_yield(nprc)

    Data.l_pdf.append(pdf)
# --------------------------
def _add_pdf_leak(sample : str) -> None:
    log.info(30 * '-')
    log.info(f'Adding: {sample}')
    log.info(30 * '-')
    cfg                   = _load_config(component='ccbar_leak')
    cfg['input']['q2bin'] = Data.q2bin
    pdf                   = cmp.get_kde(obs=Data.obs, sample=sample, cfg=cfg)

    if pdf is None:
        log.warning(f'No PDF found for leakage sample {sample}, skipping')
        return

    nleak = zfit.Parameter(f'n{sample}', 0, 0, 10_000)
    pdf.set_yield(nleak)

    Data.l_pdf.append(pdf)
# --------------------------
def _add_pdf_sig() -> None:
    log.info(30 * '-')
    log.info('Adding signal component')
    log.info(30 * '-')

    cfg  = _load_config(component='signal')
    cfg['input']['q2bin'] = Data.q2bin

    # l_nbrem parameter will not determine a selection of data
    # It will select which components to use to build the PDF
    pdf  = cmp.get_mc_reparametrized(obs=Data.obs, component_name='Signal', cfg=cfg, l_nbrem=Data.l_nbrem)
    pdf  = pdf.create_extended(Data.nsig)

    Data.l_pdf.append(pdf)
# --------------------------
def _add_pdf_mid() -> None:
    obj = MisIdPdf(obs=Data.obs, q2bin=Data.q2bin, version=Data.mid_vers)
    pdf = obj.get_pdf()

    Data.l_pdf.append(pdf)
# --------------------------
@gut.timeit
def _get_pdf() -> zpdf:
    _add_pdf_cmb()
    _add_pdf_mid()
    _add_pdf_prc(sample='Bu_Kstee_Kpi0_eq_btosllball05_DPC')
    _add_pdf_prc(sample='Bd_Kstee_eq_btosllball05_DPC')
    _add_pdf_prc(sample='Bs_phiee_eq_Ball_DPC')

    _add_pdf_leak(sample='Bu_JpsiK_ee_eq_DPC')
    _add_pdf_leak(sample='Bu_psi2SK_ee_eq_DPC')

    _add_pdf_sig()

    pdf = zfit.pdf.SumPDF(Data.l_pdf)

    return pdf
# --------------------------
@gut.timeit
def _get_data() -> zdata:
    log.info(20 * '-')
    log.info('Getting data')
    log.info(20 * '-')
    d_sel           = sel.selection(trigger=Data.trigger, q2bin=Data.q2bin, process=Data.sample)
    hsh             = hashing.hash_object([d_sel, Data.sample, Data.trigger, Data.mass])
    data_cache_path = f'{Data.cache_dir}/{hsh}.json'
    if os.path.isfile(data_cache_path):
        log.warning(f'Caching data from: {data_cache_path}')
        l_mass   = gut.load_json(data_cache_path)
        arr_mass = numpy.array(l_mass)
        data = zfit.Data.from_numpy(obs=Data.obs, array=arr_mass)

        return data

    gtr = RDFGetter(sample=Data.sample, trigger=Data.trigger)
    rdf = gtr.get_rdf()

    for cut_name, cut_expr in d_sel.items():
        log.info(f'{cut_name:<20}{cut_expr}')
        rdf = rdf.Filter(cut_expr, cut_name)

    rep = rdf.Report()
    rep.Print()

    mass = Data.mass.replace('_smr_', '_') # Real data is not smeared
    log.info(f'Using mass {mass} for real data')

    arr_mass = rdf.AsNumpy([mass])[mass]
    log.info(f'Caching data to: {data_cache_path}')
    gut.dump_json(arr_mass.tolist(), data_cache_path)

    data = zfit.Data.from_numpy(obs=Data.obs, array=arr_mass)

    return data
# --------------------------
def _get_constraints(pdf : zpdf) -> dict[str,tuple[float,float]]:
    s_par  = pdf.get_params()
    l_name = [par.name for par in s_par]

    obj    = ConstraintReader(parameters = l_name, q2bin=Data.q2bin)
    d_cns  = obj.get_constraints()

    return d_cns
# --------------------------
def _get_sensitivity() -> float:
    '''
    Returns fit sensitivity in %
    '''
    obj = FitStats(fit_dir=Data.fit_dir)
    val = obj.get_value(name='nsig', kind = 'value')
    err = obj.get_value(name='nsig', kind = 'error')

    return 100 * err/val
# --------------------------
def _get_text(data : zdata) -> str:
    arr_mass = data.to_numpy()
    mask     = (arr_mass > Data.minx) & (arr_mass < Data.maxx)
    arr_mass = arr_mass[mask]
    nentries = len(arr_mass)

    text = ''
    for name, cut in Data.d_sel.items():
        # Brem cut is too long and will be in title anyway
        if name == 'nbrem':
            continue

        text += f'\n{name}: {cut}'

    sensitivity = _get_sensitivity()
    title       = f'$\\delta={sensitivity:.2f}$%; Entries={nentries:.0f}; Brem:{Data.l_nbrem}'

    return text, title
# --------------------------
def _set_hash(cfg : dict) -> None:
    data_hash = Data.get_hash()
    Data.hsh  = hashing.hash_object([data_hash, cfg])
# --------------------------
def _initialize() -> None:
    LogStore.set_level('rx_fitter:constraint_reader' , Data.log_level)
    LogStore.set_level('rx_fitter:components'        , Data.log_level)
    LogStore.set_level('rx_fitter:rx_rare_ee'        , Data.log_level)
    LogStore.set_level('rx_calibration:fit_component', Data.log_level)

    cfg = _load_config(component='data')
    _initialize_settings(cfg=cfg)

    sel.set_custom_selection(d_cut=Data.d_sel)

    ana_dir = os.environ['ANADIR']
    sample  = Data.sample.replace('*', 'p')

    _set_hash(cfg=cfg)
    Data.fit_dir = f'{ana_dir}/fits/{sample}/{Data.trigger}/{Data.version}/{Data.q2bin}/{Data.hsh}'
    gut.dump_json(cfg, f'{Data.fit_dir}/config.yaml')
    Data.l_pdf   = []
# --------------------------
def _initialize_settings(cfg : dict) -> None:
    Data.l_nbrem = cfg['nbrem'][Data.q2bin]
    Data.mass    = cfg['input']['observable']
    Data.minx    = cfg['input']['minx']
    Data.maxx    = cfg['input']['maxx']
    Data.nbins   = cfg['input']['nbins']
    Data.trigger = cfg['input']['trigger']
    Data.sample  = cfg['input']['sample']
    if 'selection' in cfg['input']:
        Data.d_sel = cfg['input']['selection']
    else:
        Data.d_sel = {}

    Data.mid_vers = cfg['fitting']['misid']['version']

    _update_selection_with_brem()

    Data.obs = zfit.Space(Data.mass, limits=(Data.minx, Data.maxx))
# --------------------------
@gut.timeit
def _fit(pdf : zpdf, data : zdata, constraints : dict[str,tuple[float,float]]) -> Union[zres,None]:
    cfg = {
            'constraints' : constraints,
            }

    if Data.dry_run:
        log.warning('Running dry run')
        return None

    obj = Fitter(pdf, data)
    res = obj.fit(cfg=cfg)

    return res
# --------------------------
def _plot_fit(data : zdata, pdf : zpdf):
    d_leg = {
            'SumPDF_ext'                        : 'Signal',
            'exp_1_ext'                         : 'Combinatorial',
            'Bu_JpsiK_ee_eq_DPC'                : r'$B^+\to J/\psi(\to ee) K^+$',
            'Bu_Kstee_Kpi0_eq_btosllball05_DPC' : r'$B^+\to K^{*+}ee$',
            'Bd_Kstee_eq_btosllball05_DPC'      : r'$B^0\to K^{*0}ee$',
            'Bs_phiee_eq_Ball_DPC'              : r'$B_s\to \phi(\to KK)ee$',
            'suj_1_ext'                         : 'Combinatorial',
            'hypexp_1_ext'                      : 'Combinatorial',
            }

    text, title = _get_text(data)
    obj   = ZFitPlotter(data=data, model=pdf)
    obj.plot(
            nbins   =Data.nbins,
            d_leg   =d_leg,
            stacked =True,
            title   =title,
            leg_loc ='upper right',
            ext_text=text)

    obj.axs[1].set_xlabel(Data.mass)
    obj.axs[0].axvline(x=5280, linestyle='--', color='gray', label='$B^+$')

    obj.axs[1].set_ylim([-5, +5])
    obj.axs[1].plot([Data.minx, Data.maxx], [+3, +3], linestyle='--', color='red')
    obj.axs[1].plot([Data.minx, Data.maxx], [-3, -3], linestyle='--', color='red')

    plt.savefig(f'{Data.fit_dir}/fit.png')
    plt.close()
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
    if fit_result is None:
        return

    stat_utilities.save_fit(
            data   =data,
            model  =pdf,
            res    =fit_result,
            fit_dir=Data.fit_dir,
            d_const=d_cns)

    _plot_fit(data=data, pdf=pdf)
# --------------------------
if __name__ == '__main__':
    main()
