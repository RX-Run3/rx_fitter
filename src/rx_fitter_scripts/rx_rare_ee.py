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
import matplotlib.pyplot as plt
from ROOT                        import EnableImplicitMT
from zfit.core.data              import Data       as zdata
from zfit.core.basepdf           import BasePDF    as zpdf
from zfit.core.interfaces        import ZfitSpace  as zobs
from zfit.core.parameter         import Parameter  as zpar
from dmu.generic                 import hashing
from dmu.generic                 import utilities  as gut
from dmu.logging.log_store       import LogStore
from dmu.stats.zfit_plotter      import ZFitPlotter
from dmu.stats.fitter            import Fitter
from dmu.stats.utilities         import print_pdf
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
    q2bin   : str
    high_q2_trk : str  = '(q2_track > 15000000)'
    high_q2_nom : str  = '(q2       > 15500000) && (q2       < 22000000)'
    high_q2_cut : str  = f'{high_q2_trk} && {high_q2_nom}'

    l_nbrem     : list[int] = [0, 1, 2]

    brem_def    : str  = 'int(L1_HASBREMADDED_brem_track_2) + int(L2_HASBREMADDED_brem_track_2)'
    cache_dir   : str  = '/tmp/rx_fitter/cache'
    mva_cut     : str  = '(mva_cmb > 0.60) && (mva_prc > 0.80)'
    sample      : str  = 'DATA*'
    trigger     : str  = 'Hlt2RD_BuToKpEE_MVA'
    version     : str  = 'v1'
    mass        : str  = 'B_M_brem_track_2'
    minx        : int  = 4_500
    maxx        : int  = 7_000
    obs         : zobs = zfit.Space(mass, limits=(minx, maxx))
    nsig        : zpar = zfit.Parameter('nsig', 0, 0, 10_000)
    gut.TIMER_ON: bool = True

    log_level : int        = 20
    l_pdf     : list[zpdf] = []
# --------------------------------------------------------------
def _parse_args():
    parser = argparse.ArgumentParser(description='Script used to fit rare mode electron channel data for RK')
    parser.add_argument('-q', '--q2bin' , type=str, help='q2 bin', required=True, choices=['low', 'central', 'high'])
    parser.add_argument('-l', '--loglv' , type=int, help='Logging level', default=Data.log_level, choices=[10, 20, 30])
    args = parser.parse_args()

    Data.q2bin     = args.q2bin
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
def _set_logs() -> None:
    LogStore.set_level('rx_fitter:constraint_reader', Data.log_level)
# --------------------------
def _add_pdf_prc(sample : str) -> None:
    cfg                   = _load_config(component='bxhsee')
    cfg['input']['q2bin'] = Data.q2bin
    cfg['selection']      = {'bdt' : Data.mva_cut}
    if Data.q2bin == 'high':
        cfg['selection']['q2'] = Data.high_q2_cut

    pdf  = cmp.get_kde(obs=Data.obs, sample=sample, l_nbrem=Data.l_nbrem, cfg=cfg)
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
    cfg['selection']      = {'mva' : Data.mva_cut}
    if Data.q2bin == 'high':
        cfg['selection']['q2'] = Data.high_q2_cut

    pdf   = cmp.get_kde(obs=Data.obs, sample=sample, l_nbrem=Data.l_nbrem, cfg=cfg)
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
    nentries = numpy.sum(arr_mass)

    if Data.q2bin == 'high':
        return f'Entries={nentries:.0f}\n{Data.high_q2_trk}'

    return f'Entries={nentries:.0f}'
# --------------------------
@gut.timeit
def _fit(pdf : zpdf, data : zdata, constraints : dict[str,tuple[float,float]]) -> None:
    cfg = {
            'constraints' : constraints,
            }

    obj = Fitter(pdf, data)
    res = obj.fit(cfg=cfg)
    log.info(res)

    title    = _get_title()
    ext_text = _get_extra_text(data)

    obj   = ZFitPlotter(data=data, model=pdf)
    obj.plot(nbins=50, stacked=True, title=title, ext_text=ext_text)
    obj.axs[1].set_xlabel(Data.mass)
    obj.axs[0].axvline(x=5280, linestyle='--', color='gray', label='$B^+$')

    obj.axs[1].set_ylim([-5, +5])
    obj.axs[1].plot([Data.minx, Data.maxx], [+3, +3], linestyle='--', color='red')
    obj.axs[1].plot([Data.minx, Data.maxx], [-3, -3], linestyle='--', color='red')

    plt.savefig(f'fit_{Data.q2bin}.png')
    plt.close()
# --------------------------
def main():
    '''
    Start here
    '''
    _parse_args()
    _set_logs()
    EnableImplicitMT(10)

    data = _get_data()
    pdf  = _get_pdf()
    d_cns= _get_constraints(pdf)
    print_pdf(pdf=pdf, d_const=d_cns, txt_path=f'./pdf_{Data.q2bin}.txt')

    _fit(pdf=pdf, data=data, constraints=d_cns)
# --------------------------
if __name__ == '__main__':
    main()
