'''
Module with tests for the PRec class
'''

import os
import mplhep
import pytest
import matplotlib.pyplot as plt

from dmu.stats.fitter       import Fitter
from dmu.stats              import utilities as sut
from dmu.stats.zfit         import zfit
from dmu.stats.zfit_plotter import ZFitPlotter
from dmu.logging.log_store  import LogStore
from rx_selection           import selection as sel
from rx_data.rdf_getter     import RDFGetter
from rx_fitter.prec         import PRec

log=LogStore.add_logger('rx_fitter:test_prec')
#-----------------------------------------------
class Data:
    '''
    data class
    '''
    out_dir = '/tmp/tests/rx_fitter/prec'
#-----------------------------------------------
@pytest.fixture(scope='session', autouse=True)
def _initialize():
    LogStore.set_level('rx_fitter:inclusive_decay_weights' , 10)
    LogStore.set_level('rx_fitter:inclusive_sample_weights', 10)
    LogStore.set_level('rx_fitter:prec'                    , 10)

    plt.style.use(mplhep.style.LHCb2)

    RDFGetter.samples = {
        'main'       : '/home/acampove/external_ssd/Data/samples/main.yaml',
        'mva'        : '/home/acampove/external_ssd/Data/samples/mva.yaml',
        'cascade'    : '/home/acampove/external_ssd/Data/samples/cascade.yaml',
        'jpsi_misid' : '/home/acampove/external_ssd/Data/samples/jpsi_misid.yaml',
        'hop'        : '/home/acampove/external_ssd/Data/samples/hop.yaml',
        }
#-----------------------------------------------
def _plot_weight(arr_wgt, label : str, linestyle : str):
    plt.hist(arr_wgt, bins=30, label=label, histtype='step', linestyle=linestyle)
#-----------------------------------------------
def _plot_pdf(pdf, test : str, name : str, maxy : str, title : str):
    arr_mass = pdf.arr_mass
    arr_wgt  = pdf.arr_wgt
    arr_sam  = pdf.arr_sam
    arr_dec  = pdf.arr_dec

    obj = ZFitPlotter(data=arr_mass, model=pdf, weights=arr_wgt)
    obj.plot(stacked=True)

    obj.axs[0].set_title(f'#Entries: {arr_mass.size}; {title}')
    obj.axs[0].set_ylim(bottom=0, top=maxy)
    obj.axs[0].axvline(x=5080, linestyle=':')
    obj.axs[0].axvline(x=5680, linestyle=':')
    obj.axs[0].axvline(x=5280, label=r'$B^+$', color='gray', linestyle='--')

    obj.axs[1].set_ylim(-5, +5)
    obj.axs[1].axhline(y=-3, color='red')
    obj.axs[1].axhline(y=+3, color='red')
    obj.axs[1].set_label('M$(B^+)$[MeV/${}_{c^2}$]')

    out_dir = f'{Data.out_dir}/{test}'
    os.makedirs(out_dir, exist_ok=True)

    plot_path = f'{out_dir}/{name}.png'
    log.info(f'Saving to: {plot_path}')
    plt.savefig(plot_path)
    plt.close('all')

    _plot_weight(arr_sam, 'sample', '-' )
    _plot_weight(arr_dec, 'decay' , '--')
    _plot_weight(arr_wgt, 'Total' , ':' )

    plt.legend()
    plt.title(title)
    plt.savefig(f'{out_dir}/{name}_wgt.png')
    plt.close('all')

    text_path = plot_path.replace('png', 'txt')
    sut.print_pdf(pdf, txt_path=text_path)
#-----------------------------------------------
def _set_selection(d_cut : dict[str,str]) -> None:
    sel.reset_custom_selection()
    sel.set_custom_selection(d_cut = d_cut)
#-----------------------------------------------
@pytest.mark.parametrize('q2bin', ['jpsi', 'psi2', 'high'])
def test_reso(q2bin : str):
    '''
    Tests PRec building in resonant bins
    '''
    trig   = 'Hlt2RD_BuToKpEE_MVA'

    mass   = {'jpsi' : 'B_const_mass_M',
              'psi2' : 'B_const_mass_psi2S_M',
              'high' : 'B_Mass'}[q2bin]

    label  = {
            'jpsi' :r'$M_{DTF}(K^+e^+e^-)$',
            'psi2' :r'$M_{DTF}(K^+e^+e^-)$',
            'high' :r'$M(K^+e^+e^-)$'}[q2bin]

    maxy   = {'jpsi' : 10_000,
              'psi2' :  2_000,
              'high' :    150}[q2bin]

    q2     = {'jpsi' : r'$J/\psi$',
              'psi2' : r'\psi(2S)',
              'high' : r'High'}[q2bin]

    l_samp = [
            'Bu_JpsiX_ee_eq_JpsiInAcc',
            'Bd_JpsiX_ee_eq_JpsiInAcc',
            'Bs_JpsiX_ee_eq_JpsiInAcc',
            ]

    obs=zfit.Space(label, limits=(4500, 6900))

    test = f'reso/{q2bin}'

    d_wgt= {'dec' : 0, 'sam' : 0}
    obp_4=PRec(samples=l_samp, trig=trig, q2bin=q2bin, d_weight=d_wgt)
    pdf_4=obp_4.get_sum(mass=mass, name='PRec_4', obs=obs)

    title=f'$q^2$: {q2}, uncorrected'
    _plot_pdf(pdf_4, test, 'uncorrected', maxy=maxy, title=title)

    d_wgt= {'dec' : 0, 'sam' : 1}
    obp_3=PRec(samples=l_samp, trig=trig, q2bin=q2bin, d_weight=d_wgt)
    pdf_3=obp_3.get_sum(mass=mass, name='PRec_3', obs=obs)

    title=f'$q^2$: {q2}, sample weights'
    _plot_pdf(pdf_3, test, 'sample_weights', maxy=maxy, title=title)

    d_wgt= {'dec' : 1, 'sam' : 0}
    obp_2=PRec(samples=l_samp, trig=trig, q2bin=q2bin, d_weight=d_wgt)
    pdf_2=obp_2.get_sum(mass=mass, name='PRec_2', obs=obs)

    title=f'$q^2$: {q2}, decay weights'
    _plot_pdf(pdf_2, test, 'decay_weights', maxy=maxy, title=title)

    d_wgt= {'dec' : 1, 'sam' : 1}
    obp_1=PRec(samples=l_samp, trig=trig, q2bin=q2bin, d_weight=d_wgt)
    pdf_1=obp_1.get_sum(mass=mass, name='PRec_1', obs=obs)

    title=f'$q^2$: {q2}, both weights'
    _plot_pdf(pdf_1, test, 'both_weights', maxy=maxy, title=title)
#-----------------------------------------------
def test_fit():
    '''
    Tests that the PDF is fittable
    '''
    q2bin  = 'high'
    trig   = 'Hlt2RD_BuToKpEE_MVA'
    mass   = 'B_Mass'
    label  = r'$M(K^+e^+e^-)$'
    maxy   = 150
    q2     = 'High'
    l_samp = [
            'Bu_JpsiX_ee_eq_JpsiInAcc',
            'Bd_JpsiX_ee_eq_JpsiInAcc',
            'Bs_JpsiX_ee_eq_JpsiInAcc',
            ]

    obs=zfit.Space(label, limits=(4500, 6900))

    test = f'reso/fit/{q2bin}'

    d_wgt= {'dec' : 1, 'sam' : 1}
    obp=PRec(samples=l_samp, trig=trig, q2bin=q2bin, d_weight=d_wgt)
    pdf=obp.get_sum(mass=mass, name='PRec_1', obs=obs)

    title=f'$q^2$: {q2}, both weights'
    _plot_pdf(pdf, test, 'both_weights', maxy=maxy, title=title)

    nev = zfit.Parameter('nev', 0, 0, 10_000)
    pdf.set_yield(nev)

    sam = pdf.create_sampler(n=1000)

    obj = Fitter(pdf, sam)
    res = obj.fit()

    sut.save_fit(
            data   =sam,
            model  =pdf,
            res    =res,
            fit_dir=f'{Data.out_dir}/{test}',
            d_const={})
#-----------------------------------------------
@pytest.mark.parametrize('bdt_cut, name', [
    ('mva.mva_prc > 0.0 && mva.mva_cmb > 0.0', '0p0'),
    ('mva.mva_prc > 0.2 && mva.mva_cmb > 0.2', '0p2'),
    ('mva.mva_prc > 0.3 && mva.mva_cmb > 0.3', '0p2'),
    ('mva.mva_prc > 0.4 && mva.mva_cmb > 0.4', '0p2'),
    ('mva.mva_prc > 0.5 && mva.mva_cmb > 0.5', '0p5'),
    ('mva.mva_prc > 0.8 && mva.mva_cmb > 0.8', '0p8'),
    ('mva.mva_prc > 0.9 && mva.mva_cmb > 0.9', '0p9')])
@pytest.mark.parametrize('q2bin'  , ['jpsi', 'psi2'])
def test_bdt(q2bin : str, bdt_cut : str, name : str):
    '''
    Testing application of BDT cuts
    '''
    _set_selection(d_cut = {'bdt' : bdt_cut})

    obs=zfit.Space('mass', limits=(4500, 6000))
    trig   = 'Hlt2RD_BuToKpEE_MVA'
    mass   = {'jpsi' : 'B_const_mass_M', 'psi2' : 'B_const_mass_psi2S_M'}[q2bin]
    maxy   = {'jpsi' : 20_000          , 'psi2' :                  4_000}[q2bin]
    l_samp = [
            'Bu_JpsiX_ee_eq_JpsiInAcc',
            'Bd_JpsiX_ee_eq_JpsiInAcc',
            'Bs_JpsiX_ee_eq_JpsiInAcc']

    test = f'bdt/{q2bin}'

    d_wgt= {'dec' : 1, 'sam' : 1}
    obp=PRec(samples=l_samp, trig=trig, q2bin=q2bin, d_weight=d_wgt)
    pdf=obp.get_sum(mass=mass, name='PRec_1', obs=obs)

    wp    = name.replace('p', '.')
    title = f'$MVA_{{cmb}} > {wp}$ && $MVA_{{prc}} > {wp}$'

    _plot_pdf(pdf, test, name, maxy=maxy, title=title)
#-----------------------------------------------
@pytest.mark.parametrize('brem_cut, name', [
    ('nbrem == 0', 'z'),
    ('nbrem == 1', 'o'),
    ('nbrem >= 2', 't')])
def test_brem(brem_cut : str, name : str):
    '''
    Testing by brem category
    '''
    _set_selection(d_cut = {'brem' : brem_cut})

    q2bin  = 'jpsi'
    obs=zfit.Space('mass', limits=(4500, 6000))
    trig   = 'Hlt2RD_BuToKpEE_MVA'
    mass   = {'jpsi' : 'B_const_mass_M', 'psi2' : 'B_const_mass_psi2S_M'}[q2bin]
    l_samp = [
            'Bu_JpsiX_ee_eq_JpsiInAcc',
            'Bd_JpsiX_ee_eq_JpsiInAcc',
            'Bs_JpsiX_ee_eq_JpsiInAcc']

    test = f'brem/{q2bin}'

    d_wgt= {'dec' : 1, 'sam' : 1}
    obp=PRec(samples=l_samp, trig=trig, q2bin=q2bin, d_weight=d_wgt)
    pdf=obp.get_sum(mass=mass, name='PRec_1', obs=obs)

    brem  = {'z' : 0, 'o' : 1, 't' : 2}[name]
    title = f'Brem: {brem}'
    _plot_pdf(pdf, test, f'bdt_{name}', maxy=3_000, title=title)
#-----------------------------------------------
def test_cache():
    '''
    Testing caching of PDF
    '''
    q2bin  = 'jpsi'

    obs=zfit.Space('mass', limits=(4500, 6000))
    trig   = 'Hlt2RD_BuToKpEE_MVA'
    mass   = {'jpsi' : 'B_const_mass_M', 'psi2' : 'B_const_mass_psi2S_M'}[q2bin]
    maxy   = {'jpsi' : 20_000          , 'psi2' :                  4_000}[q2bin]
    l_samp = [
            'Bu_JpsiX_ee_eq_JpsiInAcc',
            'Bd_JpsiX_ee_eq_JpsiInAcc',
            'Bs_JpsiX_ee_eq_JpsiInAcc',
            ]

    test = f'cache/{q2bin}'

    d_wgt= {'dec' : 1, 'sam' : 1}
    obp=PRec(samples=l_samp, trig=trig, q2bin=q2bin, d_weight=d_wgt)
    pdf=obp.get_sum(mass=mass, name='PRec_1', obs=obs)

    _plot_pdf(pdf, test, 'cache', maxy=maxy, title='cache test')
#-----------------------------------------------
def test_extended():
    '''
    Testing that PDFs are not extended
    '''
    obs=zfit.Space('mass', limits=(4500, 6000))
    trig   = 'Hlt2RD_BuToKpEE_MVA'
    l_samp = [
            'Bu_JpsiX_ee_eq_JpsiInAcc',
            'Bd_JpsiX_ee_eq_JpsiInAcc',
            'Bs_JpsiX_ee_eq_JpsiInAcc',
            ]

    d_wgt= {'dec' : 1, 'sam' : 1}
    obp=PRec(samples=l_samp, trig=trig, q2bin='jpsi', d_weight=d_wgt)
    pdf=obp.get_sum(mass='B_M', name='PRec_1', obs=obs)

    assert pdf.is_extended is False
#-----------------------------------------------
