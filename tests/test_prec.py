'''
Module with tests for the PRec class
'''

import mplhep
import pytest
import matplotlib.pyplot as plt

from dmu.stats.fitter       import Fitter
from dmu.stats              import utilities as sut
from dmu.stats.zfit         import zfit
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
def _set_selection(d_cut : dict[str,str]) -> None:
    sel.reset_custom_selection()
    sel.set_custom_selection(d_cut = d_cut)
#-----------------------------------------------
@pytest.mark.parametrize('q2bin', ['central', 'jpsi', 'psi2', 'high'])
def test_reso(q2bin : str):
    '''
    Tests PRec building in resonant bins
    '''
    trig   = 'Hlt2RD_BuToKpEE_MVA'

    mass   = {
            'central' : 'B_Mass',
            'jpsi'    : 'B_const_mass_M',
            'psi2'    : 'B_const_mass_psi2S_M',
            'high'    : 'B_Mass'}[q2bin]

    label  = {
            'central' :r'$M(K^+e^+e^-)$',
            'jpsi'    :r'$M_{DTF}(K^+e^+e^-)$',
            'psi2'    :r'$M_{DTF}(K^+e^+e^-)$',
            'high'    :r'$M(K^+e^+e^-)$'}[q2bin]

    maxy   = {
            'central' : 1000,
            'jpsi'    : 10_000,
            'psi2'    :  2_000,
            'high'    :    150}[q2bin]

    q2     = {
            'central' : 'Central',
            'jpsi'    : r'$J/\psi$',
            'psi2'    : r'\psi(2S)',
            'high'    : r'High'}[q2bin]

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
    PRec.plot_pdf(pdf_4, name='uncorrected', maxy=maxy, title=title, out_dir=f'{Data.out_dir}/{test}')

    d_wgt= {'dec' : 0, 'sam' : 1}
    obp_3=PRec(samples=l_samp, trig=trig, q2bin=q2bin, d_weight=d_wgt)
    pdf_3=obp_3.get_sum(mass=mass, name='PRec_3', obs=obs)

    title=f'$q^2$: {q2}, sample weights'
    PRec.plot_pdf(pdf_3, name='sample_weights', maxy=maxy, title=title, out_dir=f'{Data.out_dir}/{test}')

    d_wgt= {'dec' : 1, 'sam' : 0}
    obp_2=PRec(samples=l_samp, trig=trig, q2bin=q2bin, d_weight=d_wgt)
    pdf_2=obp_2.get_sum(mass=mass, name='PRec_2', obs=obs)

    title=f'$q^2$: {q2}, decay weights'
    PRec.plot_pdf(pdf_2, name='decay_weights', maxy=maxy, title=title, out_dir=f'{Data.out_dir}/{test}')

    d_wgt= {'dec' : 1, 'sam' : 1}
    obp_1=PRec(samples=l_samp, trig=trig, q2bin=q2bin, d_weight=d_wgt)
    pdf_1=obp_1.get_sum(mass=mass, name='PRec_1', obs=obs)

    title=f'$q^2$: {q2}, both weights'
    PRec.plot_pdf(pdf_1, name='both_weights', maxy=maxy, title=title, out_dir=f'{Data.out_dir}/{test}')
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
    PRec.plot_pdf(pdf, name='both_weights', maxy=maxy, title=title, out_dir=f'{Data.out_dir}/{test}')

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

    PRec.plot_pdf(pdf, name, maxy=maxy, title=title, out_dir=f'{Data.out_dir}/{test}')
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
    PRec.plot_pdf(pdf, f'bdt_{name}', maxy=3_000, title=title, out_dir=f'{Data.out_dir}/{test}')
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

    PRec.plot_pdf(pdf, 'cache', maxy=maxy, title='cache test', out_dir=f'{Data.out_dir}/{test}')
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
    pdf=obp.get_sum(mass='B_Mass', name='PRec_1', obs=obs)

    assert pdf.is_extended is False
#-----------------------------------------------
def test_low_stats():
    '''
    Testing with low statistics sample, tight MVA
    '''
    _set_selection(d_cut = {'bdt' : 'mva_cmb > 0.9 && mva_prc > 0.9'})

    obs=zfit.Space('mass', limits=(4500, 6000))
    trig   = 'Hlt2RD_BuToKpEE_MVA'
    l_samp = [
            'Bu_JpsiX_ee_eq_JpsiInAcc',
            'Bd_JpsiX_ee_eq_JpsiInAcc',
            'Bs_JpsiX_ee_eq_JpsiInAcc',
            ]

    d_wgt= {'dec' : 1, 'sam' : 1}
    obp=PRec(samples=l_samp, trig=trig, q2bin='high', d_weight=d_wgt)
    pdf=obp.get_sum(mass='B_Mass', name='PRec_1', obs=obs)

    PRec.plot_pdf(pdf, name='pdf', title='', out_dir=f'{Data.out_dir}/low_stats')
#-----------------------------------------------
