'''
Module with tests for the PRec class
'''

import os
import ROOT
import zfit
import mplhep
import pytest
import matplotlib.pyplot as plt

from ROOT                   import EnableImplicitMT
from dmu.stats.utilities    import print_pdf
from dmu.stats.zfit_plotter import ZFitPlotter
from dmu.logging.log_store  import LogStore
from rx_data.rdf_getter     import RDFGetter
from rx_fitter.prec         import PRec as prld

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

    EnableImplicitMT(10)

    plt.style.use(mplhep.style.LHCb2)

    RDFGetter.samples = {
        'main'       : '/home/acampove/external_ssd/Data/samples/main.yaml',
        'mva'        : '/home/acampove/external_ssd/Data/samples/mva.yaml',
        'cascade'    : '/home/acampove/external_ssd/Data/samples/cascade.yaml',
        'jpsi_misid' : '/home/acampove/external_ssd/Data/samples/jpsi_misid.yaml',
        'hop'        : '/home/acampove/external_ssd/Data/samples/hop.yaml',
        }
#-----------------------------------------------
def _plot_pdf(pdf, test : str, name : str, maxy : str):
    arr_mass = pdf.arr_mass
    arr_wgt  = pdf.arr_wgt

    obj = ZFitPlotter(data=arr_mass, model=pdf, weights=arr_wgt)
    obj.plot(stacked=True, ext_text=f'{name}\n#Entries: {arr_mass.size}')
    obj.axs[0].set_ylim(bottom=0, top=maxy)
    obj.axs[0].axvline(x=5080, linestyle=':')
    obj.axs[0].axvline(x=5680, linestyle=':')

    out_dir = f'{Data.out_dir}/{test}'
    os.makedirs(out_dir, exist_ok=True)

    name      = name.replace(' ', '_')
    plot_path = f'{out_dir}/{name}.png'
    log.info(f'Saving to: {plot_path}')
    plt.savefig(plot_path)
    plt.close('all')

    text_path = plot_path.replace('png', 'txt')
    print_pdf(pdf, txt_path=text_path)
#-----------------------------------------------
@pytest.mark.parametrize('q2bin', ['jpsi', 'psi2'])
def test_reso(q2bin : str):
    '''
    Tests PRec building in resonant bins
    '''
    obs=zfit.Space('mass', limits=(4500, 6000))
    trig   = 'Hlt2RD_BuToKpEE_MVA'
    mass   = {'jpsi' : 'B_const_mass_M', 'psi2' : 'B_const_mass_psi2S_M'}[q2bin]
    maxy   = {'jpsi' : 10_000          , 'psi2' :                  2_000}[q2bin]
    bw     = {'jpsi' :  5              , 'psi2' :                     10}[q2bin]
    l_samp = [
            'Bu_JpsiX_ee_eq_JpsiInAcc',
            'Bd_JpsiX_ee_eq_JpsiInAcc',
            'Bs_JpsiX_ee_eq_JpsiInAcc',
            ]

    test = f'reso/{q2bin}'

    d_wgt= {'dec' : 0, 'sam' : 0}
    obp_4=prld(samples=l_samp, trig=trig, q2bin=q2bin, d_weight=d_wgt)
    pdf_4=obp_4.get_sum(mass=mass, name='PRec_4', obs=obs, bandwidth=bw)
    _plot_pdf(pdf_4, test,'Uncorrected', maxy=maxy)

    d_wgt= {'dec' : 0, 'sam' : 1}
    obp_3=prld(samples=l_samp, trig=trig, q2bin=q2bin, d_weight=d_wgt)
    pdf_3=obp_3.get_sum(mass=mass, name='PRec_3', obs=obs, bandwidth=bw)
    _plot_pdf(pdf_3, test,'No decay weights', maxy=maxy)

    d_wgt= {'dec' : 1, 'sam' : 0}
    obp_2=prld(samples=l_samp, trig=trig, q2bin=q2bin, d_weight=d_wgt)
    pdf_2=obp_2.get_sum(mass=mass, name='PRec_2', obs=obs, bandwidth=bw)
    _plot_pdf(pdf_2, test,'No sample weights', maxy=maxy)

    d_wgt= {'dec' : 1, 'sam' : 1}
    obp_1=prld(samples=l_samp, trig=trig, q2bin=q2bin, d_weight=d_wgt)
    pdf_1=obp_1.get_sum(mass=mass, name='PRec_1', obs=obs, bandwidth=bw)
    _plot_pdf(pdf_1, test,'Fully corrected', maxy=maxy)
#-----------------------------------------------
@pytest.mark.parametrize('bdt_cut, name', [
    ('mva.mva_prc > 0.5', '0p5'),
    ('mva.mva_prc > 0.8', '0p8'),
    ('mva.mva_prc > 0.9', '0p9')])
@pytest.mark.parametrize('q2bin'  , ['jpsi', 'psi2'])
def test_bdt(q2bin : str, bdt_cut : str, name : str):
    '''
    Testing application of BDT cuts
    '''
    obs=zfit.Space('mass', limits=(4500, 6000))
    trig   = 'Hlt2RD_BuToKpEE_MVA'
    mass   = {'jpsi' : 'B_const_mass_M', 'psi2' : 'B_const_mass_psi2S_M'}[q2bin]
    maxy   = {'jpsi' : 20_000          , 'psi2' :                  4_000}[q2bin]
    bw     = {'jpsi' :  5              , 'psi2' :                     10}[q2bin]
    l_samp = [
            'Bu_JpsiX_ee_eq_JpsiInAcc',
            'Bd_JpsiX_ee_eq_JpsiInAcc',
            'Bs_JpsiX_ee_eq_JpsiInAcc',
            ]

    test = f'bdt/{q2bin}'

    d_wgt= {'dec' : 1, 'sam' : 1}
    obp=prld(samples=l_samp, trig=trig, q2bin=q2bin, d_weight=d_wgt)
    obp.cuts = {'bdt' : bdt_cut}

    pdf=obp.get_sum(mass=mass, name='PRec_1', obs=obs, bandwidth=bw)
    _plot_pdf(pdf, test, f'bdt_{name}', maxy=maxy)
#-----------------------------------------------
def test_split_type():
    obs=zfit.Space('mass', limits=(4500, 6000))

    d_wgt   = {'dec' : 1, 'sam' : 1}
    for qsq, spl in [('jpsi', 'sam'), ('psi2', 'phy')]:
        obp=prld(samples='b*XcHs', trig='ETOS', q2bin='jpsi', dset='2018', d_weight=d_wgt, split=spl)
        pdf=obp.get_sum(mass='mass_jpsi', name='PRec', obs=obs, bandwidth=10)
        _plot_pdf(pdf, 'bdt_bxXcHs_jpsi_mass_jpsi')
#-----------------------------------------------
