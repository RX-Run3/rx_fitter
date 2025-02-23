'''
Module with tests for the PRec class
'''

import os
import ROOT
import zfit
import mplhep
import pytest
import matplotlib.pyplot as plt

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

    bp_psjp     = '(abs(Jpsi_MC_MOTHER_ID) == 100443) & (abs(Jpsi_MC_GD_MOTHER_ID) == 521) & (abs(H_MC_MOTHER_ID) == 521)'
    bd_psks     = '(abs(Jpsi_MC_MOTHER_ID) ==    511) & (abs(H_MC_MOTHER_ID) == 313) & (abs(H_MC_GD_MOTHER_ID) == 511) & (abs(Jpsi_TRUEID) == 100443)'
    bp_psks     = '(abs(Jpsi_MC_MOTHER_ID) ==    521) & (abs(H_MC_MOTHER_ID) == 323) & (abs(H_MC_GD_MOTHER_ID) == 521) & (abs(Jpsi_TRUEID) == 100443)'

    neg_bp_psjp = bp_psjp.replace('==', '!=').replace('&' , '|')
    neg_bd_psks = bd_psks.replace('==', '!=').replace('&' , '|')
    neg_bp_psks = bp_psks.replace('==', '!=').replace('&' , '|')

    bx_jpkp     = f' (abs(H_TRUEID) == 321) & (abs(Jpsi_TRUEID) == 443)  & ({neg_bp_psjp}) & ({neg_bd_psks}) & ({neg_bp_psks})'
    none        = f'((abs(H_TRUEID) != 321) | (abs(Jpsi_TRUEID) != 443)) & ({neg_bp_psjp}) & ({neg_bd_psks}) & ({neg_bp_psks})'

    d_cut            = {}
    d_cut['bp_psjp'] = bp_psjp
    d_cut['bd_psks'] = bd_psks
    d_cut['bp_psks'] = bp_psks

    d_cut['bx_jpkp'] = bx_jpkp
    d_cut['none'   ] = none
#-----------------------------------------------
@pytest.fixture(scope='session', autouse=True)
def _initialize():
    LogStore.set_level('rx_fitter:prec', 10)
    os.makedirs(Data.out_dir, exist_ok=True)

    plt.style.use(mplhep.style.LHCb2)

    RDFGetter.samples = {
        'main' : '/home/acampove/external_ssd/Data/samples/main.yaml',
        }
#-----------------------------------------------
def _plot_pdf(pdf, name : str, maxy : str):
    arr_mass = pdf.arr_mass
    arr_wgt  = pdf.arr_wgt

    obj = ZFitPlotter(data=arr_mass, model=pdf, weights=arr_wgt)
    obj.plot(stacked=True, ext_text=f'{name}\n#Entries: {arr_mass.size}')
    obj.axs[0].set_ylim(bottom=0, top=maxy)
    obj.axs[0].axvline(x=5080, linestyle=':')
    obj.axs[0].axvline(x=5680, linestyle=':')

    name      = name.replace(' ', '_')
    plot_path = f'{Data.out_dir}/{name}.png'
    log.info(f'Saving to: {plot_path}')
    plt.savefig(plot_path)
    plt.close('all')

    text_path = plot_path.replace('png', 'txt')
    print_pdf(pdf, txt_path=text_path)
#-----------------------------------------------
def test_bdt():
    obs=zfit.Space('mass', limits=(4500, 6000))

    bdt_cut = "(BDT_cmb > 0.977000) & (BDT_prc > 0.480751)"
    #bdt_cut = None

    obp=prld(proc='bdXcHs', trig='ETOS', q2bin='psi2', dset='2018')
    pdf=obp.get_pdf(mass='mass_psi2', cut=bdt_cut, name='PRec', obs=obs, use_weights=True, bandwidth=20)
    _plot_pdf(pdf, 'bdt_bpXcHs_psi2_mass_psi2')

    obp=prld(proc='bpXcHs', trig='ETOS', q2bin='psi2', dset='2018')
    pdf=obp.get_pdf(mass='mass_psi2', cut=bdt_cut, name='PRec', obs=obs, use_weights=True, bandwidth=20)
    _plot_pdf(pdf, 'bdt_bdXcHs_psi2_mass_psi2')
#-----------------------------------------------
def test_bd():
    obs=zfit.Space('mass', limits=(4000, 6000))

    obj=prld(proc='bdXcHs', trig='ETOS', q2bin='jpsi', dset='2018')
    pdf_0=obj.get_pdf(mass='mass', cut='0.0 < BDT_cmb < 0.2', obs=obs, bandwidth=20)
    pdf_1=obj.get_pdf(mass='mass', cut='0.2 < BDT_cmb < 0.4', obs=obs, bandwidth=20)
    pdf_2=obj.get_pdf(mass='mass', cut='0.4 < BDT_cmb < 0.6', obs=obs, bandwidth=20)
    pdf_3=obj.get_pdf(mass='mass', cut='0.6 < BDT_cmb < 0.8', obs=obs, bandwidth=20)
    pdf_4=obj.get_pdf(mass='mass', cut='0.8 < BDT_cmb < 1.0', obs=obs, bandwidth=20)
    pdf_5=obj.get_pdf(mass='mass', cut='1.0 < BDT_cmb < 2.0', obs=obs, bandwidth=20)

    _plot_pdf(pdf_0, 'bd_bdt_00')
    _plot_pdf(pdf_1, 'bd_bdt_01')
    _plot_pdf(pdf_2, 'bd_bdt_02')
    _plot_pdf(pdf_3, 'bd_bdt_03')
    _plot_pdf(pdf_4, 'bd_bdt_04')
    _plot_pdf(pdf_5, 'bd_bdt_05')
#-----------------------------------------------
def test_bp():
    obs=zfit.Space('mass', limits=(4000, 6000))

    obj=prld(proc='bpXcHs', trig='ETOS', q2bin='jpsi', dset='2018')
    pdf_0=obj.get_pdf(mass='mass', cut='0.0 < BDT_cmb < 0.2', obs=obs, bandwidth=20)
    pdf_1=obj.get_pdf(mass='mass', cut='0.2 < BDT_cmb < 0.4', obs=obs, bandwidth=20)
    pdf_2=obj.get_pdf(mass='mass', cut='0.4 < BDT_cmb < 0.6', obs=obs, bandwidth=20)
    pdf_3=obj.get_pdf(mass='mass', cut='0.6 < BDT_cmb < 0.8', obs=obs, bandwidth=20)
    pdf_4=obj.get_pdf(mass='mass', cut='0.8 < BDT_cmb < 1.0', obs=obs, bandwidth=20)
    pdf_5=obj.get_pdf(mass='mass', cut='1.0 < BDT_cmb < 2.0', obs=obs, bandwidth=20)

    _plot_pdf(pdf_0, 'bp_bdt_00')
    _plot_pdf(pdf_1, 'bp_bdt_01')
    _plot_pdf(pdf_2, 'bp_bdt_02')
    _plot_pdf(pdf_3, 'bp_bdt_03')
    _plot_pdf(pdf_4, 'bp_bdt_04')
    _plot_pdf(pdf_5, 'bp_bdt_05')
#-----------------------------------------------
def test_wt(q2bin='psi2', proc='bpXcHs', mass='mass'):
    obs=zfit.Space('mass', limits=(5100, 5680))

    obp=prld(proc=proc, trig='ETOS', q2bin=q2bin, dset='2018')

    bandwidth=15
    if   q2bin == 'jpsi' and mass == 'mass_psi2':
        padding = {'lowermirror' : 1}
    elif q2bin == 'psi2' and mass == 'mass_psi2':
        padding = {'uppermirror' : 1}
    elif q2bin == 'jpsi' and mass == 'mass_jpsi':
        padding = {'uppermirror' : 0.8, 'lowermirror' : 0.0}
        bandwith= 10
    else:
        padding = {}

    pdf_u=obp.get_pdf(mass=mass, use_weights=False, obs=obs, bandwidth=bandwidth, padding=padding)
    pdf_w=obp.get_pdf(mass=mass, use_weights=True , obs=obs, bandwidth=bandwidth, padding=padding)

    proc = proc.replace('*', 'x')

    _plot_pdf(pdf_u, f'{proc}_{q2bin}_{mass}_nwgt')
    _plot_pdf(pdf_w, f'{proc}_{q2bin}_{mass}_ywgt')
#-----------------------------------------------
def do_test_match(q2bin, proc, mass, weights, match):
    obs=zfit.Space('mass', limits=(5100, 5680))

    obp=prld(proc=proc, trig='ETOS', q2bin=q2bin, dset='2018')
    pdf=obp.get_pdf(mass=mass, obs=obs, use_weights=weights, padding=0.1, cut=Data.d_cut[match])
    _plot_pdf(pdf, f'{match}_{proc}_{q2bin}_{mass}', maxy=300)
#-----------------------------------------------
def test_match(weight):
    do_test_match('psi2', 'bpXcHs', 'mass_psi2', weight, 'bp_psjp')
    do_test_match('psi2', 'bdXcHs', 'mass_psi2', weight, 'bd_psks')
    do_test_match('psi2', 'bpXcHs', 'mass_psi2', weight, 'bp_psks')

    do_test_match('psi2', 'bpXcHs', 'mass_psi2', weight, 'bx_jpkp')
    do_test_match('psi2', 'bdXcHs', 'mass_psi2', weight, 'bx_jpkp')
    do_test_match('psi2', 'bpXcHs', 'mass_psi2', weight, 'none')
    do_test_match('psi2', 'bdXcHs', 'mass_psi2', weight, 'none')
#-----------------------------------------------
def test_split():
    obs=zfit.Space('mass', limits=(5100, 5680))

    obp  =prld(proc='bpXcHs', trig='ETOS', q2bin='psi2', dset='2018')
    d_pdf=obp.get_components(mass='mass_psi2', obs=obs, use_weights=True, padding=0.1)

    for name, pdf in d_pdf.items():
        _plot_pdf(pdf, f'{name}_bxXcHs_psi2_mass_psi2', maxy=300)
#-----------------------------------------------
def test_sum():
    obs=zfit.Space('mass', limits=(4250, 5680))

    for ccbar in ['jpsi', 'psi2']:
        maxy = 2000 if ccbar == 'psi2' else 8000
        for dec in [0, 1]:
            for sam in [0, 1]:
                obp=prld(proc='b*XcHs', trig='ETOS', q2bin=ccbar, dset='2018', d_weight={'dec': dec, 'sam' : sam})
                pdf=obp.get_sum(mass=f'mass_{ccbar}', name='PRec', obs=obs, padding=0.1, bandwidth=10)

                _plot_pdf(pdf, f'sum_bxXcHs_{ccbar}_mass_{ccbar}_dec{dec}_sam{sam}', maxy=maxy)
#-----------------------------------------------
def test_sum_allon():
    obs=zfit.Space('mass', limits=(4250, 5680))

    for ccbar in ['jpsi', 'psi2']:
        maxy = 3000 if ccbar == 'psi2' else 9000
        obp=prld(proc='b*XcHs', trig='ETOS', q2bin=ccbar, dset='2018', d_weight={'dec': 1, 'sam' : 1})
        pdf=obp.get_sum(mass=f'mass_{ccbar}', name='PRec', obs=obs, padding=0.1, bandwidth=10)

        _plot_pdf(pdf, f'sum_bxXcHs_{ccbar}_mass_{ccbar}', maxy=maxy)
#-----------------------------------------------
def test_brem():
    obs=zfit.Space('mass', limits=(4000, 5680))

    for nbrem in [0, 1, 2]:
        obp=prld(proc='b*XcHs', trig='ETOS', q2bin='jpsi', dset='2018')
        obp.val_dir = f'{Data.out_dir}/validation'
        obp.nbrem   = nbrem
        pdf=obp.get_sum(mass='mass_jpsi', name='PRec', obs=obs, use_weights=True, padding=0.1)

        _plot_pdf(pdf, f'sum_bxXcHs_psi2_mass_psi2_{nbrem}', maxy=3000)
#-----------------------------------------------
def test_dset():
    obs=zfit.Space('mass', limits=(4000, 5680))

    for dset in ['r1', 'r2p1', '2017', '2018']:
        obp=prld(proc='b*XcHs', trig='ETOS', q2bin='psi2', dset=dset)
        pdf=obp.get_sum(mass='mass_psi2', name='PRec', obs=obs, use_weights=True, padding=0.1)
        _plot_pdf(pdf, f'sum_bxXcHs_psi2_mass_psi2_{dset}', maxy=3000)
#-----------------------------------------------
def test_simple_jpsi():
    obs=zfit.Space('mass', limits=(4000, 6000))

    d_wgt= {'dec' : 1, 'sam' : 1}
    obp_1=prld(proc='b*', trig='ETOS', q2bin='jpsi', dset='2018', d_weight=d_wgt)
    obp_1.val_dir = Data.out_dir
    pdf_1=obp_1.get_sum(mass='mass_jpsi', name='PRec_1', obs=obs, bandwidth=10)
    _plot_pdf(pdf_1, 'Fully corrected', maxy=2500)

    d_wgt= {'dec' : 1, 'sam' : 0}
    obp_2=prld(proc='b*', trig='ETOS', q2bin='jpsi', dset='2018', d_weight=d_wgt)
    obp_2.val_dir = Data.out_dir
    pdf_2=obp_2.get_sum(mass='mass_jpsi', name='PRec_2', obs=obs, bandwidth=10)
    _plot_pdf(pdf_2, 'No sample weights', maxy=3000)

    d_wgt= {'dec' : 0, 'sam' : 1}
    obp_3=prld(proc='b*', trig='ETOS', q2bin='jpsi', dset='2018', d_weight=d_wgt)
    obp_3.val_dir = Data.out_dir
    pdf_3=obp_3.get_sum(mass='mass_jpsi', name='PRec_3', obs=obs, bandwidth=10)
    _plot_pdf(pdf_3, 'No decay weights', maxy=3000)

    d_wgt= {'dec' : 0, 'sam' : 0}
    obp_4=prld(proc='b*', trig='ETOS', q2bin='jpsi', dset='2018', d_weight=d_wgt)
    obp_4.val_dir = Data.out_dir
    pdf_4=obp_4.get_sum(mass='mass_jpsi', name='PRec_4', obs=obs, bandwidth=10)
    _plot_pdf(pdf_4, 'Uncorrected', maxy=3000)
#-----------------------------------------------
def test_simple_psi2():
    obs=zfit.Space('mass', limits=(4500, 6000))

    d_wgt= {'dec' : 1, 'sam' : 1}
    obp_1=prld(proc='b*', trig='ETOS', q2bin='psi2', dset='2018', d_weight=d_wgt)
    obp_1.val_dir = Data.out_dir
    pdf_1=obp_1.get_sum(mass='mass_psi2', name='PRec_1', obs=obs, bandwidth=10)
    _plot_pdf(pdf_1, 'Fully corrected', maxy=400)

    d_wgt= {'dec' : 1, 'sam' : 0}
    obp_2=prld(proc='b*', trig='ETOS', q2bin='psi2', dset='2018', d_weight=d_wgt)
    obp_2.val_dir = Data.out_dir
    pdf_2=obp_2.get_sum(mass='mass_psi2', name='PRec_2', obs=obs, bandwidth=10)
    _plot_pdf(pdf_2, 'No sample weights', maxy=400)

    d_wgt= {'dec' : 0, 'sam' : 1}
    obp_3=prld(proc='b*', trig='ETOS', q2bin='psi2', dset='2018', d_weight=d_wgt)
    obp_3.val_dir = Data.out_dir
    pdf_3=obp_3.get_sum(mass='mass_psi2', name='PRec_3', obs=obs, bandwidth=10)
    _plot_pdf(pdf_3, 'No decay weights', maxy=400)

    d_wgt= {'dec' : 0, 'sam' : 0}
    obp_4=prld(proc='b*', trig='ETOS', q2bin='psi2', dset='2018', d_weight=d_wgt)
    obp_4.val_dir = Data.out_dir
    pdf_4=obp_4.get_sum(mass='mass_psi2', name='PRec_4', obs=obs, bandwidth=10)
    _plot_pdf(pdf_4, 'Uncorrected', maxy=400)
#-----------------------------------------------
def test_split_type():
    obs=zfit.Space('mass', limits=(4500, 6000))

    d_wgt   = {'dec' : 1, 'sam' : 1}
    for qsq, spl in [('jpsi', 'sam'), ('psi2', 'phy')]:
        obp=prld(proc='b*XcHs', trig='ETOS', q2bin='jpsi', dset='2018', d_weight=d_wgt, split=spl)
        pdf=obp.get_sum(mass='mass_jpsi', name='PRec', obs=obs, bandwidth=10)
        _plot_pdf(pdf, 'bdt_bxXcHs_jpsi_mass_jpsi')
#-----------------------------------------------
