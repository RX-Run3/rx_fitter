'''
Script used to fit the rare mode
'''
from importlib.resources import files

import yaml
import zfit
from zfit.core.data         import Data       as zdata
from zfit.core.basepdf      import BasePDF    as zpdf
from zfit.core.interfaces   import ZfitSpace  as zobs
from zfit.core.parameter    import Parameter  as zpar
from dmu.stats.fitter       import Fitter
from dmu.logging.log_store  import LogStore
from rx_data.rdf_getter     import RDFGetter
from rx_selection           import selection  as sel
from rx_fitter              import components as cmp

log=LogStore.add_logger('rx_fitter:rx_rare_ee')
# --------------------------
class Data:
    '''
    Data class
    '''
    sample  : str  = 'DATA*'
    trigger : str  = 'Hlt2RD_BuToKpEE_MVA'
    version : str  = 'v1'
    q2bin   : str  = 'central'
    mass    : str  = 'B_M_brem_track_2'
    minx    : int  = 4500
    maxx    : int  = 6000
    obs     : zobs = zfit.Space(mass, limits=(minx, maxx))
    nsig    : zpar = zfit.Parameter('nsig', 0, 0, 10_000)
# --------------------------------------------------------------
def _load_config(component : str) -> dict:
    cfg_path = files('rx_fitter_data').joinpath(f'rare_fit/{Data.version}/rk_ee/{component}.yaml')
    with open(cfg_path, encoding='utf-8') as ifile:
        cfg = yaml.safe_load(ifile)

    return cfg
# --------------------------
def _get_pdf_cmb() -> zpdf:
    cfg  = _load_config(component = 'combinatorial')
    pdf  = cmp.get_cb(obs=Data.obs, q2bin=Data.q2bin, cfg=cfg)
    ncmb = zfit.Parameter('ncmb', 0, 0, 20_000)
    pdf  = pdf.create_extended(ncmb)

    return pdf
# --------------------------
def _get_pdf_prc(sample : str) -> zpdf:
    cfg                   = _load_config(component='bxhsee')
    cfg['input']['q2bin'] = Data.q2bin

    pdf  = cmp.get_kde(obs=Data.obs, sample=sample, nbrem=None, cfg=cfg)
    scale= zfit.Parameter(f's{sample}', 0, 0, 10)
    nprc = zfit.ComposedParameter(f'n{sample}', lambda x : x['nsig'] * x['scale'], params={'nsig' : Data.nsig, 'scale' : scale})
    pdf  = pdf.create_extended(nprc)

    return pdf
# --------------------------
def _get_pdf_sig() -> zpdf:
    cfg  = _load_config(component='signal')
    pdf  = cmp.get_mc_reparametrized(obs=Data.obs, component_name='Signal', cfg=cfg, nbrem=None)
    pdf  = pdf.create_extended(Data.nsig)

    return pdf
# --------------------------
def _get_pdf() -> zpdf:
    pdf_cmb = _get_pdf_cmb()
    pdf_pr1 = _get_pdf_prc(sample='Bu_Kstee_Kpi0_eq_btosllball05_DPC')
    pdf_pr2 = _get_pdf_prc(sample='Bd_Kstee_eq_btosllball05_DPC')
    pdf_pr3 = _get_pdf_prc(sample='Bs_phiee_eq_Ball_DPC')
    pdf_pr4 = _get_pdf_leak(sample='Bu_JpsiK_ee_eq_DPC')
    pdf_sig = _get_pdf_sig()

    pdf = zfit.pdf.SumPDF([pdf_cmb, pdf_jpsi, pdf_pr1, pdf_pr2, pdf_pr3, pdf_sig])

    return pdf
# --------------------------
def _get_data() -> zdata:
    gtr = RDFGetter(sample=Data.sample, trigger=Data.trigger)
    rdf = gtr.get_rdf()
    rdf = sel.apply_full_selection(rdf=rdf, project='RK', trigger=Data.trigger, q2bin=Data.q2bin, process=Data.sample)

    arr_mass = rdf.AsNumpy([Data.mass])[Data.mass]
    data     = zfit.Data.from_numpy(obs=Data.obs, array=arr_mass)

    return data
# --------------------------
def main():
    '''
    Start here
    '''
    pdf  = _get_pdf()
    data = _get_data()

    obj = Fitter(pdf=pdf, data=data)
    res = obj.fit()
# --------------------------
if __name__ == '__main__':
    main()
