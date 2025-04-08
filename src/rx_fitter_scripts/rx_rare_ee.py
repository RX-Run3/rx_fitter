'''
Script used to fit the rare mode
'''
from importlib.resources import files

import yaml
import zfit
from zfit.core.data         import Data       as zdata
from zfit.core.basepdf      import BasePDF    as zpdf
from zfit.core.interfaces   import ZfitSpace  as zobs
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
# --------------------------------------------------------------
def _load_config(component : str) -> dict:
    cfg_path = files('rx_fitter_data').joinpath(f'rare_fit/{Data.version}/rk_ee/{component}.yaml')
    with open(cfg_path, encoding='utf-8') as ifile:
        cfg = yaml.safe_load(ifile)

    return cfg
# --------------------------
def _get_pdf_cmb() -> zpdf:
    pdf= cmp.get_cb(obs=obs, kind=kind, cfg=cfg)
    print_pdf(pdf)
# --------------------------
def _get_pdf() -> zpdf:
    pdf_cmb = _get_pdf_cmb()

    return pdf_cmb
# --------------------------
def _get_data() -> zdata:
    gtr = RDFGetter(sample=Data.sample, trigger=Data.trigger)
    rdf = gtr.get_rdf()
    rdf = sel.apply_full_selection(rdf=rdf, project='RK', trigger=Data.trigger, q2bin=Data.q2bin, process=Data.sample)





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
