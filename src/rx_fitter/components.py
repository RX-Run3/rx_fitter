'''
Module with functions needed to provide fit components
'''
import copy

from rx_calibration.hltcalibration.fit_component import FitComponent
from rx_fitter.prec                              import PRec

# ------------------------------------
class Data:
    '''
    Data class
    '''
    cfg = {
            'out_dir': 'plots/fit',
            'fitting':
            {
                'error_method'  : 'minuit_hesse',
                'weights_column': 'weights',
                'ntries'        : 20,
                'pvalue'        : 0.02,
                },
            'plotting' :
            {
                'nbins'   : 50,
                'stacked' : True,
                },
            }
# ------------------------------------
def get_mc(name : str, q2bin : str, model : list[str]) -> FitComponent:
    wp_name        = _get_wp_name()
    cfg            = copy.deepcopy(Data.cfg)
    cfg['name']    = name
    out_dir        = cfg['out_dir']
    cfg['out_dir'] = f'{out_dir}/{q2_bin}/{name}/{wp_name}'

    rdf   = _get_rdf(sample)
    rdf   = rdf.Define('weights', '1')

    l_pdf, l_shr = _get_fitting_model(sample)
    if l_pdf == ['kde']:
        pdf = None
    else:
        mod   = ModelFactory(preffix=sample, obs = Data.obs, l_pdf = l_pdf, l_shared=l_shr)
        pdf   = mod.get_pdf()

    obj   = FitComponent(cfg=cfg, rdf=rdf, pdf=pdf, obs=Data.obs)
# ------------------------------------
def get_prc(channel : str, is_dtf : bool) -> FitComponent:
    '''
    Function returning FitComponent object for Partially reconstructed background 

    channel (str): Electron (ee) or muon (mm) channel
    is_dtf  (bool): Signals if the component was obtaining by constraining the dilepton to the J/psi mass
    '''
    
