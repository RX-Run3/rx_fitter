'''
Module with tests for MCParPdf class
'''
import os
import copy

import ROOT
import zfit
from ROOT                 import RDataFrame
from rx_data.rdf_getter   import RDFGetter
from rx_fitter.mc_par_pdf import MCParPdf
from rx_fitter            import components as cmp

# ------------------------------------------
class Data:
    '''
    data class
    '''
    fit_dir   = os.environ['FITDIR']
    cache_dir = '/tmp/tests/rx_fits/mc_par_pdf'
    cfg       = {
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

    obs = zfit.Space('B_M', limits=[4500, 6000])

    RDFGetter.samples = {
        'main'       : '/home/acampove/external_ssd/Data/samples/main.yaml',
        'mva'        : '/home/acampove/external_ssd/Data/samples/mva.yaml',
        'hop'        : '/home/acampove/external_ssd/Data/samples/hop.yaml',
        'cascade'    : '/home/acampove/external_ssd/Data/samples/cascade.yaml',
        'jpsi_misid' : '/home/acampove/external_ssd/Data/samples/jpsi_misid.yaml'}
# ------------------------------------------
def _get_rdf(cfg : dict) -> RDataFrame:
    nbrem   = cfg['nbrem']
    cuts    = {'nbrem' : f'nbrem == {nbrem}'}
    trigger = cfg['trigger']
    q2bin   = cfg['q2bin']
    sample  = cfg['name']

    return cmp.get_rdf(sample, q2bin, trigger, cuts)
# ------------------------------------------
def test_simple():
    '''
    Simplest test of MCParPdf
    '''
    cfg            = copy.deepcopy(Data.cfg)
    cfg['name'   ] = 'Bu_JpsiK_ee_eq_DPC'
    cfg['q2bin'  ] = 'jpsi'
    cfg['trigger'] = 'Hlt2RD_BuToKpEE_MVA'
    cfg['nbrem'  ] = 1
    cfg['fvers'  ] = None
    cfg['shared' ] = ['mu']
    cfg['model'  ] = ['suj', 'dscb']
    cfg['pfloat' ] = ['mu', 'sg']

    rdf   = _get_rdf(cfg=cfg)
    obj   = MCParPdf(rdf=rdf, obs=Data.obs, cfg=cfg)

    return obj.get_fcomp()
# ------------------------------------------
