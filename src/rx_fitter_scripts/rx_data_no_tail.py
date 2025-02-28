'''
Script used to fit the resonant mode in the electron channel
'''
import os
import argparse

import ROOT
import zfit
from ROOT                                        import EnableImplicitMT
from dmu.stats.model_factory                     import ModelFactory
from dmu.logging.log_store                       import LogStore
from rx_calibration.hltcalibration.fit_component import FitComponent
from rx_data.rdf_getter                          import RDFGetter
from rx_fitter                                   import components as cmp

log = LogStore.add_logger('rx_fitter:rx_data_no_tail')
# ------------------------------
class Data:
    '''
    Data class
    '''
    EnableImplicitMT(8)
    dtf_tail = 5160
    out_dir  : str
    l_model  : list[str]
    mass     = 'B_M'
    cfg      = {
            'name'   : 'data_no_tail',
            'fitting':
            {
                'error_method'  : 'minuit_hesse',
                'weights_column': 'weights',
                },
            'plotting' :
            {
                'nbins'   : 50,
                'stacked' : True,
                },
            }

    data_dir = os.environ['DATADIR']
    RDFGetter.samples = {
            'main'       : f'{data_dir}/samples/main.yaml',
            'mva'        : f'{data_dir}/samples/mva.yaml',
            'hop'        : f'{data_dir}/samples/hop.yaml',
            'cascade'    : f'{data_dir}/samples/cascade.yaml',
            'jpsi_misid' : f'{data_dir}/samples/jpsi_misid.yaml'}
# ------------------------------
def _get_cuts() -> dict[str,str]:
    d_cut = {}
    d_cut['nbrem'] = f'nbrem == {Data.nbrem}' if Data.nbrem in [0, 1] else f'nbrem >= {Data.nbrem}'
    d_cut['mass']  = f'B_const_mass_M > {Data.dtf_tail}'

    return d_cut
# ------------------------------
def _parse_args() -> None:
    parser = argparse.ArgumentParser(description='Script used to fit data in electron channel Jpsi bin, after removal of PRec')
    parser.add_argument('-b', '--nbrem' , type=int , help='Brem category', required=True, choices=[0,1,2])
    parser.add_argument('-m', '--model' , nargs='+', help='List of PDFs' , required=True, choices=['suj', 'cbl', 'cbr', 'dscb', 'gauss'])
    parser.add_argument('-o', '--odir'  , type=str , help='Directory where outputs will go', required=True)
    args = parser.parse_args()

    Data.nbrem   = args.nbrem
    Data.l_model = args.model
    Data.out_dir = args.odir
# ------------------------------
def main():
    '''
    Start here
    '''
    _parse_args()

    trigger = 'Hlt2RD_BuToKpEE_MVA'
    q2bin   = 'jpsi'
    d_cut   = _get_cuts()
    obs     = zfit.Space(Data.mass, limits=[4500,6000])
    rdf     = cmp.get_rdf(sample='DATA*', q2bin=q2bin, trigger=trigger, cuts=d_cut)

    mod     = ModelFactory(
                preffix = 'data_no_tail',
                obs     = obs,
                l_pdf   = Data.l_model,
                l_shared= ['mu', 'sg'],
                l_float = ['mu', 'sg'])
    pdf= mod.get_pdf()

    models = '_'.join(Data.l_model)
    Data.cfg['out_dir']= f'{Data.out_dir}/nbrem_{Data.nbrem:03}/{models}'
    obj= FitComponent(cfg=Data.cfg, rdf=rdf, pdf=pdf)
    obj.run()
# ------------------------------
if __name__ == '__main__':
    main()
