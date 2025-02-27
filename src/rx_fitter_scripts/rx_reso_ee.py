'''
Script used to fit the resonant mode in the electron channel
'''
import argparse

import ROOT
import zfit
from dmu.logging.log_store                   import LogStore
from rx_calibration.hltcalibration.dt_fitter import DTFitter
from rx_data.rdf_getter                      import RDFGetter
from rx_fitter                               import components as cmp

log = LogStore.add_logger('rx_fitter:rx_reso_ee')
# ------------------------------
class Data:
    '''
    Data class
    '''
    nbrem : int
    mass  : str
    rng   : str

    cfg   = {
            'error_method' : 'minuit_hesse',
            'out_dir'      : 'plots/fit/data',
            'plotting'     :
            {
                'nbins'   : 50,
                'stacked' : True,
                'd_leg'   : {
                    'Bu_JpsiK_ee_eq_DPC' : r'$B^+\to K^+J/\psi(\to e^+e^-)$',
                    'Bu_JpsiPi_ee_eq_DPC': r'$B^+\to \pi^+J/\psi(\to e^+e^-)$',
                    'combinatorial'      : 'Combinatorial',
                    }
                },
            }

    RDFGetter.samples = {
        'main'       : '/home/acampove/external_ssd/Data/samples/main.yaml',
        'mva'        : '/home/acampove/external_ssd/Data/samples/mva.yaml',
        'hop'        : '/home/acampove/external_ssd/Data/samples/hop.yaml',
        'cascade'    : '/home/acampove/external_ssd/Data/samples/cascade.yaml',
        'jpsi_misid' : '/home/acampove/external_ssd/Data/samples/jpsi_misid.yaml'}
# ------------------------------
def _get_limits() -> tuple[int,int]:
    if Data.rng == 'bcn':
        return 5050, 5600

    if Data.rng == 'wide':
        return 4500, 5600

    if Data.rng == 'largest':
        return 4500, 6000

    raise ValueError(f'Invalid range: {Data.rng}')
# ------------------------------
def _parse_args() -> None:
    parser = argparse.ArgumentParser(description='Script used to fit resonant electron mode')
    parser.add_argument('-b', '--nbrem' , type=int, help='Brem category'   , required=True, choices=[0,1,2])
    parser.add_argument('-m', '--mass'  , type=str, help='Branch with mass', required=True, choices=['B_M', 'B_const_mass_M'])
    parser.add_argument('-r', '--range' , type=str, help='Name of range'   , required=True, choices=['bcn', 'wide', 'largest'])
    args = parser.parse_args()

    Data.nbrem = args.nbrem
    Data.mass  = args.mass
    Data.rng   = args.range
# ------------------------------
def main():
    '''
    Start here
    '''
    _parse_args()

    trigger = 'Hlt2RD_BuToKpEE_MVA'
    q2bin   = 'jpsi'

    t_lim   = _get_limits()
    obs     = zfit.Space(Data.mass, limits=t_lim)

    cmp_sig = cmp.get_mc(obs = obs, sample = 'Bu_JpsiK_ee_eq_DPC' , trigger=trigger, q2bin=q2bin, nbrem=Data.nbrem)
    cmp_csp = cmp.get_mc(obs = obs, sample = 'Bu_JpsiPi_ee_eq_DPC', trigger=trigger, q2bin=q2bin, nbrem=Data.nbrem)
    cmp_prc = cmp.get_prc(obs= obs, trigger=trigger, q2bin=q2bin, nbrem=Data.nbrem)
    cmp_cmb = cmp.get_cb(obs = obs, kind='exp')

    rdf = cmp.get_rdf(sample='DATA*', q2bin=q2bin, trigger=trigger, nbrem=Data.nbrem)

    out_dir = Data.cfg['out_dir']
    Data.cfg['out_dir']= f'{out_dir}/nbrem_{Data.nbrem:03}'

    obj = DTFitter(rdf = rdf, components = [cmp_cmb, cmp_prc, cmp_csp, cmp_sig], cfg=Data.cfg)
    obj.fit()
# ------------------------------
if __name__ == '__main__':
    main()
