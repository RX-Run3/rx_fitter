'''
Script used to fit the resonant mode in the electron channel
'''
import argparse
from importlib.resources import files

import yaml
import ROOT
import zfit

from ROOT                                        import EnableImplicitMT
from dmu.generic                                 import version_management as vman
from dmu.logging.log_store                       import LogStore
from rx_calibration.hltcalibration.dt_fitter     import DTFitter
from rx_calibration.hltcalibration.fit_component import FitComponent
from rx_fitter                                   import components as cmp

log = LogStore.add_logger('rx_fitter:rx_reso_ee')
# ------------------------------
class Data:
    '''
    Data class
    '''
    nbrem : int
    mass  : str
    cfg   : dict
    vers  : str
    level : int
# ------------------------------
def _get_cuts() -> dict[str,str]:
    d_cut         = {}
    d_cut['brem'] = f'nbrem == {Data.nbrem}' if Data.nbrem in [0, 1] else f'nbrem >= {Data.nbrem}'
    d_sel         = Data.cfg['input']['selection']
    d_cut.update(d_sel)

    return d_cut
# ------------------------------
def _parse_args() -> None:
    parser = argparse.ArgumentParser(description='Script used to fit resonant electron mode')
    parser.add_argument('-b', '--nbrem' , type=int, help='Brem category'   , required=True, choices=[0,1,2])
    parser.add_argument('-m', '--mass'  , type=str, help='Branch with mass', required=True, choices=['ecalo_bias_B_M', 'B_M', 'B_const_mass_M'])
    parser.add_argument('-v', '--vers'  , type=str, help='Version of fit configuration', required=True)
    parser.add_argument('-l', '--level' , type=int, help='Logging level', default=20, choices=[10, 20, 30])
    args = parser.parse_args()

    Data.nbrem = args.nbrem
    Data.mass  = args.mass
    Data.vers  = args.vers
    Data.level = args.level
# ------------------------------
def _get_out_dir() -> str:
    q2bin   = Data.cfg['input']['q2bin'  ]
    trigger = Data.cfg['input']['trigger']
    fit_dir = Data.cfg['output']['fit_dir']

    ver_dir = f'{fit_dir}/data/{q2bin}'
    ver_dir = vman.get_last_version(dir_path=ver_dir, version_only=False)
    out_dir = f'{ver_dir}/DATA_{trigger}/{Data.mass}_{Data.nbrem}/full_model'

    return out_dir
# ------------------------------
def _load_config() -> None:
    cfg_path = files('rx_fitter_data').joinpath(f'config/{Data.vers}.yaml')
    with open(cfg_path, encoding='utf-8') as ifile:
        Data.cfg = yaml.safe_load(ifile)
# ------------------------------
def _get_limits() -> tuple[int,int]:
    return Data.cfg['fitting']['range'][Data.mass]
# ------------------------------
def _fit_data(l_cmp : list[FitComponent]) -> None:
    if not Data.cfg['fitting']['components']['data']:
        log.info('Skipping fit to data')
        return

    q2bin   = Data.cfg['input']['q2bin']
    trigger = Data.cfg['input']['trigger']
    d_cut   = _get_cuts()

    rdf     = cmp.get_rdf(sample='DATA*', q2bin=q2bin, trigger=trigger, cuts = d_cut)
    d_const = {
            'nPRec'                : [0, 1],
            'nBu_JpsiPi_ee_eq_DPC' : [0, 1],
            }

    cfg_fit = Data.cfg['fitting']['config']['data']

    cfg_fit['out_dir']  = _get_out_dir()
    obj = DTFitter(rdf  = rdf, components = l_cmp, cfg=cfg_fit)
    obj.fit(constraints = d_const)
# ------------------------------
def _get_components() -> list[FitComponent]:
    obs     = zfit.Space(Data.mass, limits=_get_limits())
    l_fcm   = []

    if Data.cfg['fitting']['components']['combinatorial']:
        kind  = Data.cfg['fitting']['config']['combinatorial']['kind']
        fcm   = cmp.get_cb(obs = obs, kind= kind)
        l_fcm.append(fcm)

    if Data.cfg['fitting']['components']['PRec']:
        fcm   = cmp.get_prc(obs= obs, nbrem=Data.nbrem, cfg=Data.cfg)
        l_fcm.append(fcm)

    if Data.cfg['fitting']['components']['Cabibbo']:
        fcm    = cmp.get_mc(obs=obs, component_name='Cabibbo', nbrem=Data.nbrem, cfg=Data.cfg)
        l_fcm.append(fcm)

    if Data.cfg['fitting']['components']['Signal']:
        fcm    = cmp.get_mc(obs=obs, component_name='Signal', nbrem=Data.nbrem, cfg=Data.cfg)
        l_fcm.append(fcm)

    for fcm in l_fcm:
        fcm.run()

    return l_fcm
# ------------------------------
def _initialize():
    _load_config()
    EnableImplicitMT(10)
    LogStore.set_level('rx_fitter:components'        , Data.level)
    LogStore.set_level('rx_calibration:fit_component', Data.level)
# ------------------------------
def main():
    '''
    Start here
    '''
    _parse_args()
    _initialize()

    l_cmp   = _get_components()
    _fit_data(l_cmp)
# ------------------------------
if __name__ == '__main__':
    main()
