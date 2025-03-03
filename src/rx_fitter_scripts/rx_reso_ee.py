'''
Script used to fit the resonant mode in the electron channel
'''
import argparse
from importlib.resources import files

import yaml
import ROOT
import zfit

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
    cpath : str

    cfg   : dict
# ------------------------------
def _get_cuts() -> dict[str,str]:
    d_cut         = {}
    d_cut['brem'] = Data.cfg['brem'][Data.nbrem]
    d_sel         = Data.cfg['input']['selection']
    d_cut.update(d_sel)

    return d_cut
# ------------------------------
def _parse_args() -> None:
    parser = argparse.ArgumentParser(description='Script used to fit resonant electron mode')
    parser.add_argument('-b', '--nbrem' , type=int, help='Brem category'   , required=True, choices=[0,1,2])
    parser.add_argument('-m', '--mass'  , type=str, help='Branch with mass', required=True, choices=['B_M', 'B_const_mass_M'])
    args = parser.parse_args()

    Data.nbrem = args.nbrem
    Data.mass  = args.mass
# ------------------------------
def _set_out_dir() -> None:
    q2bin   = Data.cfg['input']['q2bin'  ]
    trigger = Data.cfg['input']['trigger']
    fit_dir = Data.cfg['output']['fit_dir']

    ver_dir = f'{fit_dir}/data/{q2bin}'
    ver_dir = vman.get_last_version(dir_path=ver_dir, version_only=False)
    out_dir = f'{ver_dir}/DATA_{trigger}/{Data.mass}_{Data.nbrem}/full_model'

    Data.cfg['out_dir'] = out_dir
# ------------------------------
def _load_config() -> None:
    cfg_path = files('rx_fitter_data').joinpath('config/v1.yaml')
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

    _set_out_dir()
    obj = DTFitter(rdf = rdf, components = l_cmp, cfg=cfg_fit)
    obj.fit(constraints = d_const)
# ------------------------------
def _get_components() -> list[FitComponent]:
    obs     = zfit.Space(Data.mass, limits=_get_limits())
    q2bin   = Data.cfg['input']['q2bin']
    trigger = Data.cfg['input']['trigger']
    l_fcm   = []

    if Data.cfg['fiting']['components']['Signal']:
        sample = Data.cfg['fitting']['config']['Signal']['sample']
        fcm    = cmp.get_mc(obs = obs, name= sample, trigger=trigger, q2bin=q2bin, nbrem=Data.nbrem)
        l_fcm.append(fcm)

    if Data.cfg['fiting']['components']['Cabibbo']:
        sample = Data.cfg['fitting']['config']['Cabibbo']['sample']
        fcm    = cmp.get_mc(obs = obs, name= sample, trigger=trigger, q2bin=q2bin, nbrem=Data.nbrem)
        l_fcm.append(fcm)

    if Data.cfg['fiting']['components']['combinatorial']:
        kind  = Data.cfg['fitting']['config']['combinatorial']['kind']
        fcm   = cmp.get_cb(obs = obs, kind= kind)
        l_fcm.append(fcm)

    if Data.cfg['fitting']['components']['PRec']:
        d_cut = _get_cuts()
        bw    = Data.cfg['fitting']['config']['PRec']['bw']
        fcm   = cmp.get_prc(obs= obs, name= Data.nbrem, trigger=trigger, q2bin=q2bin, cuts = d_cut, bw = bw)
        l_fcm.append(fcm)

    return l_fcm
# ------------------------------
def main():
    '''
    Start here
    '''
    _parse_args()
    _load_config()

    l_cmp   = _get_components()
    _fit_data(l_cmp)
# ------------------------------
if __name__ == '__main__':
    main()
