'''
Script meant to test models to fit MC samples
'''
import argparse
from importlib.resources import files

import ROOT
import yaml
import zfit
from rx_data.rdf_getter import RDFGetter
from rx_fitter          import components as cmp

# --------------------------------
class Data:
    '''
    data class
    '''
    trigger = 'Hlt2RD_BuToKpEE_MVA'
    sample  = 'Bu_JpsiK_ee_eq_DPC'
    q2bin   = 'jpsi'

    cfg     : dict

    obs_name: str
    nbrem   : int
    model   : str
# --------------------------------
def _load_config():
    cfg_path = files('rx_fitter_data').joinpath('model_tester/reso_ee.yaml')
    with open(cfg_path, encoding='utf-8') as ifile:
        Data.cfg = yaml.safe_load(ifile)
# --------------------------------
def _get_obs():
    varname      = Data.obs_name
    [minx, maxx] = Data.cfg['binning'][varname]

    obs=zfit.Space(varname, limits=(minx, maxx))

    return obs
# --------------------------------
def _parse_args():
    parser = argparse.ArgumentParser(description='Script used to test fitting models')
    parser.add_argument('-m', '--model'  , type=str, help='Nickname of model' , required=True)
    parser.add_argument('-b', '--nbrem'  , type=int, help='Bremstrahlung category' , required=True, choices=[-1, 0, 1, 2])
    parser.add_argument('-o', '--obsname', type=str, help='Name of observable' , required=True, choices=['B_M', 'B_const_mass_M'])
    args = parser.parse_args()

    Data.model    = args.model
    Data.nbrem    = args.nbrem
    Data.obs_name = args.obsname
# --------------------------------
def main():
    '''
    Start here
    '''

    _parse_args()
    _load_config()

    obs                 = _get_obs()
    l_mod               = Data.cfg['models'][Data.model]
    RDFGetter.samples   = Data.cfg['samples']
    cmp.Data.cfg.out_dir= Data.cfg['out_dir']

    cmp_sig = cmp.get_mc(obs    = obs,
                         sample = Data.sample,
                         trigger= Data.trigger,
                         q2bin  = Data.q2bin,
                         nbrem  = Data.nbrem,
                         model  = l_mod)
    cmp_sig.run()
# --------------------------------
if __name__ == '__main__':
    main()
