'''
Script meant to test models to fit MC samples
'''
import argparse
from importlib.resources import files

import ROOT
import yaml
import zfit
from dmu.logging.log_store import LogStore
from rx_fitter             import components as cmp

log = LogStore.add_logger('rx_fitter:model_tester')
# --------------------------------
class Data:
    '''
    data class
    '''
    cfg       : dict

    obs_name  : str
    version   : str
    nbrem     : int
    model     : str
    level     : int
    selection : str
# --------------------------------
def _initialize():
    LogStore.set_level('rx_fitter:components', Data.level)
    LogStore.set_level('rx_data:rdf_getter'  , Data.level)
    _load_config()
# --------------------------------
def _load_config():
    cfg_path = files('rx_fitter_data').joinpath(f'model_tester/{Data.version}/reso_ee.yaml')
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
    parser.add_argument('-m', '--model'    , type=str, help='Nickname of model'        , required=True)
    parser.add_argument('-v', '--vers'     , type=str, help='Version of configuration' , required=True)
    parser.add_argument('-b', '--nbrem'    , type=int, help='Bremstrahlung category'   , required=True , choices=[-1, 0, 1, 2])
    parser.add_argument('-o', '--obsname'  , type=str, help='Name of observable'       , required=True , choices=['B_M', 'B_const_mass_M'])
    parser.add_argument('-s', '--selection', type=str, help='Name of selection'        , required=True , choices=['no_prc', 'default'])
    parser.add_argument('-l', '--level'    , type=int, help='Logging level'            , default=20    , choices=[10, 20, 30])
    args = parser.parse_args()

    Data.model    = args.model
    Data.version  = args.vers
    Data.nbrem    = args.nbrem
    Data.obs_name = args.obsname
    Data.level    = args.level
    Data.selection= args.selection
# --------------------------------
def main():
    '''
    Start here
    '''

    _parse_args()
    _initialize()
    component_name = Data.cfg['input']['name']
    l_model        = Data.cfg['models'][Data.model]

    Data.cfg['components'][component_name][Data.nbrem]['model'] = l_model
    d_sel = Data.cfg['input']['selection'][Data.selection]
    log.debug('Overriding selection')
    for name, expr in d_sel.items():
        log.debug(f'{name:<20}{expr}')
    Data.cfg['input']['selection'] = d_sel
    fit_dir = Data.cfg['output']['fit_dir']
    Data.cfg['output']['fit_dir'] = f'{fit_dir}/{Data.selection}'

    obs     = _get_obs()
    cmp_sig = cmp.get_mc(obs = obs, component_name=component_name, nbrem  = Data.nbrem, cfg=Data.cfg)

    cmp_sig.run()
# --------------------------------
if __name__ == '__main__':
    main()
