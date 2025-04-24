'''
Script used to make tex files from txt files corresponding to zfit PDFs
'''
import os
import glob
import argparse
from importlib.resources import files

import yaml
from dmu.stats.utilities   import pdf_to_tex
from dmu.logging.log_store import LogStore

log = LogStore.add_logger('rx_fitter:tabulate_pdfs')
# --------------------------------
class Data:
    '''
    Data class
    '''
    fit_dir     : str
    d_all_names : dict[dict,dict[str,str]]
# --------------------------------
def _parse_args():
    parser = argparse.ArgumentParser(description='Used to perform several operations on TCKs')
    parser.add_argument('-d', '--fit_dir', type=str, help='Directory where text files are stored', required=True)
    args = parser.parse_args()

    Data.fit_dir = args.fit_dir
# --------------------------------
def _kind_from_path(path : str) -> str:
    fname    = os.path.basename(path)
    [name, _]= fname.split('.')
    [_, kind]= name.split('_')

    return kind
# --------------------------------
def _initialize():
    path = files('rx_fitter_data').joinpath('names/models.yaml')
    with open(path, encoding='utf-8') as ifile:
        Data.d_all_names = yaml.safe_load(ifile)
# --------------------------------
def main():
    '''
    Start here
    '''
    _parse_args()
    _initialize()

    l_path = glob.glob(f'{Data.fit_dir}/*.txt')

    for path in l_path:
        kind = _kind_from_path(path)
        d_names = Data.d_all_names[kind]

        pdf_to_tex(path=path, d_par=d_names)
# --------------------------------
if __name__ == '__main__':
    main()
