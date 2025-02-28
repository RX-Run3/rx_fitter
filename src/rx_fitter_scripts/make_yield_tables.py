#!/usr/bin/env python3

import utils_noroot as utnr
import pandas       as pnd
import argparse
import os

log=utnr.getLogger(__name__)
#-----------------------------------------------
class data:
    data_version = 'v10.21p2'

    version = None
    out_dir = None
    year    = None 
    l_year  = ['r1', 'r2p1', '2017', '2018']

    pnd.options.display.float_format = '{:.0f}'.format
#---------------------------
def get_formatters():
    f0 = lambda x : x 
    f1 = lambda x : f'{x:.0f}'
    f2 = lambda x : f'{x:.0f}'
    f3 = lambda x : f'{x:.2f}'

    return [f0, f1, f2, f3] 
#-----------------------------------------------
def get_nsig(proc, trig):
    fit_dir   = os.environ['FITDIR']
    file_path = f'{fit_dir}/{data.version}/data/{data.data_version}/{proc}/{data.year}/pars_{trig}.json'
    utnr.check_file(file_path)

    d_par     = utnr.load_json(file_path)
    [nsig, _] = utnr.get_from_dic(d_par, 'nsig_dt')

    return nsig
#-----------------------------------------------
def get_df(etrig):
    npsi2_mm=get_nsig('psi2', 'MTOS')
    npsi2_ee=get_nsig('psi2',  etrig)

    nctrl_mm=get_nsig('ctrl', 'MTOS')
    nctrl_ee=get_nsig('ctrl',  etrig)

    df=pnd.DataFrame(columns=['Channel', 'Jpsi', 'Psi2S'])

    if etrig == 'ETOS':
        df.loc[0] = ['mTOS', nctrl_mm, npsi2_mm]
        df.loc[1] = [etrig , nctrl_ee, npsi2_ee]
    else:
        df.loc[0] = [etrig , nctrl_ee, npsi2_ee]

    df['Ratio'] = df['Jpsi'] / df['Psi2S']

    return df
#-----------------------------------------------
def get_args():
    parser = argparse.ArgumentParser(description='Used to make tables for data yields')
    parser.add_argument('-v', '--version', type =str, help='Version of fits', required=True) 
    parser.add_argument('-y', '--year'   , nargs='+', help='Year', default=data.l_year) 
    args = parser.parse_args()

    data.version = args.version
    data.l_year  = args.year

    fit_dir      = os.environ['FITDIR']
    data.out_dir = utnr.make_dir_path(f'{fit_dir}/{data.version}')
#-----------------------------------------------
def main():
    for data.year in data.l_year:
        df_all = get_df('ETOS')

        out_path = f'{data.out_dir}/data_yields_{data.year}.tex'
        log.visible(f'Saving to: {out_path}')
        df_all.to_latex(buf=open(out_path, 'w'), formatters=get_formatters(), index=False)
#-----------------------------------------------
if __name__ == '__main__':
    get_args()
    main()
#-----------------------------------------------

