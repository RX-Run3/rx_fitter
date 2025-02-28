'''
Script used to plot and tabulate mass scales and resolutions
'''
import os
import re
import json
import argparse

import numpy
import pandas            as pnd
import matplotlib.pyplot as plt

from dmu.logging.log_store import LogStore
from dmu.generic           import version_management as vman

log = LogStore.add_logger('rx_fitter:reso_scale')
#------------------------------------------
class Data:
    '''
    Data class
    '''
    fit_dir   = os.environ['FITDIR']
    mc_sample = 'Bu_JpsiK_ee_eq_DPC'
    trigger   = 'Hlt2RD_BuToKpEE_MVA'
    mass      = 'B_M'
    sgregex   = r'sg_.*_(\d)_flt'

    l_brem    = [0, 1, 2]
    l_kind    = ['data', 'mc']
#------------------------------------------
def _name_from_parname(name : str) -> str:
    if name.startswith('mu_'):
        return r'$\mu$'

    mtch = re.match(Data.sgregex, name)
    if not mtch:
        raise ValueError(f'Cannot match {name} as a width with: {Data.sgregex}')

    nsg = mtch.group(1)

    return f'$\sigma_{nsg}$'
#------------------------------------------
def _df_from_pars(d_par : dict[str,list[str]]) -> pnd.DataFrame:
    d_data = {'Parameter' : [], 'Value' : [], 'Error' : []}
    for name, [val, err] in d_par.items():
        if not name.startswith('mu_') and not name.startswith('sg_'):
            continue

        name = _name_from_parname(name)

        d_data['Parameter'].append(name)
        d_data['Value'    ].append(val)
        d_data['Error'    ].append(err)

    df = pnd.DataFrame(d_data)

    return df
#------------------------------------------
def _get_df_fit(kind : str, brem : int) -> pnd.DataFrame:
    sample   = Data.mc_sample if kind == 'mc' else 'DATA'

    inp_path = f'{Data.fit_dir}/{kind}/jpsi'
    inp_path = vman.get_last_version(dir_path=inp_path, version_only=False)
    inp_path = f'{inp_path}/{sample}_{Data.trigger}/{Data.mass}_{brem}/fit.json'

    with open(inp_path, encoding='utf-8') as ifile:
        d_par = json.load(ifile)

    df = _df_from_pars(d_par)

    return df
#------------------------------------------
def _get_df() -> pnd.DataFrame:
    l_df_kind = []
    for kind in Data.l_kind:
        l_df_brem = []
        for brem in Data.l_brem:
            df         = _get_df_fit(kind = kind, brem = brem)
            df['brem'] = brem
            l_df_brem.append(df)

        df = pnd.concat(l_df_brem, axis=0)
        df['kind'] = kind

        l_df_kind.append(df)

    df = pnd.concat(l_df_kind, axis=0)

    return df
#------------------------------------------
def _prepare_df(df : pnd.DataFrame, kind : str) -> pnd.DataFrame:
    df = df[df.kind == kind]
    df = df.drop(columns=['kind'])
    df = df.set_index('brem')

    return df
#------------------------------------------
def _subtract_df(df_1, df_2) -> pnd.DataFrame:
    sr_val  = df_1.Value - df_2.Value
    sr_err  = numpy.sqrt(df_1.Error ** 2 - df_2.Error ** 2)
    df      = pnd.DataFrame({'brem' : df_1.index, 'Value' : sr_val, 'Error' : sr_err})

    return df
#------------------------------------------
def _divide_df(df_1, df_2) -> pnd.DataFrame:
    v1 = df_1.Value
    v2 = df_2.Value

    e1 = df_1.Error
    e2 = df_2.Error

    vl = v1 / v2
    er = vl * numpy.sqrt(e1/v1 ** 2 + e2/v2 ** 2)
    df = pnd.DataFrame({'brem' : df_1.index, 'Value' : vl, 'Error' : er})

    return df
#------------------------------------------
def _scale_from_df(df : pnd.DataFrame, parameter : str) -> pnd.DataFrame:
    df_dt = _prepare_df(df, 'data')
    df_mc = _prepare_df(df, 'mc'  )

    if 'mu' in parameter:
        df_sc = _subtract_df(df_dt, df_mc)
    else:
        df_sc = _divide_df(df_dt, df_mc)

    return df_sc
#------------------------------------------
def _path_from_par(parameter : str) -> str:
    if 'mu' in parameter:
        return 'scale.png'

    if 'sigma' in parameter and '1' in parameter:
        return 'resolution_1.png'

    if 'sigma' in parameter and '2' in parameter:
        return 'resolution_2.png'

    raise ValueError(f'Parameter not a scale or resolution: {parameter}')
#------------------------------------------
def _ylim_from_par(parameter: str) -> tuple[float,float]:
    if 'mu' in parameter:
        return -20, 20

    return 0.5, 1.5
#------------------------------------------
def main():
    '''
    Starts here
    '''

    df = _get_df()
    for parameter, df_parameter in df.groupby('Parameter'):
        df_parameter = df_parameter.drop(columns=['Parameter'])
        df_scale = _scale_from_df(df_parameter, parameter)
        df_scale.plot(x='brem', y='Value', yerr='Error')

        plt.title(parameter)
        plt.grid()
        plt.legend([])

        fig_path = _path_from_par(parameter)
        t_ylim   = _ylim_from_par(parameter)

        plt.ylim(t_ylim)

        plt.savefig(fig_path)
#------------------------------------------
if __name__ == '__main__':
    main()
