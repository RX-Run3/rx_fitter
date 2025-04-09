'''
Script used to plot and tabulate mass scales and resolutions
'''
import os
import json
import glob
import argparse

import numpy
import mplhep
import jacobi
import pandas            as pnd
import matplotlib.pyplot as plt
import dmu.pdataframe.utilities as put

from dmu.logging.log_store import LogStore
from dmu.generic           import version_management as vman

log = LogStore.add_logger('rx_fitter:reso_scale')
#------------------------------------------
class Data:
    '''
    Data class
    '''
    fit_dir   = os.environ['FITDIR']
    mc_sample = 'Signal'
    trigger   = 'Hlt2RD_BuToKpEE_MVA'
    mass      = 'B_M_brem_track_2'
    sgregex   = 'sg_Signal_flt'

    l_brem    = [0, 1, 2]
    l_kind    = ['data', 'mc']
#------------------------------------------
def _name_from_parname(name : str) -> str:
    if name.startswith('mu_'):
        return r'$\mu$'

    if name.startswith('sg_'):
        return r'$\sigma$'

    if name.startswith('nSignal'):
        return 'yield'

    raise ValueError(f'Invalid parameter name {name}')
#------------------------------------------
def _is_par_needed(name : str) -> str:
    if name.startswith('nSignal'):
        return True

    if name.startswith('mu_'):
        return True

    if name.startswith('sg_'):
        return True

    return False
#------------------------------------------
def _df_from_pars(d_par : dict[str,list[str]]) -> pnd.DataFrame:
    d_data = {'Parameter' : [], 'Value' : [], 'Error' : []}
    for name, [val, err] in d_par.items():
        if not _is_par_needed(name):
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
    inp_wc   = f'{inp_path}/{sample}_{Data.trigger}/{Data.mass}_{brem}/*/fit.json'
    l_path   = glob.glob(inp_wc)
    npath    = len(l_path)
    if npath != 1:
        raise ValueError(f'No one and only one path found in {inp_wc}')

    log.debug(f'Looking for parameters in: {inp_path}')
    with open(l_path[0], encoding='utf-8') as ifile:
        d_par = json.load(ifile)

    df = _df_from_pars(d_par)

    return df
#------------------------------------------
def _get_df() -> pnd.DataFrame:
    l_df_kind = []
    for kind in Data.l_kind:
        l_df_brem = []
        for brem in Data.l_brem:
            log.debug(f'Extracting parameters for {kind}/{brem}')
            df         = _get_df_fit(kind = kind, brem = brem)
            df['brem'] = brem
            l_df_brem.append(df)

        df = pnd.concat(l_df_brem, axis=0)
        df['kind'] = kind

        l_df_kind.append(df)

    df = pnd.concat(l_df_kind, axis=0)
    df = df.reset_index(drop=True)
    df = _frac_from_yield(df)

    return df
#------------------------------------------
def _frac_from_yield(df : pnd.DataFrame) -> pnd.DataFrame:
    l_df = []
    for _, df_kind in df.groupby('kind'):
        df          = df_kind[df_kind.Parameter == 'yield']
        frac, cov   = jacobi.propagate(lambda x : x / numpy.sum(x), df.Value.values, df.Error.values)
        df['Value'] = frac
        df['Error'] = numpy.sqrt(numpy.diag(cov))
        df.Parameter= df.Parameter.replace({'yield' : 'frac'})

        for index in df.index:
            df_kind.loc[index] = df.loc[index]

        l_df.append(df_kind)

    df = pnd.concat(l_df, axis=0)

    return df
#------------------------------------------
def _parse_args():
    parser = argparse.ArgumentParser(description='Script used to calculate mass scales and resolution from fit results')
    parser.add_argument('-l', '--log_level' , type=int, help='Logging level', choices=[10, 20, 30])
    args = parser.parse_args()

    _set_log_level(args.log_level)
#------------------------------------------
def _set_log_level(level : int) -> None:
    LogStore.set_level('rx_fitter:reso_scale', level)
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
        df_sc =   _divide_df(df_dt, df_mc)

    return df_sc
#------------------------------------------
def _path_from_par(parameter : str) -> str:
    if 'mu'    in parameter:
        return 'scale.png'

    if 'sigma' in parameter:
        return 'resolution.png'

    if 'frac'  in parameter:
        return 'brem_fraction.png'

    raise ValueError(f'Parameter not a scale or resolution: {parameter}')
#------------------------------------------
def _ylim_from_par(parameter : str) -> tuple[float,float]:
    if 'mu'    in parameter:
        return -20, 0

    if 'sigma' in parameter:
        return 1.0, 2.0

    if 'frac'  in parameter:
        return 0.7, 1.3

    raise NotImplementedError(f'Invalid parameter {parameter}')
#------------------------------------------
def _ylabel_from_par(parameter : str) -> str:
    if 'mu'    in parameter:
        return r'$\mu_{Data}-\mu_{MC}$[MeV]'

    if 'sigma' in parameter:
        return r'$\sigma_{Data}/\sigma_{MC}$'

    if 'frac'  in parameter:
        return r'$frac_{Data}/frac_{MC}$'

    raise NotImplementedError(f'Invalid parameter {parameter}')
#------------------------------------------
def _caption_from_par(parameter : str) -> str:
    if 'mu'    in parameter:
        return r'$\mu$'

    if 'sigma' in parameter:
        return r'$\sigma$'

    if 'frac'  in parameter:
        return r'$f_{brem}$'

    raise NotImplementedError(f'Invalid parameter {parameter}')
#------------------------------------------
def _format_float(val : float) -> str:
    if val < 10:
        return f'{val:.2f}'

    if val < 100:
        return f'{val:.1f}'

    return f'{val:.0f}'
#------------------------------------------
def _value_from_df(row : pnd.Series) -> str:
    val = row.Value
    err = row.Error

    val = _format_float(val)
    err = _format_float(err)

    return f'${val} \pm {err}$'
#------------------------------------------
def _split_kind(df : pnd.DataFrame, name : str) -> pnd.DataFrame:
    df_dt = df[df.Sample == 'data']
    df_mc = df[df.Sample ==   'mc']

    df_dt = df_dt.reset_index(drop=True)
    df_mc = df_mc.reset_index(drop=True)

    df_dt = df_dt.drop(columns=['Sample'])
    df_mc = df_mc.drop(columns=['Sample', 'Brem'])

    df_mc = df_mc.rename(columns = {name : 'MC'})
    df_dt = df_dt.rename(columns = {name : 'Data'})

    df    = pnd.concat([df_dt, df_mc], axis=1)

    return df
#------------------------------------------
def _tabulate(df : pnd.DataFrame, name : str) -> None:
    # Brem fractions will be showin in %
    if 'frac' in name:
        df['Value'] = 100 * df.Value
        df['Error'] = 100 * df.Error

    sr_val = df.apply(_value_from_df, axis=1)
    d_data = {'Sample' : df.kind, 'Brem' : df.brem, name : sr_val}
    df     = pnd.DataFrame(d_data)
    df     = _split_kind(df, name)


    fname  = name.replace('$', '').replace('\\', '').replace('{', '').replace('}', '')

    label  = _caption_from_par(name)
    put.df_to_tex(df, f'./{fname}.tex', caption=label, column_format='ll|l')
#------------------------------------------
def _initialize() -> None:
    plt.style.use(mplhep.style.LHCb2)
#------------------------------------------
def main():
    '''
    Starts here
    '''
    _parse_args()
    _initialize()

    df = _get_df()
    for parameter, df_parameter in df.groupby('Parameter'):
        df_parameter = df_parameter.drop(columns=['Parameter'])

        _tabulate(df_parameter, parameter)

        df_scale = _scale_from_df(df_parameter, parameter)

        df_scale.plot(x='brem', y='Value', yerr='Error')

        plt.grid()
        plt.legend([])

        fig_path = _path_from_par(parameter)
        t_ylim   = _ylim_from_par(parameter)
        ylabel   = _ylabel_from_par(parameter)

        plt.ylabel(ylabel)
        plt.ylim(t_ylim)
        log.info(f'Saving to: {fig_path}')
        plt.savefig(fig_path)
#------------------------------------------
if __name__ == '__main__':
    main()
