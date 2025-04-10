'''
Module holding SignalScales class
'''
import os
import glob
import math
import json
import numpy
import jacobi
import pandas as pnd

from dmu.generic           import version_management as vman
from dmu.logging.log_store import LogStore

log = LogStore.add_logger('rx_fitter:signal_scales')
# ------------------------------------
class FitParameters:
    '''
    Class used to provide:

    Mass scales
    Mass resolutions
    Brem "resolutions"
    '''
    # -----------------------------------
    def __init__(self):
        '''
        Initializer
        '''
        self._l_brem    = [0, 1, 2]
        self._l_kind    = ['data', 'mc']
        self._fit_dir   = os.environ['FITDIR']

        self._mc_sample = 'Signal'
        self._trigger   = 'Hlt2RD_BuToKpEE_MVA'
        self._mass      = 'B_M_brem_track_2'

        self._df : pnd.DataFrame
        self._is_initialized = False
    #------------------------------------------
    def _name_from_parname(self, name : str) -> str:
        if name.startswith('mu_'):
            return r'$\mu$'

        if name.startswith('sg_'):
            return r'$\sigma$'

        if name.startswith('nSignal'):
            return 'yield'

        if name.startswith('nl_'):
            return 'nl'

        if name.startswith('nr_'):
            return 'nr'

        if name.startswith('al_'):
            return 'al'

        if name.startswith('ar_'):
            return 'ar'

        return None
    #------------------------------------------
    def _df_from_pars(self, d_par : dict[str,list[str]]) -> pnd.DataFrame:
        d_data = {'Name' : [], 'Parameter' : [], 'Value' : [], 'Error' : []}
        for name, [val, err] in d_par.items():
            nickname = self._name_from_parname(name)
            if nickname is None:
                continue

            d_data['Name'     ].append(name)
            d_data['Parameter'].append(nickname)
            d_data['Value'    ].append(val)
            d_data['Error'    ].append(err)

        df = pnd.DataFrame(d_data)

        return df
    # -----------------------------------
    def _get_df_fit(self, kind : str, brem : int) -> pnd.DataFrame:
        sample   = self._mc_sample if kind == 'mc' else 'DATA'

        inp_path = f'{self._fit_dir}/{kind}/jpsi'
        inp_path = vman.get_last_version(dir_path=inp_path, version_only=False)
        inp_wc   = f'{inp_path}/{sample}_{self._trigger}/{self._mass}_{brem}/*/parameters.json'
        l_path   = glob.glob(inp_wc)
        npath    = len(l_path)
        if npath != 1:
            raise ValueError(f'No one and only one path found in {inp_wc}')

        log.debug(f'Looking for parameters in: {inp_path}')
        with open(l_path[0], encoding='utf-8') as ifile:
            d_par = json.load(ifile)

        df = self._df_from_pars(d_par)

        log.debug(f'Parameters for kind {kind}, brem {brem}')
        log.debug(df)

        return df
    # -----------------------------------
    def _frac_from_yield(self, df : pnd.DataFrame) -> pnd.DataFrame:
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
    def _pick_common_parameters(self, df_dt : pnd.DataFrame, df_mc : pnd.DataFrame) -> tuple[pnd.DataFrame, pnd.DataFrame]:
        l_par_dt = df_dt.Parameter.unique().tolist()
        l_par_mc = df_mc.Parameter.unique().tolist()
        l_common = [ par_dt for par_dt in l_par_dt if par_dt in l_par_mc ]

        df_dt = df_dt[df_dt.Parameter.isin(l_common)]
        df_mc = df_mc[df_mc.Parameter.isin(l_common)]

        return df_dt, df_mc
    # -----------------------------------
    def _initialize(self):
        if self._is_initialized:
            return

        l_df_brem = []
        for brem in self._l_brem:
            l_df_kind = []
            for kind in self._l_kind:
                log.debug(f'Extracting parameters for {kind}/{brem}')
                df         = self._get_df_fit(kind = kind, brem = brem)
                df['kind'] = kind

                l_df_kind.append(df)

            [df_dt, df_mc] = l_df_kind
            [df_dt, df_mc] = self._pick_common_parameters(df_dt, df_mc)

            df = pnd.concat([df_dt, df_mc], axis=0)
            df['brem'] = brem

            l_df_brem.append(df)

        df = pnd.concat(l_df_brem, axis=0)
        df = df.reset_index(drop=True)
        df = self._frac_from_yield(df)

        self._is_initialized = True

        self._df = df
    # -----------------------------------
    def get_data(self) -> pnd.DataFrame:
        '''
        Returns dataframe with parameter values from fits to data and MC
        '''
        self._initialize()

        return self._df
    # -----------------------------------
    def _get_parameter_value(self, name : str, is_data : bool) -> tuple[float,float]:
        '''
        Takes name of fitting parameter, returns its value and error
        '''
        self._initialize()

        name = name.replace('scale_', '').replace('reso_', '')

        kind = 'data' if is_data else 'mc'
        df   = self._df.copy()
        df   = df[df.kind == kind ]
        df   = df[df.Name == name ]

        if len(df) != 1:
            log.info(self._df)
            log.info('--->')
            log.error(df)
            raise ValueError(f'Not found one and only one row for: {name}')

        log.debug('Using dataframe:\n')
        log.debug(df)

        val = df.Value.iloc[0]
        err = df.Error.iloc[0]

        return val, err
    # ------------------------------------
    def get_parameter_scale(self, name : str) -> tuple[float,float]:
        '''
        Takes name of scale parameter, returns tuple with value of scale and error
        '''

        if 'Signal' not in name:
            raise ValueError(f'Not a signal parameter: {name}')

        val_dt, err_dt = self._get_parameter_value(name=name, is_data= True)
        val_mc, err_mc = self._get_parameter_value(name=name, is_data=False)

        cov = [[err_dt ** 2,          0],
               [          0, err_mc **2]]

        if   'scale' in name:
            scl, var = jacobi.propagate(lambda x : x[0] - x[1], [val_dt, val_mc], cov)
        elif 'reso'  in name:
            scl, var = jacobi.propagate(lambda x : x[0] / x[1], [val_dt, val_mc], cov)
        else:
            raise ValueError(f'Neither a scale nor a resolution: {name}')

        return scl, math.sqrt(var)
    # ------------------------------------
    def get_brem_fraction(self, name : str, is_data : bool = True) -> tuple[float,float]:
        '''
        Takes name of brem fraction, returns tuple with value of scale and error for data, by default
        is_data flag can override this
        '''
        if name not in ['frac_brem_000', 'frac_brem_001', 'frac_brem_002']:
            raise ValueError(f'Parameter not a brem fraction: {name}')

        [_, _, cat ] = name.split('_')

        val, err = self._get_parameter_value(name=f'nSignal_{cat}', is_data=is_data)

        return val, err
# ------------------------------------
