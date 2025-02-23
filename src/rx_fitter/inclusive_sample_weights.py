'''
Module with Reader class used to read weights to normalize between inclusive samples
'''

import os

import pandas    as pnd
from dmu.logging.log_store import LogStore
from rx_fitter             import pdg_utils as pu

log = LogStore.add_logger('rx_fitter:inclusive_sample_weights')
#---------------------------
class Reader:
    '''
    Class used to add weights that normalize inclusive samples
    '''
    #---------------------------
    def __init__(self, df : pnd.DataFrame):
        self._df      = df
        self._fu      = 0.408
        self._fs      = 0.100

        self._bu_proc = 'Bu_JpsiX_ee_eq_JpsiInAcc'
        self._bd_proc = 'Bd_JpsiX_ee_eq_JpsiInAcc'
        self._bs_proc = 'Bs_JpsiX_ee_eq_JpsiInAcc'
    #---------------------------
    def _get_br_wgt(self, proc : str) -> float:
        '''
        Will return ratio: 

        decay file br / pdg_br 
        '''
        #--------------------------------------------
        #Decay B+sig
        #0.1596  MyJ/psi    K+           SVS ;
        #--------------------------------------------
        #Decay B0sig
        #0.1920  MyJ/psi    MyK*0        SVV_HELAMP PKHplus PKphHplus PKHzero PKphHzero PKHminus PKphHminus ;
        #--------------------------------------------
        #Decay B_s0sig
        #0.1077  MyJ/psi    Myphi        PVV_CPLH 0.02 1 Hp pHp Hz pHz Hm pHm;
        #--------------------------------------------
        if proc == self._bu_proc:
            return pu.get_bf('B+ --> J/psi(1S) K+') / 0.1596

        if proc == self._bd_proc:
            return pu.get_bf('B0 --> J/psi(1S) K0') / 0.1920

        if proc == self._bs_proc:
            return pu.get_bf('B_s()0 --> J/psi(1S) phi') / 0.1077

        raise ValueError(f'Invalid process {proc}')
    #---------------------------
    def _get_hd_wgt(self, proc : str) -> float:
        '''
        Will return hadronization fractions used as weights
        '''
        if proc in [self._bu_proc, self._bd_proc]:
            return self._fu

        if proc == self._bs_proc:
            return self._fs

        raise ValueError(f'Invalid process: {proc}')
    #---------------------------
    def _get_stats(self, path):
        proc = os.path.basename(path).replace('.json', '')
        df   = pnd.read_json(path)

        return proc, df
    #---------------------------
    def _good_rows(self, r1 : pnd.Series, r2 : pnd.Series) -> bool:
        if {r1.Polarity, r2.Polarity} != {'MagUp', 'MagDown'}:
            log.error('Latest rows are not of opposite polarities')
            return False

        if r1.Events <= 0 or r2.Events <= 0:
            log.error('Either polarity number of events is negative')
            return False

        return True
    #---------------------------
    def _get_st_wgt(self, proc : str) -> float:
        wgt = 1.
        log.warning('Assuming that all inclusive samples have same statistics')

        return wgt
    #---------------------------
    def _get_weight(self, row : pnd.Series) -> float:
        w1 = self._get_st_wgt(row.proc)
        w2 = self._get_hd_wgt(row.proc)
        w3 = self._get_br_wgt(row.proc)

        return w1 * w2 * w3
    #---------------------------
    def get_weights(self) -> pnd.Series:
        '''
        Returns:

        Pandas series with sample weights
        '''
        sr_wgt = self._df.apply(self._get_weight, axis=1)

        return sr_wgt
#---------------------------
