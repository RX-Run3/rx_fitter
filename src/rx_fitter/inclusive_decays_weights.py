'''
Module with Reader class
'''
# pylint: disable=too-many-instance-attributes, invalid-name

import pandas as pnd
from rx_fitter.pchain import PChain

#---------------------------
class Reader:
    '''
    Class used to attach decay weights to dataframe of inclusive decays
    '''
    def __init__(self, l1, l2, kp):
        self._l1_ch = l1
        self._l2_ch = l2
        self._kp_ch = kp

        self._weight= 1.0

        self._Jp_id = 443
        self._Ps_id = 100443
        self._Pi_id = 211
        self._Ph_id = 333
        self._Et_id = 221
        self._Ks_id = 313
        self._Kp_id = 321
        self._KS_id = 310

        self._Bd_id = 511
        self._Bu_id = 521
        self._Bs_id = 531

        #From https://gitlab.cern.ch/LHCb-RD/ewp-rkstz/-/blob/master/analysis/kernel/inc/ConstDef.hpp
        self._Bc         = 541
        self._Lb         = 5122
        self._L0         = 3122
        self._Kst_c      = 323
        self._K_1_1270_z = 10313
        self._K_1_1270_c = 10323
        self._K_2_1430_c = 325
        self._K_2_1430_z = 315
        self._K_1_1400_z = 20313
        self._K_1_1400_c = 20323
        self._K0         = 311
        self._KLong      = 130
        self._Pi0        = 111
        self._M          = 13
        self._E          = 11
        self._Tau        = 15
        self._P          = 2212
        self._N          = 2112
        self._Eta_prime  = 331
        self._Rho0       = 113
        self._Rho_c      = 213
        self._Omega      = 223
        self._D0         = 421
        self._Dp         = 411
        self._Ds         = 431
        self._Bst0       = 513
        self._Bst_plus   = 523
        self._Photon     = 22

        self._l_bid = [self._Bd_id, self._Bu_id, self._Bs_id]
    #---------------------------
    def _get_jpsi_wgt(self) -> float:
        weight = 1.0
        for bid in  self._l_bid:
            flg_1 = self._l1_ch.MatchUpstream( self._Ps_id, bid)
            flg_2 = self._l2_ch.MatchUpstream( self._Ps_id, bid)

            if not flg_1 and not flg_2:
                continue

            flg_3 = self._l1_ch.MatchMother(self._Jp_id)
            flg_4 = self._l2_ch.MatchMother(self._Jp_id)
            flg_5 = self._l1_ch.MatchMother(self._Ps_id)
            flg_6 = self._l2_ch.MatchMother(self._Ps_id)

            if   flg_3 or flg_4:
                #weight = 0.6254 / ( 1-0.1741) #0.75
                weight = 0.958
            elif flg_5 or flg_6:
                #weight = 1.3200 / 0.1741 #7.58
                weight = 0.771

        return weight
    #---------------------------
    def _get_brfrac_corr(self) -> float:
        weight= 1
        flg_1 = self._l1_ch.MatchUpstream(self._Ps_id, self._Bd_id)
        flg_2 = self._l2_ch.MatchUpstream(self._Ps_id, self._Bd_id)
        if flg_1 or flg_2:
            weight = 1.17

        flg_1 = self._l1_ch.MatchUpstream(self._Ps_id, self._Bu_id)
        flg_2 = self._l2_ch.MatchUpstream(self._Ps_id, self._Bu_id)
        if flg_1 or flg_2:
            weight = 1.35

        return weight
    #---------------------------
    def _either_track_has(self, pid : int) -> bool:
        flg_l1 = self._l1_ch.HasInChain(pid)
        flg_l2 = self._l2_ch.HasInChain(pid)
        flg_kp = self._kp_ch.HasInChain(pid)

        return flg_l1 or flg_l2 or flg_kp
    #---------------------------
    def _get_psi_over_jpsi(self) -> float:
        flg_ps = self._either_track_has(self._Ps_id)
        if not flg_ps:
            return 1.0

        flg_bu = self._either_track_has(self._Bu_id)
        flg_bd = self._either_track_has(self._Bd_id)
        flg_bs = self._either_track_has(self._Bs_id)

        if flg_bs:
            return ( 5.40E-4/ 1.0800E-3 ) / ( 0.0748/ 0.1077) #0.72

        if flg_bu:
            return ( 6.19E-4/ 1.0006E-3 ) / ( 0.0729/ 0.1595) #1.35

        if flg_bd:
            return ( 5.90E-4/ 1.2700E-3 ) / (0.07610/ 0.1850) #1.13

        return 1.0
    #---------------------------
    def _get_kst_wgt(self) -> float:
        weight = 1.0
        if self._kp_ch.MatchDecay( [self._Pi_id, self._KS_id             ]):
            weight *= 0.5

        if self._kp_ch.MatchDecay( [self._Kp_id, self._Ph_id, self._Bd_id]):
            weight *= 0.5/0.9974

        if self._kp_ch.MatchDecay( [self._Kp_id, self._Ph_id, self._Bu_id]):
            weight *= 0.5/0.7597

        if self._kp_ch.MatchDecay( [self._Pi_id, self._Et_id] ):
            weight *= 0.28/0.4

        if self._kp_ch.MatchDecay( [self._Kp_id, self._Ks_id] ):
            weight *= 0.66 / 0.7993

        if self._kp_ch.MatchDecay( [self._Kp_id, self._Kst_c] ):
            weight *= 0.33 / 0.4993

        if self._kp_ch.MatchDecay( [self._Kp_id, self._K_2_1430_c]):
            weight *= 0.1670/0.2485

        return weight
    #---------------------------
    def get_weight(self) -> float:
        '''
        Returns:

        wt (float): Weight for candidate
        '''
        w1 = self._get_jpsi_wgt()
        w2 = self._get_brfrac_corr()
        w3 = self._get_psi_over_jpsi()
        w4 = self._get_kst_wgt()

        wt = w1 * w2 * w3 * w4

        return wt
    #---------------------------
    @staticmethod
    def get_chain(name, row : pnd.Series) -> PChain:
        '''
        Will return an instance of a PChain for a given dataframe row
        '''
        v1 = row[f'{name}_TRUEID']
        v2 = row[f'{name}_MC_MOTHER_ID']
        v3 = row[f'{name}_MC_GD_MOTHER_ID']
        v4 = row[f'{name}_MC_GD_GD_MOTHER_ID']

        chain = PChain(v1, v2, v3, v4)

        return chain
    #---------------------------
    @staticmethod
    def read_weight(row : pnd.Series, p1 : str, p2 : str, p3 : str) -> float:
        '''
        This method will return the BR weights

        Parameters
        -------------------
        row: Row in pandas dataframe holding PDG id information for the three particles.
        pX: Name of particle in dataframe, e.g. L1, L2, H

        Returns
        -------------------
        wgt (float): Weight for candidate associated to input row
        '''
        p1_ch = Reader.get_chain(p1, row)
        p2_ch = Reader.get_chain(p2, row)
        p3_ch = Reader.get_chain(p3, row)

        obj = Reader(p1_ch, p2_ch, p3_ch)
        wgt = obj.get_weight()

        return wgt
#---------------------------
