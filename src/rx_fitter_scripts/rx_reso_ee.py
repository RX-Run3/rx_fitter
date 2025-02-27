'''
Script used to fit the resonant mode in the electron channel
'''
import zfit
from rx_fitter import components as cmp

# ------------------------------
def main():
    '''
    Start here
    '''

    obs=zfit.Space('B_const_mass_M', limits=(4500, 6000))
    trigger = 'Hlt2RD_BuToKpEE_MVA'

    cmp_sig = cmp.get_mc(obs = obs, sample = 'Bu_JpsiK_ee_eq_DPC' , trigger, q2bin='jpsi') 
    cmp_csp = cmp.get_mc(obs = obs, sample = 'Bu_JpsiPi_ee_eq_DPC', trigger, q2bin='jpsi')
    cmb_prc = cmp.get_prc(obs = obs, trigger=trigger, q2bin='jpsi')
    cmb_cmb = cmp.get_cb(kind='expo')
# ------------------------------
if __name__ == '__main__':
    main()
