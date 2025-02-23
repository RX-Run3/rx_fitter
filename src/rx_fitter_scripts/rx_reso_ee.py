'''
Script used to fit the resonant mode in the electron channel
'''
from rx_fitter import components as cmp

# ------------------------------
def main():
    '''
    Start here
    '''
    cmp_sig = cmp.get_mc(name = 'Bu_JpsiK_ee_eq_DPC', q2bin='jpsi', model=['cbl', 'cbl', 'cbr'])
    cmp_prc = cmp.get_mc(name = 'Bu_JpsiPi_ee_eq_DPC', q2bin='jpsi', model=['cbl', 'cbl', 'cbr'])
    cmb_cmb = cmp.get_cb(kind='expo')
    cmb_prc = cmp.get_prc()
# ------------------------------
if __name__ == '__main__':
    main()
