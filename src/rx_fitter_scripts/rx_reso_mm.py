'''
Script used to fit the resonant mode in the electron channel
'''
from rx_calibration.hltcalibration.dt_fitter     import DTFitter

from rx_fitter import components as cmp
from rx_fitter import datasets   as dst

# ------------------------------
class Data:
    '''
    Data class
    '''
    cfg = {
            'error_method' : 'minuit_hesse',
            'out_dir'      : 'plots/fit',
            'plotting'     :
            {
                'nbins'   : 50,
                'stacked' : True,
                },
            }
# ------------------------------
def main():
    '''
    Start here
    '''
    cmp_sig = cmp.get_mc(name = 'Bu_JpsiK_mm_eq_DPC' , q2bin='jpsi', model=['cbl', 'cbl', 'cbr'])
    cmp_prc = cmp.get_mc(name = 'Bu_JpsiPi_mm_eq_DPC', q2bin='jpsi', model=['cbl', 'cbl', 'cbr'])
    cmp_cmb = cmp.get_cb(kind='expo')

    rdf     = dst.get_rdf(kind='DATA*', trigger='')

    obj = DTFitter(rdf = rdf, components = [cmp_cmb, cmp_prc, cmp_sig], cfg=Data.cfg)
    obj.fit()
# ------------------------------
if __name__ == '__main__':
    main()
