'''
Script used to fit the resonant mode in the electron channel
'''
import zfit
from dmu.logging.log_store                   import LogStore
from rx_data.rdf_getter                      import RDFGetter
from rx_calibration.hltcalibration.dt_fitter import DTFitter
from rx_fitter                               import components as cmp

log = LogStore.add_logger('rx_fitter:rx_reso_ee')
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
                'd_leg'   : {
                    'SumPDF_ext'     : r'$B^+\to K^+\mu^+\mu^-$',
                    'Exponential_ext': 'Combinatorial',
                    }
                },
            }
# ------------------------------
def main():
    '''
    Start here
    '''

    trigger = 'Hlt2RD_BuToKpEE_MVA'
    q2bin   = 'jpsi'
    nbrem   = 1
    obs     = zfit.Space('B_const_mass_M', limits=(5000, 6000))

    cmp_sig = cmp.get_mc(obs = obs, sample = 'Bu_JpsiK_ee_eq_DPC' , trigger=trigger, q2bin=q2bin, nbrem=nbrem)
    cmp_csp = cmp.get_mc(obs = obs, sample = 'Bu_JpsiPi_ee_eq_DPC', trigger=trigger, q2bin=q2bin, nbrem=nbrem)
    cmp_prc = cmp.get_prc(obs= obs, trigger=trigger, q2bin=q2bin, nbrem=nbrem)
    cmp_cmb = cmp.get_cb(obs = obs, kind='exp')

    gtr = RDFGetter(sample='DATA*', trigger=trigger)
    rdf = gtr.get_rdf()

    obj = DTFitter(rdf = rdf, components = [cmp_cmb, cmp_prc, cmp_csp, cmp_sig], cfg=Data.cfg)
    obj.fit()
# ------------------------------
if __name__ == '__main__':
    main()
