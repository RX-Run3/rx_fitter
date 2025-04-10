'''
Script holding ConstraintReader class
'''
from dmu.logging.log_store       import LogStore
from rx_efficiencies.decay_names import DecayNames as dn
from rx_fitter.signal_scales     import FitParameters
from rx_fitter.prec_scales       import PrecScales

log=LogStore.add_logger('rx_fitter:constraint_reader')
# -------------------------------------------------------------
class ConstraintReader:
    '''
    Class meant to provide constraints for fitting model
    '''
    # -------------------------------------------------------------
    def __init__(self, parameters : list[str], mva_cut : str, q2bin : str):
        '''
        Parameters: List of parameter names as in the PDF
        mva_cut   : Needed to override default MVA selection
        q2bin     : q2 bin
        '''

        self._l_par   = parameters
        self._mva_cut = mva_cut
        self._q2bin   = q2bin

        self._d_const = {}
        self._signal  = 'bpkpee' # This is the signal decay nickname, needed for PRec scales constraints
    # -------------------------------------------------------------
    def _add_signal_constraints(self) -> None:
        obj = FitParameters()

        for par in self._l_par:
            if 'Signal' in par:
                val, err = obj.get_parameter_scale(name=par)
            elif par.startswith('frac_brem_'):
                val, err = obj.get_brem_fraction(name=par)
            else:
                continue

            log.debug(f'Adding constrint for: {par}')

            self._d_const[par] = val, err
    # -------------------------------------------------------------
    def _proc_from_par(self, par_name : str) -> str:
        sample = par_name[1:] # Parameter name is expected to look like sSAMPLE_NICKNAME
        decay  = dn.nic_from_sample(sample)

        return decay
    # -------------------------------------------------------------
    def _add_prec_constraints(self) -> None:
        d_cut    = {'mva' : self._mva_cut}

        for par in self._l_par:
            if 'Signal' in par:
                continue

            if par.startswith('frac_brem_'):
                continue

            log.debug(f'Adding constrint for: {par}')

            process  = self._proc_from_par(par)
            obj      = PrecScales(proc=process, q2bin=self._q2bin, d_cut=d_cut)
            val, err = obj.get_scale(signal=self._signal)

            self._d_const[par] = val, err
    # -------------------------------------------------------------
    def get_constraints(self) -> dict[str,tuple[float,float]]:
        '''
        Returns dictionary with constraints, i.e.

        Key  : Name of fitting parameter
        Value: Tuple with mu and error
        '''
        self._add_signal_constraints()
        self._add_prec_constraints()

        return self._d_const
# -------------------------------------------------------------
