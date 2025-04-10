'''
Script holding ConstraintReader class
'''
from dmu.logging.log_store   import LogStore
from rx_fitter.signal_scales import FitParameters

log=LogStore.add_logger('rx_fitter:constraint_reader')
# -------------------------------------------------------------
class ConstraintReader:
    '''
    Class meant to provide constraints for fitting model
    '''
    # -------------------------------------------------------------
    def __init__(self, parameters : list[str], mva_cut : str):
        '''
        Parameters: List of parameter names as in the PDF
        mva_cut   : Needed to override default MVA selection
        '''

        self._l_par   = parameters
        self._mva_cut = mva_cut
        self._d_const = {}
    # -------------------------------------------------------------
    def _add_signal_constraints(self) -> None:
        l_par = [ par for par in self._l_par if 'Signal' in par ]
        obj = FitParameters()

        for par in l_par:
            val, err = obj.get_parameter_scale(name=par)

            self._d_const[par] = val, err
    # -------------------------------------------------------------
    def _add_brem_constraints(self) -> None:
        pass
    # -------------------------------------------------------------
    def _add_prec_constraints(self) -> None:
        pass
    # -------------------------------------------------------------
    def get_constraints(self) -> dict[str,tuple[float,float]]:
        '''
        Returns dictionary with constraints, i.e.

        Key  : Name of fitting parameter
        Value: Tuple with mu and error
        '''
        self._add_signal_constraints()
        self._add_brem_constraints()
        self._add_prec_constraints()

        return self._d_const
# -------------------------------------------------------------
