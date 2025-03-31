'''
Module with fitting models 
'''

import zfit
from dmu.stats.zfit_models  import HypExp
from dmu.stats.zfit_models  import ModExp

from zfit.core.interfaces   import ZfitSpace as zobs
from zfit.core.basepdf      import BasePDF   as zpdf

# ---------------------------------------------
def _get_hypexp(obs : zobs) -> zpdf:
    mu = zfit.Parameter('mu',  5000,  4000,  6000)
    ap = zfit.Parameter('ap', 0.020,     0,  0.10)
    bt = zfit.Parameter('bt', 0.002, 0.001, 0.003)

    pdf= HypExp(obs=obs, mu=mu, alpha=ap, beta=bt)

    return pdf
# ---------------------------------------------
def _get_modexp(obs : zobs) -> zpdf:
    mu = zfit.Parameter('mu',  4500,  4000,  6000)
    ap = zfit.Parameter('ap', 0.020,     0,   0.1)
    bt = zfit.Parameter('bt', 0.002, 0.001, 0.005)

    pdf= ModExp(obs=obs, mu=mu, alpha=ap, beta=bt)

    return pdf
# ---------------------------------------------
def get_pdf(obs : zobs, name : str) -> zpdf:
    '''
    Function returning a zfit PDF from observable and name.
    Raises NotImplementedError if PDF is missing
    '''
    if name == 'HypExp':
        return _get_hypexp(obs=obs)

    if name == 'ModExp':
        return _get_modexp(obs=obs)

    raise NotImplementedError(f'Cannot find {name} PDF')
# ---------------------------------------------
