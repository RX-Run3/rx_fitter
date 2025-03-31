'''
Module with fitting models 
'''

import zfit
from dmu.stats.zfit_models  import HypExp
from dmu.stats.zfit_models  import ModExp

from zfit.core.interfaces   import ZfitSpace as zobs
from zfit.core.basepdf      import BasePDF   as zpdf

# ---------------------------------------------
def _get_suj(obs : zobs) -> zpdf:
    mu  = zfit.Parameter('mu', 5000, 4000, 6000)
    lb  = zfit.Parameter('lb',  100,   10, 1000)
    dl  = zfit.Parameter('dl',  2.5,    1,   10)
    gm  = zfit.Parameter('gm',  -10,  -20,   20)

    dl.floating = False
    gm.floating = False

    pdf = zfit.pdf.JohnsonSU(mu=mu, lambd=lb, gamma=gm, delta=dl, obs=obs, name='SUJohnson')

    return pdf
# ---------------------------------------------
def _get_pol2(obs : zobs) -> zpdf:
    a   = zfit.Parameter('a', -0.005, -0.95, 0.00)
    b   = zfit.Parameter('b',  0.000, -0.95, 0.95)
    pdf = zfit.pdf.Chebyshev(obs=obs, coeffs=[a, b], name='Chebyshev 2nd')

    return pdf
# ---------------------------------------------
def _get_pol3(obs : zobs) -> zpdf:
    a   = zfit.Parameter('a', -0.005, -0.95, 0.00)
    b   = zfit.Parameter('b',  0.000, -0.95, 0.95)
    c   = zfit.Parameter('c',  0.000, -0.95, 0.95)
    pdf = zfit.pdf.Chebyshev(obs=obs, coeffs=[a, b, c], name='Chebyshev 3rd')

    return pdf
# ---------------------------------------------
def _get_exponential(obs : zobs) -> zpdf:
    c  = zfit.Parameter('c', -0.002, -0.003, 0.0)
    pdf= zfit.pdf.Exponential(obs=obs, lam=c)

    return pdf
# ---------------------------------------------
def _get_hypexp(obs : zobs) -> zpdf:
    mu = zfit.Parameter('mu',  5000,   4000,  6000)
    ap = zfit.Parameter('ap', 0.020,      0,  0.10)
    bt = zfit.Parameter('bt', 0.002, 0.0001, 0.003)

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

    if name == 'Exp':
        return _get_exponential(obs=obs)

    if name == 'Pol2':
        return _get_pol2(obs=obs)

    if name == 'Pol3':
        return _get_pol3(obs=obs)

    if name == 'SUJohnson':
        return _get_suj(obs=obs)

    raise NotImplementedError(f'Cannot find {name} PDF')
# ---------------------------------------------
