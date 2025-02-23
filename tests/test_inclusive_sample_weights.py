'''
Module with tests for the Reader class of the inclusive_sample_weights module
'''


import numpy
import pandas as pnd

from rx_fitter.inclusive_sample_weights import Reader

#------------------------------------
def _get_df() -> pnd.DataFrame:
    d_data         = {'proc' : [], 'b' : []}
    d_data['proc'] = numpy.random.choice(['Bu_JpsiX_ee_eq_JpsiInAcc', 'Bd_JpsiX_ee_eq_JpsiInAcc', 'Bs_JpsiX_ee_eq_JpsiInAcc'], size=10)
    d_data['b']    = numpy.random.normal(0, 1, size=10)

    return pnd.DataFrame(d_data)
#------------------------------------
def test_simple():
    '''
    Simplest test
    '''
    df            = _get_df()
    obj           = Reader(df)
    df['wgt_sam'] = obj.get_weights()
    df            = df.sort_values(by='proc')

    print(df)
#------------------------------------
