'''
Script used to do fits to data
'''
from dataclasses         import dataclass

from ROOT                    import RDataFrame
from rx_data.rdf_getter      import RDFGetter
from rx_selection.selection  import load_selection_config

# ---------------------------------
@dataclass
class Data:
    '''
    Class used to share attributes
    '''
    trigger  = 'Hlt2RD_BuToKpMuMu_MVA'

    RDFGetter.samples_dir = '/home/acampove/Data/RX_run3/NO_q2_bdt_mass_Q2_central_VR_v1'

    chanel  : str
    trigger : str
    q2_bin  : str
# ---------------------------------
def _get_rdf() -> RDataFrame:
    gtr = RDFGetter(sample='DATA_24_Mag*_24c*', trigger=Data.trigger)
    rdf = gtr.get_rdf()

    return rdf
# ---------------------------------
def _apply_selection(rdf : RDataFrame) -> RDataFrame:
    cfg = load_selection_config()
    cut = cfg['q2_common'][Data.q2_bin]
    rdf = rdf.Filter(cut)

    return rdf
# ---------------------------------
def _fit():
    pass
# ---------------------------------
def main():
    '''
    Process starts here
    '''
    rdf = _get_rdf()
    rdf = _apply_selection(rdf)

    _fit()
# ---------------------------------
if __name__ == '__main__':
    main()
