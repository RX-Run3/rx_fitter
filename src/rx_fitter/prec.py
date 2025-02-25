'''
Module containing PRec
'''

from ROOT     import RDataFrame
import zfit
import numpy
import pandas as pnd

from zfit.core.basepdf     import BasePDF   as zpdf
from rx_selection          import selection as sel
from rx_data.rdf_getter    import RDFGetter
from dmu.logging.log_store import LogStore

from rx_fitter.inclusive_decays_weights import Reader as inclusive_decays_weights
from rx_fitter.inclusive_sample_weights import Reader as inclusive_sample_weights

log=LogStore.add_logger('rx_fitter:prec')
#-----------------------------------------------------------
class PRec:
    '''
    Class used to calculate the PDF associated to the partially reconstructed background
    '''
    # pylint: disable=too-many-instance-attributes
    #-----------------------------------------------------------
    def __init__(self, samples : list[str], trig : str, q2bin : str, d_weight : dict[str,int]):
        '''
        Parameters:
        -------------------------
        samples (str): MC samples
        trig (str): HLT2 trigger.
        q2bin(str): q2 bin
        d_weight (dict): Dictionary specifying which weights to use, e.g. {'dec' : 1, 'sam' : 1}
        '''

        self._l_sample = samples
        self._trig     = trig
        self._q2bin    = q2bin
        self._d_wg     = d_weight

        self._name     : str
        self._df       : pnd.DataFrame
        self._d_fstat  = {}

        self._nbrem : int = None
        self._d_cut       = None
        self._d_match     = None
        self._initialized = False
    #-----------------------------------------------------------
    def _initialize(self):
        if self._initialized:
            return

        self._check_valid(self._q2bin, ['low', 'central', 'jpsi', 'psi2', 'high'], 'q2bin')
        self._check_weights()

        self._d_match = self._get_match_str()
        self._df      = self._get_df()

        self._initialized = True
    #-----------------------------------------------------------
    @property
    def cuts(self) -> dict[str,str]:
        '''
        Getter for cuts
        '''
        return self._d_cut

    @cuts.setter
    def cuts(self, value : dict[str,str]):
        '''
        These cuts will override whatever is taken from the selection module
        '''
        self._d_cut = value
    #-----------------------------------------------------------
    def _get_df(self) -> pnd.DataFrame:
        '''
        Returns dataframe with masses and weights
        '''
        d_df         = self._get_samples_df()
        d_df         = { sample : self._add_dec_weights(df) for sample, df in d_df.items() }
        df           = pnd.concat(d_df.values(), axis=0)
        df           = self._add_sam_weights(df)

        arr_wgt      = df.wgt_dec.to_numpy() * df.wgt_sam.to_numpy()
        df['wgt_br'] = self._normalize_weights(arr_wgt)

        return df
    #-----------------------------------------------------------
    def _need_var(self, name : str) -> bool:
        if name.endswith('ID'):
            return True

        if name.startswith('B_const_mass'):
            return True

        if name in ['B_M']:
            return True

        return False
    #-----------------------------------------------------------
    def _filter_rdf(self, rdf : RDataFrame, sample : str) -> RDataFrame:
        d_sel = sel.selection(project='RK', analysis='EE', q2bin=self._q2bin, process=sample)
        if self._d_cut is not None:
            log.warning('Overriding default selection')
            d_sel.update(self._d_cut)

        for name, expr in d_sel.items():
            if name == 'mass':
                continue

            rdf = rdf.Filter(expr, name)

        rep = rdf.Report()
        rep.Print()

        return rdf
    #-----------------------------------------------------------
    def _get_samples_df(self) -> dict[str,pnd.DataFrame]:
        '''
        Returns dataframes for each sample
        '''
        d_df = {}
        for sample in self._l_sample:
            gtr        = RDFGetter(sample=sample, trigger=self._trig)
            rdf        = gtr.get_rdf()
            rdf        = self._filter_rdf(rdf, sample)
            l_var      = [ name.c_str() for name in rdf.GetColumnNames() if self._need_var( name.c_str() )]
            data       = rdf.AsNumpy(l_var)
            df         = pnd.DataFrame(data)
            df         = self._filter_by_brem(df)
            df['proc'] = sample

            d_df[sample] = df

        return d_df
    #-----------------------------------------------------------
    def _add_dec_weights(self, df : pnd.DataFrame) -> pnd.DataFrame:
        dec = self._d_wg['dec']

        if   dec == 1:
            log.debug('Adding decay weights')
            df['wgt_dec'] = df.apply(inclusive_decays_weights.read_weight, args=('L1', 'L2', 'H'), axis=1)
        elif dec == 0:
            log.warning('Not using decay weights')
            df['wgt_dec'] = 1.
        else:
            raise ValueError(f'Invalid value of wgt_dec: {dec}')

        arr_wgt      = df.wgt_dec.to_numpy()
        arr_wgt      = self._normalize_weights(arr_wgt)
        df['wgt_dec']= arr_wgt

        return df
    #-----------------------------------------------------------
    def _add_sam_weights(self, df):
        sam = self._d_wg['sam']

        if   sam == 1:
            log.debug('Adding sample weights')
            obj           = inclusive_sample_weights(df)
            df['wgt_sam'] = obj.get_weights()
        elif sam == 0:
            log.warning('Not using sample weights')
            df['wgt_sam'] = 1.
        else:
            raise ValueError(f'Invalid value of wgt_sam: {sam}')

        return df
    #-----------------------------------------------------------
    @property
    def nbrem(self):
        '''
        Number of brem photons
        '''
        return self._nbrem

    @nbrem.setter
    def nbrem(self, value):
        if value not in [0, 1, 2]:
            log.error(f'Invalid nbrem value of: {value}')
            raise ValueError

        self._nbrem = value
    #-----------------------------------------------------------
    def _check_weights(self):
        try:
            [(k1, v1), (k2, v2)] = self._d_wg.items()
        except:
            log.error(f'Cannot extract two weight flags from: {self._d_wg}')
            raise

        if ([k1, k2] != ['dec', 'sam'])  and ([k1, k2] != ['sam', 'dec']):
            raise ValueError(f'Invalid weight keys: {k1}, {k2}')

        if (v1 not in [0, 1]) or (v2 not in [0, 1]):
            raise ValueError(f'Invalid weight values: {v1}, {v2}')
    #-----------------------------------------------------------
    def _check_valid(self, var, l_var, name):
        if var not in l_var:
            log.error(f'Value for {name}, {var}, is not valid')
            raise ValueError
    #-----------------------------------------------------------
    def _get_match_str(self):
        if   self._q2bin == 'jpsi':
            d_match = self._get_match_str_jpsi()
        elif self._q2bin == 'psi2':
            d_match = self._get_match_str_psi2()
        else:
            raise ValueError(f'Invalid q2bin: {self._q2bin}')

        return d_match
    #-----------------------------------------------------------
    def _get_match_str_jpsi(self):
        bd          = '(abs(B_TRUEID) == 511)'
        bp          = '(abs(B_TRUEID) == 521)'
        bs          = '(abs(B_TRUEID) == 531)'

        d_cut                                  = {}
        d_cut[r'$B_d\to c\bar{c}(\to ee)H_s$'] = bd
        d_cut[r'$B^+\to c\bar{c}(\to ee)H_s$'] = bp
        d_cut[r'$B_s\to c\bar{c}(\to ee)H_s$'] = bs

        return d_cut
    #-----------------------------------------------------------
    def _get_match_str_psi2(self):
        bd          = '(abs(B_TRUEID) == 511)'
        bp_psjp     = '(abs(B_TRUEID) == 521) & (abs(Jpsi_TRUEID) == 443) & (abs(Jpsi_MC_MOTHER_ID) == 100443) & (abs(Jpsi_MC_GD_MOTHER_ID) == 521) & (abs(H_MC_MOTHER_ID) == 521)'
        bs          = '(abs(B_TRUEID) == 531)'

        neg_bp_psjp = bp_psjp.replace('==', '!=').replace('&' , '|')
        bp_ex       = f'(abs(B_TRUEID) == 521) & ({neg_bp_psjp})'

        d_cut       = {}
        d_cut[r'$B^+\to \psi(2S)(\to J/\psi X)H_{s}$'] = bp_psjp
        d_cut[r'$B^+\to c\bar{c}(\to ee)H_s$']         = bp_ex
        d_cut[r'$B_d\to c\bar{c}(\to ee)H_s$']         = bd
        d_cut[r'$B_s\to c\bar{c}(\to ee)H_s$']         = bs

        return d_cut
    #-----------------------------------------------------------
    def _get_match_str_psi2_large(self) -> dict[str,str]:
        '''
        Returns dictionary needed to split mix of MC inclusive samples
        '''
        # pylint: disable=too-many-locals

        bp_psjp     = '(abs(Jpsi_MC_MOTHER_ID) == 100443) & (abs(Jpsi_MC_GD_MOTHER_ID) == 521) & (abs(H_MC_MOTHER_ID) == 521)'
        bd_psks     = '(abs(Jpsi_MC_MOTHER_ID) ==    511) & (abs(H_MC_MOTHER_ID) == 313) & (abs(H_MC_GD_MOTHER_ID) == 511) & (abs(Jpsi_TRUEID) == 100443)'
        bp_psks     = '(abs(Jpsi_MC_MOTHER_ID) ==    521) & (abs(H_MC_MOTHER_ID) == 323) & (abs(H_MC_GD_MOTHER_ID) == 521) & (abs(Jpsi_TRUEID) == 100443)'

        neg_bp_psjp = bp_psjp.replace('==', '!=').replace('&' , '|')
        neg_bd_psks = bd_psks.replace('==', '!=').replace('&' , '|')
        neg_bp_psks = bp_psks.replace('==', '!=').replace('&' , '|')

        bp_jpkp     = '(abs(B_TRUEID) == 521) & (abs(H_TRUEID) == 321) & (abs(Jpsi_TRUEID) == 443)'
        bd_jpkp     = '(abs(B_TRUEID) == 511) & (abs(H_TRUEID) == 321) & (abs(Jpsi_TRUEID) == 443)'

        bp_jpkp_ex  = f'({bp_jpkp}) & ({neg_bp_psjp}) & ({neg_bd_psks}) & ({neg_bp_psks})'
        bd_jpkp_ex  = f'({bd_jpkp}) & ({neg_bp_psjp}) & ({neg_bd_psks}) & ({neg_bp_psks})'

        neg_bp_jpkp = bp_jpkp.replace('==', '!=').replace('&' , '|')
        neg_bd_jpkp = bd_jpkp.replace('==', '!=').replace('&' , '|')


        bs          = '(abs(B_TRUEID) == 531)'
        neg_bs      = '(abs(B_TRUEID) != 531)'

        none        = f'({neg_bp_jpkp}) & ({neg_bd_jpkp}) & ({neg_bp_psjp}) & ({neg_bd_psks}) & ({neg_bp_psks}) & ({neg_bs})'

        d_cut            = {}
        d_cut['bp_psjp'] = bp_psjp
        d_cut['bp_psks'] = bp_psks
        d_cut['bp_jpkp'] = bp_jpkp_ex

        d_cut['bd_psks'] = bd_psks
        d_cut['bd_jpkp'] = bd_jpkp_ex

        d_cut['bs']      = bs

        d_cut['unmatched'] = none

        return d_cut
    #-----------------------------------------------------------
    def _get_match_str_psi2_all(self) -> dict[str,str]:
        d_cut           = {}
        d_cut['jpsi']   = '(Jpsi_TRUEID == 443)'
        d_cut['nojpsi'] = '(Jpsi_TRUEID != 443)'

        return d_cut
    #-----------------------------------------------------------
    def _filter_by_brem(self, df):
        if self._nbrem is None:
            return df

        log.debug(f'Applying nbrem = {self._nbrem} requirement')
        df = df[df.nbrem == self._nbrem] if self._nbrem < 2 else df[df.nbrem >= 2]
        df = df.reset_index(drop=True)

        return df
    #-----------------------------------------------------------
    def _print_wgt_stat(self, arr_wgt):
        l_wgt = arr_wgt.tolist()
        s_wgt = set(l_wgt)

        log.debug('-' * 20)
        log.debug(f'{"Frequency":<10}{"Weight":>10}')
        for wgt in s_wgt:
            nwgt = numpy.count_nonzero(wgt == arr_wgt)
            log.debug(f'{nwgt:<10}{wgt:>10.3}')
    #-----------------------------------------------------------
    def _normalize_weights(self, arr_wgt):
        tot_wgt = arr_wgt.sum()
        num_wgt = arr_wgt.shape[0]
        fact    = num_wgt / tot_wgt
        arr_wgt = fact * arr_wgt

        return arr_wgt
    #-----------------------------------------------------------
    def _get_df_id(self, df : pnd.DataFrame) -> pnd.DataFrame:
        l_col = [
                'L1_TRUEID',
                'L2_TRUEID',
                'Jpsi_TRUEID',
                'Jpsi_MC_MOTHER_ID',
                'Jpsi_MC_GD_MOTHER_ID',
                'H_TRUEID',
                'H_MC_MOTHER_ID',
                'H_MC_GD_MOTHER_ID',
                'B_TRUEID'
                ]

        df = df[l_col]

        return df.reset_index(drop=True)
    #-----------------------------------------------------------
    def _filter_mass(self, df : pnd.DataFrame, mass : str, obs):
        ([[minx]], [[maxx]]) = obs.limits

        cut   = f'({minx} < {mass}) & ({mass} < {maxx})'
        log.debug(f'Applying: {cut}')
        inum  = df.shape[0]
        df    = df.query(cut)
        fnum  = df.shape[0]

        self._d_fstat[cut] = inum, fnum

        return df
    #-----------------------------------------------------------
    def _filter_cut(self, cut : str) -> pnd.DataFrame:
        if cut is None:
            return self._df

        log.info(f'Applying cut: {cut}')
        inum = self._df.shape[0]
        df   = self._df.query(cut)
        fnum = df.shape[0]

        self._d_fstat[cut] = inum, fnum

        return df
    #-----------------------------------------------------------
    def _get_pdf(self, mass : str, cut : str, **kwargs) -> zpdf:
        '''
        Will take the mass, with values in:

        mass: Non constrained B mass
        mass_jpsi: Jpsi constrained B mass
        mass_psi2: Psi2S constrained B mass

        The observable.

        Optional arguments:
        Cut

        **kwargs: These are all arguments for KDE1DimFFT

        and it will return a KDE1DimFFT PDF.
        '''
        df = self._filter_cut(cut)
        df = self._filter_mass(df, mass, kwargs['obs'])

        log.info(f'Using mass: {mass} for component {kwargs["name"]}')
        arr_mass = df[mass].to_numpy()
        arr_wgt  = df.wgt_br.to_numpy()
        df_id    = self._get_df_id(df)

        self._print_cutflow()

        pdf          = zfit.pdf.KDE1DimFFT(arr_mass, weights=arr_wgt, **kwargs)
        pdf.arr_mass = arr_mass
        pdf.arr_wgt  = arr_wgt
        pdf.df_id    = df_id

        return pdf
    #-----------------------------------------------------------
    def _print_cutflow(self) -> None:
        log.debug('-' * 50)
        log.debug(f'{"Cut":<30}{"Total":<20}{"Passed":<20}')
        log.debug('-' * 50)
        for cut, (inum, fnum) in self._d_fstat.items():
            log.debug(f'{cut:<30}{inum:<20}{fnum:<20}')
        log.debug('-' * 50)
    #-----------------------------------------------------------
    def get_sum(self, mass : str, name='unnamed', **kwargs) -> zpdf:
        '''Provides extended PDF that is the sum of multiple KDEs representing PRec background

        Parameters:
        mass (str) : Defines which mass constrain to use, choose between "B_M", "B_const_mass_M", "B_const_mass_psi2S_M"
        name (str) : PDF name
        **kwargs: Arguments meant to be taken by zfit KDE1DimFFT

        Returns:
        zfit.pdf.SumPDF instance
        '''
        self._name = name
        self._initialize()

        d_pdf     = { name : self._get_pdf(mass, cut, name=name, **kwargs) for name, cut in self._d_match.items()}
        l_pdf     = list(d_pdf.values())
        l_wgt_yld = [ sum(pdf.arr_wgt) for pdf in l_pdf ]
        l_frc     = [ wgt_yld / sum(l_wgt_yld) for wgt_yld in l_wgt_yld ]
        l_yld     = [ zfit.param.Parameter(f'f_{pdf.name}', frc, 0, 1) for pdf, frc in zip(l_pdf, l_frc)]
        for yld in l_yld:
            yld.floating = False
        l_df_id   = [ pdf.df_id for pdf in l_pdf ]

        pdf          = zfit.pdf.SumPDF(l_pdf, fracs=l_yld)
        nor          = zfit.param.Parameter('nprc', sum(l_wgt_yld), 0, 1000_000)
        pdf          = pdf.create_extended(nor, name=self._name)

        l_arr_mass   = [ pdf.arr_mass for pdf in l_pdf ]
        l_arr_wgt    = [ pdf.arr_wgt  for pdf in l_pdf ]

        pdf.arr_mass = numpy.concatenate(l_arr_mass)
        pdf.arr_wgt  = numpy.concatenate(l_arr_wgt )
        pdf.df_id    = pnd.concat(l_df_id, ignore_index=True)

        return pdf
#-----------------------------------------------------------
