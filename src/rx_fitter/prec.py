'''
Module containing PRec
'''
import os
import copy
import json
from typing     import Union
from contextlib import contextmanager

import numpy
import pandas            as pnd
import matplotlib.pyplot as plt

from dmu.stats.zfit         import zfit
from dmu.generic            import hashing
from dmu.logging.log_store  import LogStore
from dmu.stats.zfit_plotter import ZFitPlotter
from dmu.stats.utilities    import is_pdf_usable
from dmu.stats              import utilities as sut

from zfit.core.parameter   import Parameter as zpar
from zfit.core.basepdf     import BasePDF   as zpdf
from rx_selection          import selection as sel
from rx_data.rdf_getter    import RDFGetter
from ROOT                  import RDataFrame

from rx_fitter.inclusive_decays_weights import Reader as inclusive_decays_weights
from rx_fitter.inclusive_sample_weights import Reader as inclusive_sample_weights

log=LogStore.add_logger('rx_fitter:prec')
#-----------------------------------------------------------
class PRec:
    '''
    Class used to calculate the PDF associated to the partially reconstructed background
    '''
    use_cache = True # Use cached if found
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
        self._d_wg     = copy.deepcopy(d_weight)

        self._name     : str
        self._df       : pnd.DataFrame
        self._d_fstat  = {}

        self._d_match         = self._get_match_str()
        self._l_mass          = ['B_Mass', 'B_Mass_smr', 'B_const_mass_M', 'B_const_mass_psi2S_M']
        self._min_entries     = 40 # Will not build KDE if fewer entries than this are found
        self._min_isj_entries = 500 #if Fewer entries than this, switch from ISJ to FFT
        self._initialized     = False
    #-----------------------------------------------------------
    def _initialize(self):
        if self._initialized:
            return

        self._check_valid(self._q2bin, ['low', 'central', 'jpsi', 'psi2', 'high'], 'q2bin')
        self._check_weights()

        self._df      = self._get_df()

        self._initialized = True
    #-----------------------------------------------------------
    def _get_df(self) -> pnd.DataFrame:
        '''
        Returns dataframe with masses and weights
        '''
        d_df         = self._get_samples_df()
        d_df         = { sample : self._add_dec_weights(sample, df) for sample, df in d_df.items() }
        df           = pnd.concat(d_df.values(), axis=0)
        df           = self._add_sam_weights(df)

        if len(df) == 0:
            return df

        arr_wgt      = df.wgt_dec.to_numpy() * df.wgt_sam.to_numpy()
        df['wgt_br'] = self._normalize_weights(arr_wgt)

        return df
    #-----------------------------------------------------------
    def _need_var(self, name : str) -> bool:
        needed = False

        if name.endswith('ID'):
            needed = True

        if name in self._l_mass:
            needed = True

        if needed:
            log.debug(f'Picking up {name}')

        return needed
    #-----------------------------------------------------------
    def _filter_rdf(self, rdf : RDataFrame, sample : str) -> RDataFrame:
        d_sel = sel.selection(trigger=self._trig, q2bin=self._q2bin, process=sample)
        for name, expr in d_sel.items():
            if name == 'mass':
                continue

            log.debug(f'{name:<20}{expr}')
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
            df['proc'] = sample

            d_df[sample] = df

        return d_df
    #-----------------------------------------------------------
    def _add_dec_weights(self, sample : str, df : pnd.DataFrame) -> pnd.DataFrame:
        if len(df) == 0:
            return df

        dec = self._d_wg['dec']

        if   dec == 1:
            log.debug(f'Adding decay weights to: {sample}')
            df['wgt_dec'] = df.apply(inclusive_decays_weights.read_weight, args=('L1', 'L2', 'H'), axis=1)
        elif dec == 0:
            log.warning(f'Not using decay weights in: {sample}')
            df['wgt_dec'] = 1.
        else:
            raise ValueError(f'Invalid value of wgt_dec: {dec}')

        arr_wgt      = df.wgt_dec.to_numpy()
        arr_wgt      = self._normalize_weights(arr_wgt)
        df['wgt_dec']= arr_wgt

        return df
    #-----------------------------------------------------------
    def _add_sam_weights(self, df : pnd.DataFrame) -> pnd.DataFrame:
        if len(df) == 0:
            return df

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

        arr_wgt      = df.wgt_sam.to_numpy()
        arr_wgt      = self._normalize_weights(arr_wgt)
        df['wgt_sam']= arr_wgt

        return df
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
    def _get_match_str(self) -> dict[str,str]:
        if   self._q2bin == 'jpsi':
            d_match = self._get_match_str_jpsi()
        elif self._q2bin == 'psi2':
            d_match = self._get_match_str_psi2()
        elif self._q2bin in ['central', 'high']:
            d_match = self._get_match_str_psi2()
        else:
            raise ValueError(f'Invalid q2bin: {self._q2bin}')

        return d_match
    #-----------------------------------------------------------
    def _get_match_str_jpsi(self) -> dict[str,str]:
        bd          = '(abs(B_TRUEID) == 511)'
        bp          = '(abs(B_TRUEID) == 521)'
        bs          = '(abs(B_TRUEID) == 531)'

        d_cut                                  = {}
        d_cut[r'$B_d\to c\bar{c}(\to ee)H_s$'] = bd
        d_cut[r'$B^+\to c\bar{c}(\to ee)H_s$'] = bp
        d_cut[r'$B_s\to c\bar{c}(\to ee)H_s$'] = bs

        return d_cut
    #-----------------------------------------------------------
    def _get_match_str_psi2(self) -> dict[str,str]:
        bd          = '(abs(B_TRUEID) == 511)'
        bp_psjp     = '(abs(B_TRUEID) == 521) & (abs(Jpsi_TRUEID) == 443) & (abs(Jpsi_MC_MOTHER_ID) == 100443) & (abs(Jpsi_MC_GD_MOTHER_ID) == 521) & (abs(H_MC_MOTHER_ID) == 521)'
        bs          = '(abs(B_TRUEID) == 531)'

        neg_bp_psjp = bp_psjp.replace('==', '!=').replace('&' , '|')
        bp_ex       = f'(abs(B_TRUEID) == 521) & ({neg_bp_psjp})'

        d_cut       = {}
        d_cut[r'$B^+\to \psi(2S)(\to J/\psi+X)H_{s}$'] = bp_psjp
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
            log.debug('Not applying any cut')
            return self._df

        log.info(f'Applying cut: {cut}')
        inum = self._df.shape[0]
        df   = self._df.query(cut)
        fnum = df.shape[0]

        self._d_fstat[cut] = inum, fnum

        return df
    #-----------------------------------------------------------
    def _get_identifier(self, mass : str, cut : str, **kwargs) -> str:
        cwargs = copy.deepcopy(kwargs)
        del cwargs['obs']

        swgt = json.dumps(self._d_wg , sort_keys=True)
        scwg = json.dumps(cwargs     , sort_keys=True)

        l_d_sel   = [ sel.selection(trigger=self._trig, q2bin=self._q2bin, process=sample) for sample in self._l_sample ]
        l_element = [
                swgt,
                self._trig,
                self._q2bin,
                mass,
                scwg,           # Stringified keyword arguments
                self._l_sample, # ccbar cocktail sample names
                l_d_sel,        # list of selections, one for each ccbar cocktail sample
                cut]

        hsh  = hashing.hash_object(l_element)

        return hsh
    #-----------------------------------------------------------
    def _path_from_identifier(self, identifier : str) -> str:
        dir_path = '/tmp/cache/prec'
        os.makedirs(dir_path, exist_ok=True)

        return f'{dir_path}/pdf_{identifier}.json'
    #-----------------------------------------------------------
    def _drop_before_saving(self, df : pnd.DataFrame) -> pnd.DataFrame:
        l_needed = self._l_mass + ['wgt_br', 'wgt_dec', 'wgt_sam']
        l_drop   = [ name for name in df.columns if  name not in l_needed ]
        df       = df.drop(l_drop, axis=1)

        return df
    #-----------------------------------------------------------
    def _get_pdf(self, mass : str, cut : str, **kwargs) -> Union[zpdf,None]:
        '''
        Will take the mass, with values in:

        mass: Non constrained B mass
        mass_jpsi: Jpsi constrained B mass
        mass_psi2: Psi2S constrained B mass

        The observable.

        Optional arguments:
        Cut

        **kwargs: These are all arguments for KDE1DimISJ or KDE1DimFFT
        '''
        identifier = self._get_identifier(mass, cut, **kwargs)
        cache_path = self._path_from_identifier(identifier)

        if os.path.isfile(cache_path) and PRec.use_cache:
            log.warning(f'Cached PDF found, loading: {cache_path}')
            log.debug(f'Cut: {cut}')
            df = pnd.read_json(cache_path)
            if len(df) == 0:
                return None
        else:
            self._initialize()
            if len(self._df) == 0:
                return None

            if PRec.use_cache:
                log.info('Cached PDF not found, calculating it')
            else:
                log.warning('Caching turned off, recalculating PDF')

            df = self._filter_cut(cut)
            df = self._filter_mass(df, mass, kwargs['obs'])
            log.info(f'Using mass: {mass} for component {kwargs["name"]}')
            self._print_cutflow()
            df=self._drop_before_saving(df)
            df.to_json(cache_path, indent=4)

        arr_mass     = df[mass].to_numpy()
        nentries     = len(arr_mass)
        name         = kwargs['name']
        if nentries < self._min_entries:
            log.warning(f'Found fewer than {self._min_entries}: {nentries}, skipping PDF {name}')
            return None

        log.info(f'Building PDF with {nentries} entries for {name}')

        if nentries < self._min_isj_entries:
            log.info('Using FFT KDE for low statistics sample')
            pdf      = zfit.pdf.KDE1DimFFT(arr_mass, weights=df.wgt_br.to_numpy(), **kwargs)
        else:
            log.info('Using ISJ KDE for high statistics sample')
            if 'bandwidth' in kwargs: # ISJ does not accept this argument
                del kwargs['bandwidth']

            pdf      = zfit.pdf.KDE1DimISJ(arr_mass, weights=df.wgt_br.to_numpy(), **kwargs)

        pdf.arr_mass = arr_mass
        pdf.arr_wgt  = df.wgt_br.to_numpy()
        pdf.arr_sam  = df.wgt_sam.to_numpy()
        pdf.arr_dec  = df.wgt_dec.to_numpy()

        if not is_pdf_usable(pdf):
            return None

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
    def _frac_from_pdf(self, pdf : zpdf, frc : float) -> zpar:
        name = pdf.name
        name = name.replace(r' ', '_')
        name = name.replace(r'$', '_')
        name = name.replace(r'{', '_')
        name = name.replace(r'}', '_')
        name = name.replace(r'(', '_')
        name = name.replace(r')', '_')
        name = name.replace(r'>', '_')
        name = name.replace(r'-', '_')
        name = name.replace('\\', '_')
        name = name.replace(r'/', '_')
        name = name.replace(r'+', '_')
        name = name.replace(r'^', '_')

        par  = zfit.param.Parameter(f'f_{name}', frc, 0, 1)

        return par
    #-----------------------------------------------------------
    def get_sum(self, mass : str, name='unnamed', **kwargs) -> Union[zpdf,None]:
        '''Provides extended PDF that is the sum of multiple KDEs representing PRec background

        Parameters:
        mass (str) : Defines which mass constrain to use, choose between "B_M", "B_const_mass_M", "B_const_mass_psi2S_M"
        name (str) : PDF name
        **kwargs: Arguments meant to be taken by zfit KDE1DimFFT

        Returns:
        zfit.pdf.SumPDF instance
        '''
        self._name = name

        # These cuts are not meant to override the selection, they are used to classify the fully selected data
        # into physically meaningful categories, such that they can be used later, together
        d_pdf     = { name : self._get_pdf(mass, cut, name=name, **kwargs) for name, cut in self._d_match.items()}
        d_pdf     = { name : pdf                                           for name, pdf in d_pdf.items() if pdf is not None}

        l_pdf     = list(d_pdf.values())
        l_wgt_yld = [ sum(pdf.arr_wgt) for pdf in l_pdf ]
        l_frc     = [ wgt_yld / sum(l_wgt_yld) for wgt_yld in l_wgt_yld ]
        l_yld     = [ self._frac_from_pdf(pdf=pdf, frc=frc) for pdf, frc in zip(l_pdf, l_frc)]
        for yld in l_yld:
            yld.floating = False

        if   len(l_pdf) >= 2:
            pdf   = zfit.pdf.SumPDF(l_pdf, fracs=l_yld, name='ccbar PRec')
        elif len(l_pdf) == 1:
            [pdf] = l_pdf
        else:
            log.warning('No PDF can be built with dataset')
            return None

        l_arr_mass   = [ pdf.arr_mass for pdf in l_pdf ]
        l_arr_wgt    = [ pdf.arr_wgt  for pdf in l_pdf ]
        l_arr_sam    = [ pdf.arr_sam  for pdf in l_pdf ]
        l_arr_dec    = [ pdf.arr_dec  for pdf in l_pdf ]

        pdf.arr_mass = numpy.concatenate(l_arr_mass)
        pdf.arr_wgt  = numpy.concatenate(l_arr_wgt )
        pdf.arr_dec  = numpy.concatenate(l_arr_dec )
        pdf.arr_sam  = numpy.concatenate(l_arr_sam )

        return pdf
    #-----------------------------------------------------------
    @staticmethod
    def plot_pdf(
            pdf     : zpdf,
            name    : str,
            title   : str,
            out_dir : str,
            maxy    : float = None) -> None:
        '''
        Utility method, meant to plot PDF after it was built

        Parameters
        ------------------
        pdf    : PDF
        neme   : used to name the PNG file as {name}.png
        title  : Title for plots, will be appended after number of entries
        maxy   : Will be used to plot fit properly in case labels overlap
        out_dir: Directory where plots will go
        '''

        if pdf is None:
            log.warning(f'PDF {name} not build, not plotting')
            return

        arr_mass = pdf.arr_mass
        arr_wgt  = pdf.arr_wgt
        arr_sam  = pdf.arr_sam
        arr_dec  = pdf.arr_dec

        obj = ZFitPlotter(data=arr_mass, model=pdf, weights=arr_wgt)
        obj.plot(stacked=True)

        obj.axs[0].set_title(f'#Entries: {arr_mass.size}; {title}')

        if maxy is not None:
            obj.axs[0].set_ylim(bottom=0, top=maxy)

        obj.axs[0].axvline(x=5080, linestyle=':')
        obj.axs[0].axvline(x=5680, linestyle=':')
        obj.axs[0].axvline(x=5280, label=r'$B^+$', color='gray', linestyle='--')

        obj.axs[1].set_ylim(-5, +5)
        obj.axs[1].axhline(y=-3, color='red')
        obj.axs[1].axhline(y=+3, color='red')
        obj.axs[1].set_label('M$(B^+)$[MeV/${}_{c^2}$]')

        os.makedirs(out_dir, exist_ok=True)

        plot_path = f'{out_dir}/{name}.png'
        log.info(f'Saving to: {plot_path}')
        plt.savefig(plot_path)
        plt.close('all')

        plt.hist(arr_sam, bins=30, label='sample', histtype='step', linestyle='-' )
        plt.hist(arr_dec, bins=30, label='decay' , histtype='step', linestyle='--')
        plt.hist(arr_wgt, bins=30, label='Total' , histtype='step', linestyle=':' )

        plt.legend()
        plt.title(title)
        plt.savefig(f'{out_dir}/{name}_wgt.png')
        plt.close('all')

        text_path = plot_path.replace('png', 'txt')
        sut.print_pdf(pdf, txt_path=text_path)
    #-----------------------------------------------------------
    @staticmethod
    @contextmanager
    def apply_setting(use_cache : bool):
        '''
        Used to override default behaviour

        use_cache : If False (default is True) will recalculate the PDF
        '''
        old_val = PRec.use_cache
        try:
            PRec.use_cache = use_cache
            yield
        finally:
            PRec.use_cache = old_val
#-----------------------------------------------------------
