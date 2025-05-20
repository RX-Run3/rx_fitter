'''
Module containing class ModelScales and helper functions
'''

import os
import math
from importlib.resources import files

import yaml
import numpy
import pandas                   as pnd
import jacobi                   as jac

from dmu.logging.log_store                 import LogStore
from dmu.generic                           import hashing
from dmu.generic.version_management        import get_last_version
from rx_selection                          import selection  as sel
from rx_efficiencies.decay_names           import DecayNames as dn
from rx_efficiencies.efficiency_calculator import EfficiencyCalculator

log=LogStore.add_logger('rx_fitter:prec_scales')
#------------------------------------------
class PrecScales:
    '''
    Class used to calculate scale factor between yields of partially reconstructed component and signal
    '''
    #------------------------------------------
    def __init__(self, proc : str, q2bin : str):
        '''
        proc : Nickname of decay process, nicknames are in the DecayNames class
        q2bin: Needed to apply correct selection to get correct efficiencies and scales
        '''
        self._proc        = proc
        self._q2bin       = q2bin

        self._d_frbf      : dict
        self._df_eff      : pnd.DataFrame
        self._trigger     = 'Hlt2RD_BuToKpEE_MVA'
        self._initialized = False

        self._hash        = self._get_hash()
    #------------------------------------------
    def _get_hash(self) -> str:
        project = {'Hlt2RD_BuToKpEE_MVA' : 'RK'}[self._trigger]
        process = dn.sample_from_decay(self._proc)

        d_sel   = sel.selection(trigger=self._trigger, q2bin=self._q2bin, process=process)
        hsh     = hashing.hash_object([self._proc, self._q2bin, d_sel])

        return hsh
    #------------------------------------------
    def _check_arg(self, l_val, val, name):
        if val not in l_val:
            raise ValueError(f'{name} {val} not allowed')

        log.debug(f'{name:<20}{"->":20}{val:<20}')
    #------------------------------------------
    def _initialize(self):
        if self._initialized:
            return

        log.debug('Initializing')
        self._load_fractions()
        self._load_efficiencies()

        self._initialized = True
    #------------------------------------------
    def _load_fractions(self):
        log.debug('Getting hadronization fractions and branching ratios')

        frbf_dir  = files('rx_efficiencies_data').joinpath('prec_sf')
        frbf_path = get_last_version(dir_path=frbf_dir, version_only=False)
        frbf_path = f'{frbf_path}/fr_bf.yaml'

        log.debug(f'Picking up branching fractions from: {frbf_path}')
        with open(frbf_path, encoding='utf-8') as ifile:
            self._d_frbf = yaml.safe_load(ifile)
    #------------------------------------------
    def _calculate_efficiencies(self, yaml_path : str) -> None:
        log.debug('Efficiencies not found, calculating them')
        out_dir     = os.path.dirname(yaml_path)
        obj         = EfficiencyCalculator(q2bin=self._q2bin)
        obj.out_dir = out_dir
        df          = obj.get_stats()
        d_data      = df.to_dict()
        with open(yaml_path, 'w', encoding='utf-8') as ofile:
            yaml.safe_dump(d_data, ofile)
    #------------------------------------------
    def _load_efficiencies(self):
        log.debug('Getting efficiencies')

        eff_dir  = files('rx_efficiencies_data').joinpath('prec_sf')
        eff_path = get_last_version(dir_path=eff_dir, version_only=False)
        eff_path = f'{eff_path}/efficiencies_{self._q2bin}/{self._hash}/data.yaml'

        if not os.path.isfile(eff_path):
            self._calculate_efficiencies(yaml_path=eff_path)

        with open(eff_path, encoding='utf-8') as ifile:
            data = yaml.safe_load(ifile)

        df = pnd.DataFrame(data)

        self._df_eff = df
    #------------------------------------------
    def _get_fr(self, proc : str) -> float:
        '''
        Returns hadronization fraction for given process
        '''
        if   proc.startswith('bp'):
            fx = 'fu'
        elif proc.startswith('bd'):
            fx = 'fd'
        elif proc.startswith('bs'):
            fx = 'fs'
        else:
            raise ValueError(f'Cannot find hadronization fraction for: {proc}')

        fx = self._d_frbf['fr'][fx]

        return fx
    #------------------------------------------
    def _mult_brs(self, l_br : list[tuple[float,float]]):
        log.debug('Multiplying branching fractions')

        l_br_val = [ float(br[0]) for br in l_br ] # These numbers come from YAML files
        l_br_err = [ float(br[1]) for br in l_br ] # when using "e" in scientific notation, these numbers are made into strings
        br_cov   = numpy.diag(l_br_err) ** 2
        val, var = jac.propagate(math.prod, l_br_val, br_cov)
        err      = math.sqrt(var)

        return val, err
    #------------------------------------------
    def _get_br(self, proc : str) -> tuple[float,float]:
        log.debug(f'Calculating BR for {proc}')

        l_dec = dn.subdecays_from_decay(proc)
        l_bf  = [ self._d_frbf['bf'][dec] for dec in l_dec ]

        return self._mult_brs(l_bf)
    #------------------------------------------
    def _get_ef(self, proc : str):
        log.debug(f'Calculating efficiencies for {proc}')

        df = self._df_eff
        df = df[df.Process == proc]

        if len(df) != 1:
            print(df)
            raise ValueError('Expected one and only one row for process {proc}')

        passed, total = df.iloc[-1][['Passed', 'Total']]

        eff = passed / total
        err = math.sqrt(eff * (1 - eff) / total)

        return eff, err
    #------------------------------------------
    def _print_vars(self, l_tup : list[tuple[float,float]], proc : str) -> None:
        log.debug('')
        log.debug(f'Decay: {proc}')
        log.debug('-' * 20)
        log.debug(f'{"Var":<20}{"Value":<20}{"Error":<20}')
        log.debug('-' * 20)
        for (val, err), name in zip(l_tup, ['fr', 'br', 'eff']):
            log.debug(f'{name:<20}{val:<20.3e}{err:<20.3e}')
        log.debug('-' * 20)
    #------------------------------------------
    def get_scale(self, signal : str) -> tuple[float,float]:
        '''
        Returns scale factor k and error, meant to be used in:

        Nprec = k * Nsignal

        reparametrization, during fit.

        Parameters
        -----------------------
        signal: String representing signal WRT which to scale, e.g. bpkpee
        '''
        self._initialize()

        fr_bkg = self._get_fr(self._proc)
        br_bkg = self._get_br(self._proc)
        ef_bkg = self._get_ef(self._proc)

        fr_sig = self._get_fr(signal)
        br_sig = self._get_br(signal)
        ef_sig = self._get_ef(signal)

        l_tup = [fr_sig, br_sig, ef_sig, fr_bkg, br_bkg, ef_bkg]
        l_val = [ tup[0] for tup in l_tup]
        l_err = [ tup[1] for tup in l_tup]
        cov   = numpy.diag(l_err) ** 2

        self._print_vars(l_tup[:3], proc=    signal)
        self._print_vars(l_tup[3:], proc=self._proc)

        val, var = jac.propagate(lambda x : (x[3] * x[4] * x[5]) / (x[0] * x[1] * x[2]), l_val, cov)
        val = float(val)
        err = math.sqrt(var)

        return val, err
#------------------------------------------
