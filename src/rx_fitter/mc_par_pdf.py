'''
Module with class MCParPdf
'''
# pylint: disable=too-many-positional-arguments, too-many-function-args, too-many-arguments, too-many-locals, too-many-instance-attributes

import os
import json

from ROOT                                        import RDataFrame
from dmu.logging.log_store                       import LogStore
from dmu.stats.model_factory                     import ModelFactory
from dmu.generic                                 import version_management as vman
from zfit.core.basepdf                           import ZfitPDF            as zpdf
from zfit.core.interfaces                        import ZfitSpace          as zobs
from rx_calibration.hltcalibration.fit_component import FitComponent, load_fit_component

log = LogStore.add_logger('rx_fitter:mc_par_pdf')
# ---------------------------------------
class MCParPdf:
    '''
    Class intended to provide FitComponent instances from fits to MC
    '''
    # ---------------------------------------
    def __init__(
            self,
            rdf    : RDataFrame,
            obs    : zobs,
            cfg    : dict) -> FitComponent:

        self._rdf    = rdf
        self._obs    = obs
        self._cfg    = cfg
        self._mass   = obs.obs[0]
        self._fit_dir= os.environ['FITDIR']

        self._sample = cfg['name'   ]
        self._q2bin  = cfg['q2bin'  ]
        self._trigger= cfg['trigger']
        self._nbrem  = cfg['nbrem'  ]
        self._fvers  = cfg['fvers'  ]

        self._shared = ['mu'] if cfg['shared'] is None else cfg['shared']
        self._model  = self._get_model(cfg['model'])

        if self._sample == 'Bu_JpsiK_ee_eq_DPC':
            self._pfloat = ['mu', 'sg'] if cfg['pfloat'] is None else cfg['pfloat']
        else:
            self._pfloat = []

        self._cfg['out_dir'] = self._get_pars_dir()
    # ---------------------------------------
    def _get_pars_dir(self, version : str = None) -> str:
        model_name = '_'.join(self._model)
        init_dir = f'{self._fit_dir}/mc/{self._q2bin}'
        fnal_dir = f'{self._sample}_{self._trigger}/{self._mass}_{self._nbrem}/{model_name}'

        if version is not None:
            log.debug(f'Using user defined version of fit in: {init_dir}')
            return f'{init_dir}/{version}/{fnal_dir}'

        if not os.path.isdir(init_dir):
            init_dir = f'{init_dir}/v1'
            log.info(f'No fitting path found, making first version of fit directory in: {init_dir}')
            return f'{init_dir}/{fnal_dir}'

        init_dir = vman.get_last_version(dir_path=init_dir, version_only=False)
        if self._cfg['create']:
            init_dir = vman.get_next_version(init_dir)
            log.info(f'Creating new version of fit in: {init_dir}')
        else:
            log.info(f'Using latest version of fit in: {init_dir}')

        return f'{init_dir}/{fnal_dir}'
    # ---------------------------------------
    def _get_model(self, model : list[str]) -> list[str]:
        if model is not None:
            return model

        log.info('Model not passed, will pick default')

        is_sig = self._sample  == 'Bu_JpsiK_ee_eq_DPC'
        is_trg = self._trigger == 'Hlt2RD_BuToKpEE_MVA'
        is_jps = self._q2bin   == 'jpsi'
        is_brm = self._nbrem   in [0, 1, 2]

        if is_sig and is_jps and is_brm and is_trg:
            return {
                    0 : ['suj', 'suj' ],
                    1 : ['suj', 'dscb'],
                    2 : ['suj', 'dscb']}[self._nbrem]

        if self._sample == 'Bu_JpsiPi_ee_eq_DPC':
            return {
                    0 : ['suj'],
                    1 : ['suj'],
                    2 : ['suj']}[self._nbrem]

        raise ValueError(f'Cannot assign default model for: {self._sample}/{self._q2bin}/{self._trigger}/{self._nbrem}')
    # ------------------------------------
    def _fix_tails(self, pdf : zpdf, fix_dir : str) -> zpdf:
        json_path = f'{fix_dir}/fit.json'
        if not os.path.isfile(json_path):
            log.debug('No tail fixing path found, leaving all parameters to float')
            return pdf

        log.info(40 * '-')
        log.info(f'Fixing parameters with: {json_path}')
        log.info(40 * '-')
        s_par = pdf.get_params()

        with open(json_path, encoding='utf-8') as ifile:
            d_par = json.load(ifile)

        for par in s_par:
            if par.name not in d_par:
                continue

            [val, _] = d_par[par.name]

            par.set_value(val)

            log.info(f'{par.name:<30}{"--->":<10}{val:.3f}')
            par.floating = False

        return zpdf
    # ---------------------------------------
    def get_fcomp(self) -> FitComponent:
        '''
        Returns instance of FitComponent
        '''
        log.debug(f'Bulding model: {self._model}')
        mod   = ModelFactory(preffix=self._sample, obs=self._obs, l_pdf=self._model, l_shared=self._shared, l_float=self._pfloat)
        pdf   = mod.get_pdf()

        obj   = load_fit_component(cfg=self._cfg, pdf=pdf)
        if obj is not None:
            log.info('Will load PDF from cached parameters file')
            return obj

        fix_dir = self._get_pars_dir(self._fvers)
        pdf     = self._fix_tails(pdf=pdf, fix_dir=fix_dir)

        obj     = FitComponent(cfg=self._cfg, rdf=self._rdf, pdf=pdf, obs=self._obs)

        return obj
# ---------------------------------------
