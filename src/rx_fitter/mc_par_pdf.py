'''
Module with class MCParPdf
'''
# pylint: disable=too-many-positional-arguments, too-many-function-args, too-many-arguments, too-many-locals, too-many-instance-attributes

import os
import copy

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
    Class intended to provide zfit PDF instances from fits to MC, thin wrapper around FitComponent class
    '''
    # ---------------------------------------
    def __init__(
            self,
            rdf    : RDataFrame,
            obs    : zobs,
            cfg    : dict) -> FitComponent:

        self._rdf    = rdf
        self._obs    = obs
        self._cfg    = copy.deepcopy(cfg)
        self._mass   = obs.obs[0]

        self._component_name = self._cfg['component_name']

        self._q2bin       = self._cfg['q2bin'  ]
        self._trigger     = self._cfg['trigger']
        self._nbrem       = self._cfg['nbrem'  ]
        self._model       = self._cfg['model'  ]

        self._out_dir        = cfg['output' ]['out_dir']
        self._cfg['out_dir'] = self._get_pars_dir()
        self._cfg['name'   ] = self._component_name
    # ---------------------------------------
    def _get_pars_dir(self, version : str = None) -> str:
        model_name = '_'.join(self._model)
        init_dir   = f'{self._out_dir}/{self._q2bin}'
        fnal_dir   = f'{self._component_name}_{self._trigger}/{self._mass}_{self._nbrem}/{model_name}'

        if version is not None:
            pars_dir = f'{init_dir}/{version}/{fnal_dir}'
            log.debug(f'Using user defined version of fit in: {pars_dir}')
            return pars_dir

        if not os.path.isdir(init_dir):
            init_dir = f'{init_dir}/v1'
            log.info(f'No output directory found, making first version of fit directory in: {init_dir}')
            return f'{init_dir}/{fnal_dir}'

        log.debug(f'Looking for latest version in: {init_dir}')
        init_dir = vman.get_last_version(dir_path=init_dir, version_only=False)
        log.info(f'Will fit and save to: {init_dir}')

        return f'{init_dir}/{fnal_dir}'
    # ------------------------------------
    def get_pdf(self, must_load_pars : bool = False) -> zpdf:
        '''
        Returns instance of zfit PDF

        must_load_pars (bool): Will use must_load_pars of fit component. If RDF is missing and if this flag is true, will raise NoFitDataFoundException
        '''
        log.debug(f'Bulding model: {self._model}')
        d_rep = None
        if 'reparametrize' in self._cfg:
            d_rep = self._cfg['reparametrize']
            log.info(f'Reparametrizing PDF: {d_rep}')
        else:
            log.debug('Not reparametrizing PDF')

        preffix = f'{self._component_name}_{self._nbrem:03}'

        mod   = ModelFactory(
                obs     = self._obs,
                preffix = preffix,
                l_pdf   = self._cfg['model' ],
                d_rep   = d_rep,
                l_shared= self._cfg['shared'],
                l_float = self._cfg['pfloat'])

        pdf   = mod.get_pdf()

        obj   = load_fit_component(cfg=self._cfg, pdf=pdf)
        if obj is not None:
            log.info('Will load PDF from cached parameters file')
            return pdf

        log.debug('No fixing version found, using original PDF for fit component object')
        obj = FitComponent(cfg=self._cfg, rdf=self._rdf, pdf=pdf, obs=self._obs)
        pdf = obj.get_pdf(must_load_pars)

        return pdf
# ---------------------------------------
