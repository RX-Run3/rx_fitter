'''
Module with functions intended to interface with the PDG API
'''

import pdg

from dmu.logging.log_store import LogStore

log=LogStore.add_logger('scripts:pdg_utils')
#-------------------------------------------------------
def get_bf(decay : str) -> float:
    '''
    Returns branching fraction for a given decay
    '''
    api    = pdg.connect()
    mother = decay.split('-->')[0].replace(' ', '')
    for bf in api.get_particle_by_name(mother).exclusive_branching_fractions():
        if bf.is_limit:
            continue

        if bf.description == decay:
            return bf.value

    log.error(f'Cannot find BF for decay: {decay}')
    raise ValueError('Make sure your pdg>=0.1.2')
#-------------------------------------------------------
