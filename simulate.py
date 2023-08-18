from ctypes import CDLL, c_double, c_int
import numpy as np
from state import Population

lib = CDLL('./libtrachoma.so')
advance = lib.apply_rules
advance.restype = None

lib.set_base_periods.argtypes = [
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1)
] * 3


def simulate(task: list, events: list):
    for p, beta in task:
        for e, e_next in zip(events[:-1], events[1:]):
            advance(p, e_next[0] - e[0], c_double(beta))
            e_next[1](p)


def set_base_periods(p, rng):
    base_periods = (
        np.array([p.av_I_duration] * p.population_size, dtype=np.int32),
        rng.poisson(
            lam=p.av_ID_duration, size=p.population_size
        ).astype(np.int32),
        rng.poisson(
            lam=p.av_D_duration, size=p.population_size
        ).astype(np.int32),
    )
    lib.set_base_periods(*base_periods)
    return base_periods


def set_infection_parameters(p):
    lib.set_infection_parameters.restype = None
    lib.set_infection_parameters.argtypes = [c_double] * 4
    lib.set_infection_parameters(
        c_double(p.v1),
        c_double(p.v2),
        c_double(p.phi),
        c_double(p.epsilon),
    )


def set_background_mortality(p):
    rate = 1. - np.exp(- 1. / p.tau)
    lib.set_background_mortality.restype = None
    lib.set_background_mortality.argtypes = [c_double]
    lib.set_background_mortality(
        c_double(rate)
    )


def set_groups(p):
    groups = p.groups + [p.max_age]
    array_type = c_double * len(groups)
    lib.set_groups.restype = None
    lib.set_groups.argtypes = [array_type, c_int]
    lib.set_groups(array_type(*groups), len(groups))
