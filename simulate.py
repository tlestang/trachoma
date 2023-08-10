from ctypes import CDLL
import numpy as np
from state import Population

lib = CDLL('./libtrachoma.so')
advance = lib.apply_rules
advance.restype = None

lib.set_base_periods.argtypes = [
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1)
] * 3


def simulate(p: Population, events: list):
    for e, e_next in zip(events[:-1], events[1:]):
        nsteps = int((e_next[1] - e[1]).days / 7)
        advance(p, nsteps)
        e_next(p)


def set_base_periods(p, rng):
    base_periods = (
        np.array([p.av_I_duration] * p.popsize, dtype=np.float64),
        rng.poisson(lam=p.av_ID_duration, size=p.popsize),
        rng.poisson(lam=p.av_D_duration, size=p.popsize),
    )
    lib.set_base_periods(*base_periods)
    return base_periods
