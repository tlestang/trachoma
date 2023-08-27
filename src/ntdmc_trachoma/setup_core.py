import ctypes
from numpy import int32, exp
from numpy.ctypeslib import ndpointer

from .parameters import InfectionParameters, BasePeriods


def set_base_periods(lib, base_periods: BasePeriods):
    lib.set_base_periods.argtypes = [
        ndpointer(dtype=int32, ndim=1)
    ] * 3
    lib.set_base_periods(*base_periods)
    return base_periods


def set_infection_parameters(lib, p: InfectionParameters):
    lib.set_infection_parameters.restype = None
    lib.set_infection_parameters.argtypes = [ctypes.c_double] * 4
    lib.set_infection_parameters(
        ctypes.c_double(p.v1),
        ctypes.c_double(p.v2),
        ctypes.c_double(p.phi),
        ctypes.c_double(p.epsilon),
    )


def set_bgd_mortality(lib, tau: float):
    lib.set_background_mortality.restype = None
    lib.set_background_mortality.argtypes = [ctypes.c_double]
    lib.set_background_mortality(
        ctypes.c_double(1. - exp(- 1. / tau))
    )


def set_groups(lib, groups: list[int]):
    lib.set_groups.restype = None
    array_type = ctypes.c_int * len(groups)
    lib.set_groups.argtypes = [array_type, ctypes.c_int]
    lib.set_groups(array_type(*groups), len(groups))
