from ctypes import CDLL, c_double, c_int
from importlib import util
import numpy as np


sharedlib_path = util.find_spec("ntdmc_trachoma.libtrachoma")
lib = CDLL(sharedlib_path.origin)
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


def set_groups(p):
    groups = p.groups + [p.max_age]
    array_type = c_double * len(groups)
    lib.set_groups.restype = None
    lib.set_groups.argtypes = [array_type, c_int]
    lib.set_groups(array_type(*groups), len(groups))
