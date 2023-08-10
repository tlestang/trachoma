from ctypes import CDLL
from state import Population

lib = CDLL('./libtrachoma.so')
advance = lib.apply_rules
advance.restype = None


def simulate(p: Population, events: list):
    for e, e_next in zip(events[:-1], events[1:]):
        nsteps = int((e_next[1] - e[1]).days / 7)
        advance(p, nsteps)
        e_next(p)
