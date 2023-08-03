import numpy as np
from numpy.ctypeslib import ndpointer

from ctypes import CDLL, c_int
lib = CDLL('./libtrachoma.so')
lib.apply_rules.restype = None
argtypes = [
    ndpointer(dtype=np.uint8, ndim=1),
    ndpointer(dtype=np.uint8, ndim=1),
    ndpointer(dtype=np.uint8, ndim=1),
    ndpointer(dtype=np.int32, ndim=1),
    ndpointer(dtype=np.int32, ndim=1),
    ndpointer(dtype=np.float64, ndim=1),
    c_int,
]
lib.apply_rules.argtypes = argtypes

ages = np.sort(
    np.array([
        1824, 388, 2032, 1446, 1323, 2222, 2889, 846, 586,
        2840, 2650, 2850, 203, 1400, 1917, 1327,
    ], dtype=np.int32)
)
bactld = np.array(
    [
        0.54589385, 0.64449697, 0.94586293, 0.89743933, 0.74142636,
        0.42041218, 0.30355981, 0.76096054, 0.62761568, 0.88280128,
        0.16847746, 0.41943052, 0.56054159, 0.09329043, 0.67442931,
        0.16027857,
    ], dtype=np.float64
)
diseased = np.array(
    [0, 0, 1, 1, 1, 1, 1, 0   , 0, 0, 0, 1, 1, 0, 1, 0]
).astype(np.bool_)
infected = np.array(
    [1, 0, 0, 0, 1, 0, 1, 1   , 1, 1, 0, 0, 1, 0, 1, 0]
).astype(np.bool_)
latent = np.array(
    [1, 0, 0, 0, 0, 0, 0, 1   , 1, 1, 1, 0, 0, 1, 0, 1]
).astype(np.bool_)
clock = np.array([2] * 16, dtype=np.int32);
clock[9] = 0

inf = np.packbits(infected)
lat = np.packbits(latent)
dis = np.packbits(diseased)
lib.apply_rules(inf, dis, lat, clock, ages, bactld, 16)

