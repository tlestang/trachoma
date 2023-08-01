from ctypes import CDLL, POINTER, c_int
import numpy as np
from numpy.ctypeslib import ndpointer

diseased = np.array(
    [0, 0, 1, 1, 1, 1, 0, 0   , 0, 0, 0, 1, 1, 1, 1, 0]
).astype(np.bool_)
infected = np.array(
    [1, 0, 0, 0, 1, 0, 1, 0   , 1, 1, 0, 0, 1, 0, 1, 0]
).astype(np.bool_)
latent = np.array(
    [1, 0, 0, 0, 0, 0, 0, 1   , 1, 1, 1, 0, 0, 1, 0, 1]
).astype(np.bool_)
clock = np.array([2] * 16, dtype=np.int32);
clock[9] = 0

lib = CDLL('./libapply.so')
lib.apply_rules.restype = None

argtypes = [
    ndpointer(dtype=np.uint8, ndim=1),
    ndpointer(dtype=np.uint8, ndim=1),
    ndpointer(dtype=np.uint8, ndim=1),
    ndpointer(dtype=np.int32, ndim=1),
    c_int,
]
lib.apply_rules.argtypes = argtypes

inf = np.packbits(infected)
lat = np.packbits(latent)
dis = np.packbits(diseased)
lib.apply_rules(inf, dis, lat, clock, 16)

