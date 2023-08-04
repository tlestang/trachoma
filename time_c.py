
import ctypes
import numpy as np
import setup

lib = ctypes.CDLL('./libtrachoma.so')
set_arrays = lib.set_arrays
set_arrays.restype = None
set_arrays.argtypes = (
    [np.ctypeslib.ndpointer(dtype=np.uint8, ndim=1)] * 4 +
    [np.ctypeslib.ndpointer(dtype=np.int32, ndim=1)] * 3 +
    [np.ctypeslib.ndpointer(dtype=np.float64, ndim=1)] * 2
)
set_times = lib.set_times
set_times.restype = None
set_times.argtypes = [np.ctypeslib.ndpointer(dtype=np.int32, ndim=1)] * 3

inf = np.packbits(setup.infected)
dis = np.packbits(setup.diseased)
lat = np.packbits(setup.latent)
new_i = np.zeros(len(lat), dtype=np.uint8)
set_arrays(
    inf, dis, lat, new_i, setup.clock, setup.ages,
    setup.infection_counter, setup.bact_load, setup.prob
)
set_times(setup.D_base, setup.ID_base, setup.latent_base)

lib.apply_rules.restype = None
lib.apply_rules.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.c_int]
lib.apply_rules(setup.n, 512, 1000)
