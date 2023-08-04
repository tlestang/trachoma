import numpy as np
import init
import periods
import infection
from parameters import AV_D_DURATION, AV_ID_DURATION

import ctypes

n = 2048
rng = np.random.default_rng()
ages = np.sort(init.ages(n, 60 * 52, 20 * 52, rng)).astype(np.int32)

clock = np.zeros(n, dtype=np.int32) - 1
bact_load = np.zeros(n, dtype=np.float64)
prob = np.zeros(n, dtype=np.float64)
infection_counter = np.zeros(n, dtype=np.int32)
latent = init.infected(n, 0.01, rng)
# latent = np.array([False, True, False, False, False, False, False, False, False,
#                   False, False, True, False, True, False, False])
infected = latent.copy()
diseased = np.zeros(n, dtype=np.bool_)

latent_base = np.zeros(n, dtype=np.int32) + 1
ID_base = rng.poisson(lam=AV_ID_DURATION, size=n).astype(np.int32)
D_base = rng.poisson(lam=AV_D_DURATION, size=n).astype(np.int32)

clock[latent] = periods.latent_time(latent_base[latent])
infection_counter[latent] = 1
bact_load[latent] = infection.get_load(infection_counter[latent])


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

inf = np.packbits(infected)
dis = np.packbits(diseased)
lat = np.packbits(latent)
new_i = np.zeros(len(lat), dtype=np.uint8)
set_arrays(
    inf, dis, lat, new_i, clock, ages,
    infection_counter, bact_load, prob
)
set_times(D_base, ID_base, latent_base)

lib.apply_rules.restype = None
lib.apply_rules.argtypes = [ctypes.c_int, ctypes.c_int]
lib.apply_rules(n, 256)
