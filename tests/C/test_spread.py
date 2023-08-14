import ctypes
from ctypes import c_int
import numpy as np
import numpy.testing as nptest

from state import Population


lib = ctypes.CDLL("./libtrachoma.so")


def test_spread():
    ages = np.array(range(16), dtype=np.int32)
    latent_base = np.array([2] * 16)
    rng = np.random.default_rng()
    pop = Population(ages, latent_base, rng)

    # lat+inf: 0, 4, 9, 14
    # inf+dis: 2, 10
    # dis: 5, 6, 15
    # sus: 1, 3, 7, 8, 11, 12, 13
    # clock: 4, 14, 2, 6
    # new_i: 1, 8, 12, 6
    pop.lat = np.array(
        [
            1, 0, 0, 0, 1, 0, 0, 0,
            0, 1, 0, 0, 0, 0, 1, 0,
        ]
    )
    pop.inf = np.array(
        [
            1, 0, 1, 0, 1, 0, 0, 0,
            0, 1, 1, 0, 0, 0, 1, 0,
        ]
    )
    pop.dis = np.array(
        [
            0, 0, 1, 0, 0, 1, 1, 0,
            0, 0, 1, 0, 0, 0, 0, 1,
        ]
    )
    pop.clock = np.array(
        [
            1, 1, 0, 1, 0, 1, 0, 1,
            1, 1, 1, 1, 1, 1, 0, 1,
        ]
    ).astype(np.int32)
    new_i = np.array(
        [
            0, 1, 0, 0, 0, 0, 1, 0,
            1, 0, 0, 0, 1, 0, 0, 0,
        ]
    )
    pop.bact_load = np.zeros(16, dtype=np.float64)
    pop.bact_load[pop.inf] = 1.
    pop.count = np.zeros(16, dtype=np.int32)
    pop.count += pop.lat + pop.inf + pop.dis

    # expected lat+inf: 0, 9, 6, 1, 8, 12
    # expected inf+dis: 14, 4, 10
    # expected dis: 2, 5, 15
    # expected sus: 3, 7, 11, 13
    expected_lat = np.array(
        [
            1, 1, 0, 0, 0, 0, 1, 0,
            1, 1, 0, 0, 1, 0, 0, 0,
        ]
    )
    expected_inf = np.array(
        [
            1, 1, 0, 0, 1, 0, 1, 0,
            1, 1, 1, 0, 1, 0, 1, 0,
        ]
    )
    expected_dis = np.array(
        [
            0, 0, 1, 0, 1, 1, 0, 0,
            0, 0, 1, 0, 0, 0, 1, 1,
        ]
    )
    expected_count = pop.count + new_i
    expected_clock = pop.clock - 1
    latent_base_period = np.array([2] * 16, dtype=np.int32)
    ID_base_period = latent_base_period.copy()
    D_base_period = latent_base_period.copy()
    lib.set_base_periods.restype = None
    lib.set_base_periods.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.int32, ndim=1)
    ] * 3
    lib.set_base_periods(
        latent_base_period, ID_base_period, D_base_period
    )
    for i in [4, 14]:
        expected_clock[i] = lib.setdtime(
            2, c_int(pop.count[i]), c_int(ages[i])
        )
    expected_clock[2] = lib.setidtime(2, c_int(pop.count[2]), c_int(ages[2]))
    for i in [1, 6, 8, 12]:
        expected_clock[i] = lib.setlatenttime(
            2, c_int(pop.count[i]), c_int(ages[i])
        )

    lib.spread.restype = None
    lib.spread(pop, ctypes.c_ubyte(np.packbits(new_i)[0]), ctypes.c_int(0))
    lib.spread(pop, ctypes.c_ubyte(np.packbits(new_i)[1]), ctypes.c_int(1))

    nptest.assert_array_equal(expected_lat, pop.lat)
    nptest.assert_array_equal(expected_inf, pop.inf)
    nptest.assert_array_equal(expected_dis, pop.dis)
    nptest.assert_array_equal(expected_count, pop.count)
    nptest.assert_array_equal(expected_clock, pop.clock)
