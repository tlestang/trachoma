import ctypes
import numpy as np
import numpy.testing as nptest

from state import Population


lib = ctypes.CDLL("./libtrachoma.so")


def test_remove_individual():
    ages = np.array(range(16))
    latent_base = np.array([2] * 16)
    rng = np.random.default_rng()
    pop = Population(ages, latent_base, rng)

    pop.inf = np.array(
        [
            1, 0, 1, 1, 1, 0, 1, 0,
            0, 1, 1, 1, 0, 0, 0, 1,
        ]
    )
    pop.dis = np.array(
        [
            1, 0, 0, 0, 1, 0, 1, 0,
            1, 0, 1, 0, 1, 0, 0, 1,
        ]
    )
    pop.lat = np.array(
        [
            1, 0, 0, 0, 0, 1, 1, 0,
            0, 1, 0, 1, 0, 1, 1, 0,
        ]
    )

    expected_inf = np.array(
        [
            0, 0, 1, 0, 1, 1, 0, 1,
            0, 0, 1, 1, 1, 0, 0, 1,
        ]
    ).astype(np.bool_)

    expected_dis = np.array(
        [
            0, 0, 1, 0, 0, 1, 0, 1,
            0, 1, 0, 1, 0, 0, 0, 1,
        ]
    ).astype(np.bool_)

    expected_lat = np.array(
        [
            0, 0, 1, 0, 0, 0, 1, 1,
            0, 0, 1, 0, 1, 1, 1, 0,
        ]
    ).astype(np.bool_)

    pop.clock = np.array(range(16), dtype=np.int32)
    expected_clock = [
        -1, -1, 0, 1, 3, 4, 5, 6,
        7, 8, 9, 10, 11, 13, 14, 15
    ]
    pop.count = np.array(range(16), dtype=np.int32)
    expected_count = [0, 0] + expected_clock[2:]
    pop.bact_load = np.array(range(16)) * 0.1
    expected_bact_load = np.array(
        [0, 0] + expected_clock[2:]
    ) * 0.1

    lib.remove_indiv.restype = None
    lib.remove_indiv(pop, 2)
    lib.remove_indiv(pop, 12)

    nptest.assert_array_equal(expected_inf, pop.inf)
    nptest.assert_array_equal(expected_dis, pop.dis)
    nptest.assert_array_equal(expected_lat, pop.lat)
    nptest.assert_array_equal(expected_clock, pop.clock)
    nptest.assert_array_equal(expected_count, pop.count)
    nptest.assert_array_equal(expected_bact_load, pop.bact_load)
