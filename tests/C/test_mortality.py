import ctypes
from importlib import util
import numpy as np
import numpy.testing as nptest

from ntdmc_trachoma.state import Population


LIBTRACHO_PATH = util.find_spec(
    "ntdmc_trachoma.libtrachoma"
).origin
lib = ctypes.CDLL(LIBTRACHO_PATH)


def test_remove_individual():
    ages = np.array(range(16))
    latent_base = np.array([2] * 16)
    rng = np.random.default_rng()
    pop = Population(ages, latent_base, rng, lib)

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


def test_old_age_mortality():
    ages = np.array(range(24))
    latent_base = np.array([2] * 24)
    rng = np.random.default_rng()
    pop = Population(ages, latent_base, rng, lib)

    pop.inf = np.array(
        [
            1, 0, 1, 1, 1, 0, 1, 0,
            0, 1, 1, 1, 0, 0, 0, 1,
            1, 1, 1, 1, 1, 1, 1, 1
        ]
    )

    expected_ages = np.array([0] * 11 + list(range(13)))
    expected_inf = np.array(
        [
            0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 1, 0, 1, 1, 1,
            0, 1, 0, 0, 1, 1, 1, 0,
        ]
    )
    lib.old_age_mortality.restype = None
    lib.old_age_mortality(pop, 11)

    nptest.assert_array_equal(expected_inf, pop.inf)
    nptest.assert_array_equal(expected_ages, pop.ages)
