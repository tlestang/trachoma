import ctypes
from ctypes import c_int, c_ubyte
import numpy as np
import numpy.testing as nptest

from state import Population


lib = ctypes.CDLL("./libtrachoma.so")

# We'll be using the set<period>time functions provided by the C
# library (periods.c). So we need to initialise the base period
# arrays before we make the calls.
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

RNG = np.random.default_rng()


def mk_bitarray(indices):
    """Return boolean array with True at locations given by 'indices'
    and False elsewhere.
    """
    a = np.zeros(16, dtype=np.bool_)
    a[list(indices)] = True
    return a


def set_expected_clock(
        pop, transition, latent, infected, new_infections
):
    # Decrement everyone's clock.  Individuals transitioning to
    # another stage will have their value overriden anyway.
    expected_clock = pop.clock - 1

    # Finally, it's time to set the clock value for individual
    # transitioning to the next stage.
    for i in transition & latent:
        expected_clock[i] = lib.setidtime(
            2, c_int(pop.count[i]), c_int(pop.ages[i])
        )
    for i in transition & (infected - latent):
        expected_clock[i] = lib.setdtime(
            2, c_int(pop.count[i]), c_int(pop.ages[i])
        )
    for i in new_infections:
        expected_clock[i] = lib.setlatenttime(
            2, c_int(pop.count[i]), c_int(pop.ages[i])
        )

    return expected_clock


def test_spread():
    ages = np.array(range(16), dtype=np.int32)
    pop = Population(ages, latent_base_period, RNG)

    latent = {0, 4, 9, 14}
    infected = {0, 2, 4, 9, 10, 14}
    diseased = {2, 5, 6, 10, 15}
    new_infections = {1, 6, 8, 12}
    pop.lat = mk_bitarray(latent)
    pop.inf = mk_bitarray(infected)
    pop.dis = mk_bitarray(diseased)

    new_i = mk_bitarray(new_infections)

    transition = {4, 14, 2, 6}
    pop.clock = (
        RNG.integers(low=1, high=11, size=16, dtype=np.int32) *
        (1 - mk_bitarray(transition).astype(np.int32))
    )

    pop.bact_load = np.zeros(16, dtype=np.float64)
    pop.bact_load[pop.inf] = 1.
    pop.count = np.zeros(16, dtype=np.int32)
    pop.count += pop.lat + pop.inf + pop.dis

    expected_lat = mk_bitarray([0, 9, 6, 1, 8, 12])
    expected_inf = mk_bitarray([0, 9, 6, 1, 8, 12, 14, 4, 10])
    expected_dis = mk_bitarray([14, 4, 10, 2, 5, 15])

    expected_clock = set_expected_clock(
        pop, transition, latent, infected, new_infections
    )
    expected_count = pop.count + new_i

    lib.spread.restype = None
    lib.spread(pop, c_ubyte(np.packbits(new_i)[0]), c_int(0))
    lib.spread(pop, c_ubyte(np.packbits(new_i)[1]), c_int(1))

    nptest.assert_array_equal(expected_lat, pop.lat)
    nptest.assert_array_equal(expected_inf, pop.inf)
    nptest.assert_array_equal(expected_dis, pop.dis)
    nptest.assert_array_equal(expected_count, pop.count)
    nptest.assert_array_equal(expected_clock, pop.clock)
