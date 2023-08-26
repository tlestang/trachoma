import ctypes
from ctypes import c_double, c_int
from importlib import util
import numpy as np
from numpy.ctypeslib import ndpointer

LIBTRACHO_PATH = util.find_spec(
    "ntdmc_trachoma.libtrachoma"
).origin
lib = ctypes.CDLL(LIBTRACHO_PATH)


def test_infection_prob():
    # COMPUTE INFECTION PROBABLITY USING C CORE ####
    PHI = 1
    EPSILON = 0.5
    V1 = 1.0
    V2 = 3
    lib.set_infection_parameters.restype = None
    lib.set_infection_parameters.argtypes = [c_double] * 4
    lib.set_infection_parameters(V1, V2, PHI, EPSILON)

    groups = (c_int * 3)(*[468, 780, 3121])
    lib.set_groups.restype = None
    lib.set_groups.argtypes = [c_int * 3]
    lib.set_groups(groups, 3)

    ages = np.array(
        [5 * 52] * 20 + [12 * 52] * 50 + [45 * 52] * 30,
        dtype=np.int32
    )
    rng = np.random.default_rng()
    bact_load = np.concatenate((
        rng.normal(loc=0.1, scale=0.1, size=20),
        rng.normal(loc=0.2, scale=0.1, size=50),
        rng.normal(loc=0.3, scale=0.1, size=30),
        )
    )
    BETA = 1.
    n = 100

    lib.get_infection_prob.restype = None
    lib.get_infection_prob.argtypes = [
        ndpointer(dtype=np.int32, ndim=1),
        ndpointer(dtype=np.float64, ndim=1),
        c_int,
        c_double,
        ndpointer(dtype=np.float64, ndim=1),
    ]
    prob = np.zeros(len(ages), dtype=np.float64)
    lib.get_infection_prob(ages, bact_load, n, BETA, prob)

    # COMPUTE EXPECTED INFECTION PROBABILITY ####
    prevLambda = [
        BETA * (
            V1 * bact_load[0:20].mean() + V2 * (bact_load[0:20].mean()**(PHI + 1))
        ),
        BETA * (
            V1 * bact_load[20:70].mean() + V2 * (bact_load[20:70].mean()**(PHI + 1))
        ),
        BETA * (
            V1 * bact_load[70:100].mean() + V2 * (bact_load[70:100].mean()**(PHI + 1))
        ),
    ]
    demog_matrix = np.array([0.2, 0.5, 0.3])
    social_mixing = (
        (EPSILON * np.diag(np.ones(3)) + (1 - EPSILON)) * demog_matrix
    )
    P = 1. - np.exp(- np.dot(social_mixing, prevLambda))
    expected_prob = np.array(
        [P[0]] * 20 + [P[1]] * 50 + [P[2]] * 30
    )

    np.testing.assert_allclose(prob, expected_prob)
