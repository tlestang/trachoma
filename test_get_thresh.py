from ctypes import CDLL, POINTER, c_int
lib = CDLL('./libtrachoma.so')

import numpy as np
from numpy.ctypeslib import ndpointer
import infection

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
prob = np.zeros(16, dtype=np.float64)

probpy = infection.getlambdaStep(ages, bactld)

lib.get_infection_prob.restype = None

argtypes = [
    ndpointer(dtype=np.int32, ndim=1),
    ndpointer(dtype=np.float64, ndim=1),
    ndpointer(dtype=np.float64, ndim=1),
    c_int,
]
lib.get_infection_prob.argtypes = argtypes

lib.get_infection_prob(ages, bactld, prob, 16)
