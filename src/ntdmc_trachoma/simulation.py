import ctypes
from importlib import util
import json

import numpy as np

from .parameters import AverageDurations


LIBTRACHO_PATH = util.find_spec(
    "ntdmc_trachoma.libtrachoma"
).origin


class Simulation:

    def __init__(self, popsize, groups, ages, parameters_filepath):
        self.rng = np.random.default_rng()
        self.groups = groups
        self.lib = ctypes.CDLL(LIBTRACHO_PATH)

        with open(parameters_filepath, 'r') as f:
            p = json.load(f)
        self.set_base_periods(AverageDurations(p["durations"]))

    def set_base_periods(self, avg_durations):
        base_periods = (
            np.array([avg_durations.I] * self.popsize, dtype=np.int32),
            self.rng.poisson(
                lam=avg_durations.ID, size=self.popsize,
            ).astype(np.int32),
            self.rng.poisson(
                lam=avg_durations.D, size=self.popsize
            ).astype(np.int32),
        )

        self.lib.set_base_periods.argtypes = [
            np.ctypeslib.ndpointer(dtype=np.int32, ndim=1)
        ] * 3
        self.lib.set_base_periods(*base_periods)
        return None

    #TODO: Add getter for base periods, read arrays from C lib
