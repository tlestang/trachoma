import ctypes
from importlib import util
import json
from pathlib import Path

import numpy as np

from .parameters import AverageDurations, InfectionParameters


LIBTRACHO_PATH = util.find_spec(
    "ntdmc_trachoma.libtrachoma"
).origin


class Simulation:

    def __init__(self, parameters_filepath: Path):
        self.rng = np.random.default_rng()
        self.lib = ctypes.CDLL(LIBTRACHO_PATH)
        self.load_parameters(parameters_filepath)

    def load_parameters(self, parameters_filepath: Path):
        with parameters_filepath.open('r') as f:
            p = json.load(f)
        self.popsize = p["population"]["size"]
        self.set_base_periods(AverageDurations(**p["durations"]))
        self.set_infection_parameters(
            InfectionParameters(**p["infection"])
        )
        self.set_bgd_mortality(p["bgd_mortality_rate"])
        self.set_groups(p["population"]["groups"])

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
        # FIXME: What if arrays in base_periods are garbage collected?
        self.lib.set_base_periods(*base_periods)
        return None

    def set_infection_parameters(self, p: InfectionParameters):
        self.lib.set_infection_parameters.restype = None
        self.lib.set_infection_parameters.argtypes = [ctypes.c_double] * 4
        self.lib.set_infection_parameters(
            ctypes.c_double(p.v1),
            ctypes.c_double(p.v2),
            ctypes.c_double(p.phi),
            ctypes.c_double(p.epsilon),
    )

    def set_bgd_mortality(self, tau: float):
        self.lib.set_background_mortality.restype = None
        self.lib.set_background_mortality.argtypes = [ctypes.c_double]
        self.lib.set_background_mortality(
            ctypes.c_double(1. - np.exp(- 1. / tau))
        )

    def set_groups(self, groups: list[int]):
        self.lib.set_groups.restype = None
        array_type = ctypes.c_int * len(groups)
        self.lib.set_groups.argtypes = [array_type, ctypes.c_int]
        self.lib.set_groups(array_type(*groups), len(groups))

    #TODO: Add getter for base periods, read arrays from C lib
