import ctypes
from importlib import util
import json
from pathlib import Path

import numpy as np

from .parameters import AverageDurations, InfectionParameters
import ntdmc_trachoma.setup_core as setup_core


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

        setup_core.set_base_periods(
            self.lib,
            AverageDurations(**p["durations"]),
            self.popsize, self.rng
        )
        setup_core.set_infection_parameters(
            self.lib,
            InfectionParameters(**p["infection"])
        )
        setup_core.set_bgd_mortality(self.lib, p["bgd_mortality_rate"])
        setup_core.set_groups(self.lib, p["population"]["groups"])

    #TODO: Add getter for base periods, read arrays from C lib
