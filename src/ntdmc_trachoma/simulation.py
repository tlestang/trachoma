from collections import namedtuple
from copy import deepcopy
import ctypes
from importlib import util
import json
from pathlib import Path

import numpy as np

from .parameters import (
    AverageDurations, InfectionParameters,
    PopulationParameters, BasePeriods
)
import ntdmc_trachoma.setup_core as setup_core
import ntdmc_trachoma.process_scenario_definition as scenario
import ntdmc_trachoma.init as init
from .state import Population


LIBTRACHO_PATH = util.find_spec(
    "ntdmc_trachoma.libtrachoma"
).origin


class Simulation:
    def __init__(self, parameters_filepath: Path):
        self.rng = np.random.default_rng()
        self.lib = ctypes.CDLL(LIBTRACHO_PATH)
        params = self.load_parameters(parameters_filepath)
        self.base_periods = BasePeriods(
            latent=np.array(
                [params["durations"].I] * params["pop"].size, dtype=np.int32
            ),
            ID=self.rng.poisson(
                lam=params["durations"].ID, size=params["pop"].size,
            ).astype(np.int32),
            D=self.rng.poisson(
                lam=params["durations"].D, size=params["pop"].size
            ).astype(np.int32),
        )
        ages = init.ages(
            params["pop"].size,
            params["pop"].max_age,
            params["pop"].average_age,
            self.rng,
        )
        self.pop = Population(
            ages=ages,
            latent_base=self.base_periods.latent,
            rng=self.rng,
            lib=self.lib,
        )

        setup_core.set_base_periods(self.lib, self.base_periods)
        setup_core.set_infection_parameters(self.lib, params["inf"])
        setup_core.set_bgd_mortality(
            self.lib, params["pop"].bgd_mortality_rate
        )
        setup_core.set_groups(self.lib, params["pop"].groups)

    def load_parameters(self, parameters_filepath: Path):
        with parameters_filepath.open('r') as f:
            p = json.load(f)
        return {
            "pop": PopulationParameters(**p["population"]),
            "inf": InfectionParameters(**p["infection"]),
            "durations": AverageDurations(**p["durations"]),
        }

    #TODO: Add getter for base periods, read arrays from C lib

    def simulate(self, scenario_filepath: Path, betavals: list[float]):
        events = scenario.process(scenario_filepath)
        for betaval in betavals:
            pop = deepcopy(self.pop)
            for e, e_next in zip(events[:-1], events[1:]):
                self.lib.advance(
                    pop,
                    e_next[0] - e[0],
                    ctypes.c_double(betaval)
                )
                e_next[1](pop)
