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
from ntdmc_trachoma.state import Population
from ntdmc_trachoma.output import Output


LIBTRACHO_PATH = util.find_spec(
    "ntdmc_trachoma.libtrachoma"
).origin


class Simulation:
    """Encapsulate top-level facilities for managing a trachoma simulation.

    The ``Simulation`` class provides is the primary entry point to
    handling a trachoma simulation over a range of beta parameters
    values.  It is responsible for reading parameters, loading and
    initialising the underlying C library, as well as providing an
    interface for writing simulation results on disk.

    :param parameters_filepath: The location of a valid parameter
        file. See :doc:`/parameters`.
    :type parameters_filepath: pathlib.Path

    Example
    ~~~~~~~

    >>> from pathlib import Path
    >>> from ntdmc_trachoma import Simulation
    >>> sim = Simulation(Path("parameters.json"))
    >>> betavals = [0.2, 0.3, 0.4]
    >>> sim.simulate(Path("scenario.json"), betavals, record=False)

    """
    def __init__(self, parameters_filepath: Path):
        self.rng = np.random.default_rng()
        self.lib = ctypes.CDLL(LIBTRACHO_PATH)
        params = self.load_parameters(parameters_filepath)

        # We set base periods as a instance attribute to make sure a
        # reference the arrays is kept, because pointers in the C
        # library will be set to point to them and we don't want the
        # arrays to be garbage collected.
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
        # TODO: Register groups as instance attribute otherwise it's
        # gcollected.  Should group array not allocated at the C level
        # instead?
        self.groups = setup_core.set_groups(self.lib, params["pop"].groups)

    def load_parameters(self, parameters_filepath: Path):
        with parameters_filepath.open('r') as f:
            p = json.load(f)
        return {
            "pop": PopulationParameters(**p["population"]),
            "inf": InfectionParameters(**p["infection"]),
            "durations": AverageDurations(**p["durations"]),
        }


    #TODO: Add getter for base periods, read arrays from C lib
    def simulate(self, scenario_filepath: Path, betavals: list[float], record=False):
        """Simulate model over a given scenario, for a range of
        spread parameter values.

        This method does not modify the Simulation's instance ``pop``
        attribute.  Instead, a copy is made before passing it down to
        the function stepping the model, which does alter this copy.
        On the contrary, the `output` attribute is modified, adding
        new values generated at each iteration.  In other words, the
        ``simulate`` method has no side effects other than recording
        simulation ouput in the `Simulation.output` attribute.

        :param scenario_filepath: location of a valid scenario
            file. See :doc:`/scenarios`.
        :type scenario_filepath: pathlib.Path
        :param betavals: The list of beta parameter values to run the
            simulation for.
        :type betavals: list[float]
        """
        # TODO: Even with record=False we need to record the final
        # state of the population
        events = scenario.process_scenario_definition(scenario_filepath)
        nsteps = sum(
            [
                e_next[0] - e[0]
                for e, e_next in zip(events[:-1], events[1:])
            ]
        )
        self.output = (
            Output(
                popsize=len(self.pop.ages),
                max_records=nsteps * len(betavals),
            )
        ) if record else None
        for betaval in betavals:
            pop = deepcopy(self.pop)
            for e, e_next in zip(events[:-1], events[1:]):
                self.lib.apply_rules(
                    pop,
                    self.output,
                    e_next[0] - e[0],
                    ctypes.c_double(betaval)
                )
                e_next[1](pop)
