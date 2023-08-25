import ctypes
from copy import deepcopy
import json
from pathlib import Path
from importlib import util

import numpy as np

import ntdmc_trachoma.init as init
from ntdmc_trachoma.state import Population
import ntdmc_trachoma.simulate as simulate
from ntdmc_trachoma.parameters import AverageDurations, InfectionParameters
from ntdmc_trachoma.process_scenario_definition import process_scenario_definition
import ntdmc_trachoma.setup_core as setup_core

SCENARIO_PATH = Path("scenario.json")
BETA_PATH = Path("k_values.txt")
PARAMETERS_PATH = Path("parameters.json")


def main():
    events = process_scenario_definition(SCENARIO_PATH)
    with PARAMETERS_PATH.open('r') as f:
        p = json.load(f)

    rng = np.random.default_rng()

    pop_params = p["population"]
    ages = init.ages(
        pop_params["size"], pop_params["max_age"],
        pop_params["mean_age"], rng
    )

    LIBTRACHO_PATH = util.find_spec(
        "ntdmc_trachoma.libtrachoma"
    ).origin
    lib = ctypes.CDLL(LIBTRACHO_PATH)
    base_periods = setup_core.set_base_periods(
        lib,
        AverageDurations(**p["durations"]),
        p["population"]["size"], rng
    )
    setup_core.set_infection_parameters(
            lib,
            InfectionParameters(**p["infection"])
        )
    setup_core.set_bgd_mortality(lib, p["bgd_mortality_rate"])
    setup_core.set_groups(lib, p["population"]["groups"])

    p = Population(ages, latent_base=base_periods[0], rng=rng, lib=lib)

    with BETA_PATH.open() as f:
        kvalues = [
            float(line) for line in f
        ]
    NTHREADS = 1
    task_size = len(kvalues) // NTHREADS
    tasks = [
        [(deepcopy(p), k) for k in kvalues[i * task_size:(i + 1) * task_size]]
        for i in range(0, NTHREADS)
    ]
    if (len(kvalues) % NTHREADS):
        tasks[-1].extend(
            [(deepcopy(p), k) for k in kvalues[NTHREADS * task_size:]]
        )
    simulate.simulate(tasks[0], events)

if __name__ == "__main__":
    main()
