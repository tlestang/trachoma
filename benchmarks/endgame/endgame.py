import sys
sys.path.append("../../")
from copy import deepcopy
from pathlib import Path
from threading import Thread

import numpy as np

import init
from state import Population
import simulate
from process_scenario_definition import process_scenario_definition
from parameters import get_params

SCENARIO_PATH = Path("scenario.json")
BETA_PATH = Path("k_values.txt")
PARAMETERS_PATH = Path("parameters.json")

rng = np.random.default_rng()


def main():
    # We start by reading in the parameters and setting some global
    # pointers on the C side.
    p = get_params(PARAMETERS_PATH)
    base_periods = simulate.set_base_periods(p, rng)
    simulate.set_infection_parameters(p)
    simulate.set_background_mortality(p)
    simulate.set_groups(p)

    # The next step is to create a Population object
    ages = init.ages(
        p.population_size, p.max_age, p.mean_age, rng
    )
    pop = Population(ages, latent_base=base_periods[0], rng=rng)

    # In this script we use a fixed number of threads NTHREADS and
    # distribute the work (beta values) across threads
    with BETA_PATH.open() as f:
        kvalues = [
            float(line) for line in f
        ]
    NTHREADS = 1
    task_size = len(kvalues) // NTHREADS  # Nb of beta values per thread

    # tasks is a list of NTHREADS lists that contain the population
    # object and values of beta to work with.  Note that the
    # population is not reset between beta values, which is
    # incorrect at the moment.
    tasks = [
        [
            (deepcopy(pop), k)
            for k in kvalues[i * task_size:(i + 1) * task_size]
        ]
        for i in range(0, NTHREADS)
    ]
    if (len(kvalues) % NTHREADS):
        # If the total number of beta values is not a multiple of
        # NTHREADS, we just add them to the last thread's work.
        tasks[-1].extend(
            [(deepcopy(p), k) for k in kvalues[NTHREADS * task_size:]]
        )

    # The last remaining bit of pre-processing is reading the
    # scenarion file and prepare the list of events (MDA, vaccine...).
    events = process_scenario_definition(SCENARIO_PATH)

    # Finally we're ready to launch the threads.
    threads = [
        Thread(target=simulate.simulate, args=(task, events))
        for task in tasks
    ]
    for t in threads:
        t.start()


if __name__ == "__main__":
    main()
