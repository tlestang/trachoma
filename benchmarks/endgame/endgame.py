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


def main():
    events = process_scenario_definition(SCENARIO_PATH)

    p = get_params(PARAMETERS_PATH)

    rng = np.random.default_rng()

    ages = init.ages(
        p.population_size, p.max_age, p.mean_age, rng
    )
    base_periods = simulate.set_base_periods(p, rng)
    simulate.set_infection_parameters(p)
    simulate.set_background_mortality(p)
    simulate.set_groups(p)
    p = Population(ages, latent_base=base_periods[0], rng=rng)

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
    threads = [
        Thread(target=simulate.simulate, args=(task, events))
        for task in tasks
    ]
    for t in threads:
        t.start()


if __name__ == "__main__":
    main()
