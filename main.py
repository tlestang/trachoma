import argparse
from copy import deepcopy
from pathlib import Path
from threading import Thread

import numpy as np

import init
from state import Population
import simulate
from process_scenario_definition import process_scenario_definition
from parameters import get_params


def parse_args():
    parser = argparse.ArgumentParser(
        prog="Trachoma simulation",
        description="Simulation of the dynamics of trachoma disease"
    )
    parser.add_argument(
        "-s", "--scenario", type=Path, required=True, dest="scenario_path"
    )
    parser.add_argument(
        "-p", "--parameters", type=Path, required=True, dest="parameters_path"
    )
    parser.add_argument(
        "-k", type=Path, required=True, dest="k_path"
    )
    return parser.parse_args()


def main(args):
    events = process_scenario_definition(args.scenario_path)

    p = get_params(args.parameters_path)

    rng = np.random.default_rng()

    ages = init.ages(
        p.population_size, p.max_age, p.mean_age, rng
    )
    base_periods = simulate.set_base_periods(p, rng)
    simulate.set_infection_parameters(p)
    simulate.set_background_mortality(p)
    simulate.set_groups(p)
    p = Population(ages, latent_base=base_periods[0], rng=rng)

    with args.k_path.open() as f:
        kvalues = [
            float(line) for line in f
        ]
    pop_objs = [deepcopy(p) for k in kvalues]
    threads = [
        Thread(target=simulate.simulate, args=(p, events, k))
        for k, p in zip(kvalues, pop_objs)
    ]
    for t in threads:
        t.start()


if __name__ == "__main__":
    args = parse_args()
    main(args)
