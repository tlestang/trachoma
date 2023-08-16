from pathlib import Path
from datetime import datetime, timedelta
from numpy.random import default_rng
from mda import start_sim, end_sim
from parameters import get_params
import init
import simulate
from state import Population

p = get_params(Path("./benchmark_parameters"))

rng = default_rng()

ages = init.ages(
    p.population_size, p.max_age, p.mean_age, rng
)
base_periods = simulate.set_base_periods(p, rng)
simulate.set_infection_parameters(p)
simulate.set_background_mortality(p)
simulate.set_groups(p)
pop = Population(ages, latent_base=base_periods[0], rng=rng)

events = [
    (start_sim, datetime(2008, 1, 1)),
    (end_sim, datetime(2008, 1, 1) + timedelta(weeks=520)),
]


simulate(pop, events, 0.21)
