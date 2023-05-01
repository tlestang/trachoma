import numpy as np
import init
import parameters as p
from infection import get_new_infections, get_load
import periods


rng = np.random.default_rng()

ages = init.ages(p.POP_SIZE, p.MAX_AGE, p.MEAN_AGE, rng)

# Baseline disease periods
latent_period_base = np.zeros(p.POP_SIZE) + 2
ID_period_base=rng.poisson(lam=p.AV_ID_DURATION, size=p.POP_SIZE)
D_period_base=rng.poisson(lam=p.AV_D_DURATION, size=p.POP_SIZE)

newinf = init.infected(p.POP_SIZE, fraction=0.1, rng=rng)
infected = infected | newinf
infection_counter[newinf] = 1


def compose(f, g):
    return lambda x: f(g(x))

def identity(f)

from functools import reduce

reduce(compose, [f]* n)


for t in range(0,10000):
    transition = np.logical_not(clock.astype(np.bool_))
    new_s = diseased & ~infected & transition
    new_d = infected & diseased & transition
    new_id = infected & ~diseased & transition
    new_i = get_new_infections(~infected, ages, bact_load, rng)

    clock[new_i] = periods.latent_time(latent_period_base[new_i])
    clock[new_id] = periods.id_time(ID_period_base[new_id], infection_counter[new_id])
    clock[new_d] = periods.d_time(D_period_base[new_d], infection_counter[new_d], ages[new_d])

    diseased = diseased & ~new_s | new_d
    infected = infected & ~new_d | new_i

    bact_load[new_d] = 0
    bact_load[new_i] = get_load(infection_counter[new_i])

    # housekeeping
    infection_counter[new_i] += 1
    clock -= -1
