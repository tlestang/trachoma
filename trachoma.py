import numpy as np

def init_ages(params, max_age, mean_age):
    ages = np.arange(max_age) + 1
    age_prob =  (
        (1. - np.exp(-1. / mean_age)) *
        np.exp( - ages / mean_age)
    )
    age_prob[-1] = 1 - age_prob[:-1].sum()

    return rng.choice(
        ages, p=age_prob,
        size=pop_size,
        replace=True
    )

def init_infected(pop_size, fraction, rng):
    I = np.zeros(pop_size, dtype=np.bool_)
    T_latent = np.zeros(pop_size)
    ninf = int(fraction * pop_size)
    infected_id = rng.integers(low=0, high=pop_size, size=ninf)
    I[infected_id] = True
    T_latent[infected_id] = LATENT_PERIOD
    return I, T_latent


LATENT_PERIOD = 2
POP_SIZE = 200
NWEEKS = 52
MAX_AGE = 60 * NWEEKS
MEAN_AGE = 20 * NWEEKS

rng = np.random.default_rng()

ages = init_ages(pop_size, MAX_AGE, MEAN_AGE, rng)
I, T_latent = init_infected(pop_size, fraction, rng)
