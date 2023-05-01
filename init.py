import numpy as np


def ages(pop_size, max_age, mean_age, rng):
    """Returns age array"""
    ages = np.arange(max_age) + 1
    age_prob = (
        (1. - np.exp(-1. / mean_age)) *
        np.exp(- ages / mean_age)
    )
    age_prob[-1] = 1 - age_prob[:-1].sum()

    return rng.choice(
        ages, p=age_prob,
        size=pop_size,
        replace=True
    )


def infected(pop_size, fraction, rng):
    newinf = np.zeros(pop_size, dtype=np.bool_)
    ninf = int(fraction * pop_size)
    infected_id = rng.choice(range(0, pop_size), size=ninf, replace=False)
    newinf[infected_id] = True
    return newinf


