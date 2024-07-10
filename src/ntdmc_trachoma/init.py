import numpy as np


# TODO: Document initial age distribution
def ages(pop_size, max_age, mean_age, rng):
    """Returns age array"""
    allowed_age_values = 1 + np.arange(max_age)

    ages_normed = allowed_age_values / mean_age
    weights = np.empty(len(ages_normed))
    weights[:-1] = np.exp(-ages_normed[:-1]) - np.exp(-ages_normed[1:])
    weights[-1] = 1. - np.sum(weights[:-1])

    return rng.choice(
        allowed_age_values, p=weights,
        size=pop_size,
        replace=True
    )


# TODO: Document initial infection
def infected(pop_size, fraction, rng):
    newinf = np.zeros(pop_size, dtype=np.bool_)
    ninf = int(fraction * pop_size)
    infected_id = rng.choice(range(0, pop_size), size=ninf, replace=False)
    newinf[infected_id] = True
    return newinf
