import numpy as np
from ntdmc_trachoma.state import Population


def init_ages(pop_size, max_age, mean_age, rng):
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


AGES = np.array(
    [13, 19, 11, 23,  1, 59, 27, 60,  7, 11, 13,  3, 10,  2, 21, 34]
)
POPSIZE = len(AGES)


def test_seed_infection():
    pop = Population(AGES)

    latent_period_func = lambda count: (count + 2)
    bact_load_func = lambda count: (count + 0.7)

    pop.seed_infection(
        k=4,
        latent_period_func=latent_period_func,
        bact_load_func=bact_load_func,
    )

    # Here we assume the correct fraction of individuals have been
    # selected for infection.

    expected_clock = np.zeros(POPSIZE, dtype=np.int32) - 1
    expected_clock[pop.lat] = 2
    np.testing.assert_equal(pop.clock, expected_clock)

    expected_bact_load = np.zeros(POPSIZE, dtype=np.float64)
    expected_bact_load[pop.lat] = 0.7
    np.testing.assert_equal(pop.bact_load, expected_bact_load)
    
    



