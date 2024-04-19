import numpy as np
from ntdmc_trachoma.state import Population

AGES = np.array(
    [13, 19, 11, 23, 1, 59, 27, 60, 7, 11, 13, 3, 10, 2, 21, 34]
)
POPSIZE = len(AGES)

rng = np.random.default_rng()
latent_base_periods = rng.poisson(lam=2, size=POPSIZE).astype(np.int32)

def test_seed_infection():
    pop = Population(AGES)

    latent_period_func = lambda count, base: (count + 2) * base
    bact_load_func = lambda count: (count + 0.7)

    pop.seed_infection(
        k=4,
        latent_periods=latent_base_periods,
        latent_period_func=latent_period_func,
        bact_load_func=bact_load_func,
    )

    # Here we assume the correct fraction of individuals have been
    # selected for infection.

    expected_clock = np.zeros(POPSIZE, dtype=np.int32) - 1
    expected_clock[pop.lat] = 2 * latent_base_periods[pop.lat]
    np.testing.assert_equal(pop.clock, expected_clock)

    expected_bact_load = np.zeros(POPSIZE, dtype=np.float64)
    expected_bact_load[pop.lat] = 0.7
    np.testing.assert_equal(pop.bact_load, expected_bact_load)
    
    



