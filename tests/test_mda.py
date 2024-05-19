from unittest.mock import patch
import numpy as np
from ntdmc_trachoma.state import Population
from ntdmc_trachoma.mda import MDA, MDA_data

rng = np.random.default_rng()
ages = rng.integers(low=1, high=60, size=16)
a, b = (0.2, 0.3)


def test_draw_tment_prob():
    pop = Population(ages)
    current_ages = pop.ages[np.argsort(pop.indexes)]
    data = MDA_data(
        ages=current_ages,  # Sorted according to indexes
        treatment_probability=rng.beta(a, b, size=16),
        beta_dist_params=(a, b),
    )
    # Pretend that individuals 5, 9, 13, 14 and 15 were reset since
    # the last MDA event.
    pop.indexes = np.array(
        [13, 14, 15, 5, 9, 0, 1, 2, 3, 4, 6, 7, 8, 10, 11, 12],
        dtype=np.int32,
    )
    # Reorder array according to resets and age all individuals by
    # arbitrary duration (e.g. 1 week)
    pop.ages = pop.ages[pop.indexes] + 1
    # Set age of reset individuals to something smaller than their
    # previous ages before being reset, e.g. 0
    pop.ages[0:5] = 0
    with patch("ntdmc_trachoma.mda.np.random.default_rng") as mock:
        expected_tment_probs = rng.beta(a, b, size=5)
        mock.return_value.beta.return_value = expected_tment_probs
        mda = MDA(1., 1., 1.)
    mda.a, mda.b = (a, b)
    # Remember that tment_prob is sorted according to increasing
    # individual indexes and not according to ages.
    tment_prob = mda.draw_tment_prob(
        pop.ages[np.argsort(pop.indexes)], data
    )
    np.testing.assert_allclose(
        tment_prob[[5, 9, 13, 14, 15]], expected_tment_probs,
    )


def test_draw_tment_prob_parameter_change():
    #  In this test we don't reset any individuals since
    #  treatment probabilties are entirely redrawn.
    pop = Population(ages)
    current_ages = pop.ages[np.argsort(pop.indexes)]
    data = MDA_data(
        ages=current_ages,  # Sorted according to indexes
        #  Set tment prob to arbitrary array, it only matters that we can
        #  check that individuals are sorted in the same order according to
        # to redrawn probabilities.
        treatment_probability=rng.permutation(pop.indexes).astype(np.float64),
        beta_dist_params=(a, b),
    )
    mda = MDA(0.7, 0.9, 1.)  # Beta dist params different from a, b above
    tment_prob = mda.draw_tment_prob(current_ages, data)
    #  New values for treatment probabilties should be different,
    #  because redrawn from new values of a and b.
    with np.testing.assert_raises(AssertionError):
        np.testing.assert_allclose(tment_prob, data.treatment_probability)
    #  Order according to treatment probability should be kept the same.
    np.testing.assert_array_equal(
        np.argsort(tment_prob),
        np.argsort(data.treatment_probability),
    )
