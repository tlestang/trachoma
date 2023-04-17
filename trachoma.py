import numpy as np


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


def init_infected(pop_size, fraction, rng):
    newinf = np.zeros(pop_size, dtype=np.bool_)
    ninf = int(fraction * pop_size)
    infected_id = rng.choice(range(0, pop_size), size=ninf, replace=False)
    newinf[infected_id] = True
    return newinf


def getlambdaStep(ages, bact_load):

    y_children = ages < 9 * 52
    o_children = (ages >= 9 * 52) & (ages < 15 * 52)
    adults = ages >= 15 * 52

    totalLoad = np.array([
        np.sum(bact_load[y_children]) / len(y_children),
        np.sum(bact_load[o_children]) / len(o_children),
        np.sum(bact_load[adults]) / len(adults)
    ])
    # [lambda1, lambda2, lambda3]
    prevLambda = bet * (V_1 * totalLoad + V_2 * (totalLoad ** (PHI + 1)))

    a = len(y_children)/POP_SIZE
    b = len(o_children)/POP_SIZE
    c = len(adults)/POP_SIZE
    epsm = 1 - EPSILON
    A = 1 - np.exp([
        - prevLambda[0]*a + prevLambda[1]*epsm*b + prevLambda[2]*epsm*c,
        - prevLambda[0]*a*epsm + prevLambda[1]*b + prevLambda[2]*epsm*c,
        - prevLambda[0]*a*epsm + prevLambda[1]*epsm*b + prevLambda[2]*c,
    ])
    returned = np.ones(params['N'])
    returned[y_children] = A[0]
    returned[o_children] = A[1]
    returned[adults] = A[2]
    return returned


def get_new_infections(susceptibles, ages, bact_load, rng):
    prob = getlambdaStep(ages, bact_load)
    target_size = np.count_nonzero(susceptibles)
    newinf = np.zeros(ages.size, dtype=np.bool_)
    newinf[susceptibles] = rng.uniform(size=nsusceptibles) < prob[susceptibles]
    return newinf


V_1 = 1
V_2 = 2.6
PHI = 1.4
T_LATENT = 2
T_ID = 2
T_D = 2
POP_SIZE = 200
NWEEKS = 52
MAX_AGE = 60 * NWEEKS
MEAN_AGE = 20 * NWEEKS

rng = np.random.default_rng()

infected = np.zeros(POP_SIZE, dtype=np.bool_)
diseased = np.zeros(POP_SIZE, dtype=np.bool_)
clock = np.zeros(POP_SIZE) - 1
infection_counter = np.zeros(POP_SIZE)

newinf = init_infected(POP_SIZE, initial_fraction_infected, rng)
infected = infected | newinf
infection_counter[newinf] = 1

ages = init_ages(POP_SIZE, MAX_AGE, MEAN_AGE, rng)

for t in timesteps:
    transition = np.logical_not(clock.astype(np.bool_))
    new_s = diseased & ~infected & transition
    new_d = infected & diseased & transition
    new_id = infected & ~diseased & transition
    newinf = get_new_infections(~I, ages, bact_load, rng)

    T[new_i] = T_LATENT
    T[new_id] = T_ID
    T[new_d] = T_D

    diseased = diseased & ~new_s | new_d
    infected = infected & ~new_d | new_i

    # housekeeping
    infection_counter[newinf] += 1
    clock -= -1


