import numpy as np


def init_ages(params, max_age, mean_age):
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
    I = np.zeros(pop_size, dtype=np.bool_)
    T_latent = np.zeros(pop_size)
    ninf = int(fraction * pop_size)
    infected_id = rng.integers(low=0, high=pop_size, size=ninf)
    I[infected_id] = True
    T_latent[infected_id] = LATENT_PERIOD
    return I, T_latent


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


def get_new_infections(I, ages, bact_load, rng):
    prob = getlambdaStep(ages, bact_load)
    target = ~I
    target_size = np.count_nonzero(target)
    newinf = np.zeros(POP_SIZE, dtype=np.bool_)
    newinf[target] = rng.uniform(size=target_size) < prob[target]
    return newinf



V_1 = 1
V_2 = 2.6
PHI = 1.4
LATENT_PERIOD = 2
POP_SIZE = 200
NWEEKS = 52
MAX_AGE = 60 * NWEEKS
MEAN_AGE = 20 * NWEEKS

rng = np.random.default_rng()
infection_counter = np.zeros(POP_SIZE)

ages = init_ages(POP_SIZE, MAX_AGE, MEAN_AGE, rng)

I = np.zeros(POP_SIZE, dtype=np.bool_)
D = np.zeros(POP_SIZE, dtype=np.bool_)
T = np.zeros(POP_SIZE) - 1

init_infected(I, T)
infection_counter[I] = 1

for t in timesteps:

    new_s = D & ~I & ~T
    new_d = I & D & ~T
    new_id = I & ~D & ~T
    new_i = get_new_infections(I, ages, bact_load, rng)

    T[new_i] = T_latent
    T[new_id] = T_ID
    T[new_d] = T_D

    D = D & ~new_s | new_d
    I = I & ~new_d | new_i

    # housekeeping
    infection_counter[new_i] += 1
    T -= -1


