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

def getlambdaStep(ages, bact_load):

    y_children = ages < 9 * 52
    o_children = (ages >= 9 * 52) & (ages < 15 * 52)
    adults = ages >= 15 * 52

    totalLoad = np.array([np.sum(bact_load[y_children]) / len(y_children),
    np.sum(bact_load[o_children]) / len(o_children), np.sum(bact_load[adults]) / len(adults)])
    prevLambda = bet * (V_1 * totalLoad + V_2 * (totalLoad ** (PHI + 1)))

    a = len(y_children)/POP_SIZE
    b = len(o_children)/POP_SIZE
    c = len(adults)/POP_SIZE
    epsm = 1 - EPSILON
    A = [
        prevLambda[0]*a + prevLambda[1]*epsm*b + prevLambda[2]*epsm*c,
        prevLambda[0]*a*epsm + prevLambda[1]*b + prevLambda[2]*epsm*c,
        prevLambda[0]*a*epsm + prevLambda[1]*epsm*b + prevLambda[2]*c,
    ]
    returned = np.ones(params['N'])
    returned[y_children] = A[0]
    returned[o_children] = A[1]
    returned[adults] = A[2]
    return returned


def infected_mask(IndI, ages, bact_load, rng):
    lam = getlambdaStep(ages, bact_load)
    infection_prob = 1 - np.exp(-lam)    
    return (
        (rng.uniform(size=POP_SIZE) < infection_prob) &
        ~IndI
    )

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

I = [[] for i in range(latent_period)]
ID = [[] for i in range(ID_period)]
D = [[] for i in range(D_period)]

I[-1] = list(rng.integers(low=0, high=POP_SIZE, size=ninf))
infection_counter[I[-1]] += 1

for t in timesteps:

    iD = t % D_period
    iID = t % ID_period
    iI = t % latent_period
    
    D[iD] = ID[iID]; IndI[D[iD]] = False # ID -> D
    bact_load[D[iD]] = get_bact_load()
    ID[iID] = I[iI] # I -> ID

    infected = infected_mask(IndI, ages, bact_load)
    I[iI] = ids[infected]
    IndI[infected] = True
    infection_counter[infected] += 1




