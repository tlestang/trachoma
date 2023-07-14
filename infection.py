import numpy as np
from parameters import PHI, BET, V_1, V_2, POP_SIZE, EPSILON


def getlambdaStep(ages, bact_load):

    y_children = ages < 9 * 52
    o_children = (ages >= 9 * 52) & (ages < 15 * 52)
    adults = ages >= 15 * 52

    totalLoad = np.array([
        np.sum(bact_load[y_children]) / np.count_nonzero(y_children),
        np.sum(bact_load[o_children]) / np.count_nonzero(o_children),
        np.sum(bact_load[adults]) / np.count_nonzero(adults)
    ])
    # [lambda1, lambda2, lambda3]
    prevLambda = BET * (V_1 * totalLoad + V_2 * (totalLoad ** (PHI + 1)))

    a = np.count_nonzero(y_children) / POP_SIZE
    b = np.count_nonzero(o_children) / POP_SIZE
    c = np.count_nonzero(adults) / POP_SIZE
    epsm = 1 - EPSILON
    A = 1 - np.exp([
        - prevLambda[0]*a - prevLambda[1]*epsm*b - prevLambda[2]*epsm*c,
        - prevLambda[0]*a*epsm - prevLambda[1]*b - prevLambda[2]*epsm*c,
        - prevLambda[0]*a*epsm - prevLambda[1]*epsm*b - prevLambda[2]*c,
    ])
    returned = np.ones(POP_SIZE)
    returned[y_children] = A[0]
    returned[o_children] = A[1]
    returned[adults] = A[2]
    return returned


def get_new_infections(susceptibles, ages, bact_load, rng):
    prob = getlambdaStep(ages, bact_load)
    target_size = np.count_nonzero(susceptibles)
    newinf = np.zeros(ages.size, dtype=np.bool_)
    newinf[susceptibles] = np.random.uniform(size=target_size) < prob[susceptibles]
    return newinf


def get_load(No_Inf):

    '''
    Function to scale bacterial load according to infection history
    '''

    b1 = 1
    ep2 = 0.114

    return b1 * np.exp((No_Inf - 1) * - ep2)

