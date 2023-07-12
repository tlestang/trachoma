import pickle
import numpy as np
from trachoma.trachoma_functions import stepF_fixed, Set_inits, Seed_infection
from test_params import params, demog, bet
import init

vals = Set_inits(params=params, demog=demog, sim_params={"N_MDA": 0})

rng = np.random.default_rng()
vals['Age'] = init.ages(params['N'], demog['max_age'], demog['mean_age'], rng)
vals = Seed_infection(params=params, vals=vals)

nt = 100
data = np.zeros((nt, 2))

for i in range(nt):
    vals = stepF_fixed(vals, params, demog, bet)

    children_ages_1_9 = np.logical_and(vals['Age'] < 10 * 52, vals['Age'] >= 52)
    n_children_ages_1_9 = np.count_nonzero(children_ages_1_9)
    n_true_diseased_children_1_9 = np.count_nonzero(vals['IndD'][children_ages_1_9])
    n_true_infected_children_1_9 = np.count_nonzero(vals['IndI'][children_ages_1_9])
    data[i, 0] = n_true_diseased_children_1_9 / n_children_ages_1_9
    data[i, 1] = n_true_infected_children_1_9 / n_children_ages_1_9 # 

with open("reference.csv", "w") as f:
    np.savetxt(f, data, delimiter=",", fmt="%.3f")
