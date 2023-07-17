import numpy as np
import pickle
import init
import parameters as p
from state import Population
from infection import get_load

from trachoma.trachoma_functions import stepF_fixed

from test_params import params, demog, bet

rng = np.random.default_rng()
ages = init.ages(10, p.MAX_AGE, p.MEAN_AGE, rng)
pop = Population(ages)
pop.infected = np.array(
    [True, False, True, False, False, True, True, False, False, True]
)
pop.diseased = np.array(
    [True, False, True, True,  True,  False,False,False, False, False]
)
pop.clock = np.array(
    [1,    1,     2,    3,     1,     1,    2,    1,     1,     1]
)
pop.infection_counter = (pop.infected | pop.diseased).astype(int)
pop.bact_load[pop.infected] = get_load(pop.infection_counter[pop.infected])

vals = {}
vals['IndD'] = pop.diseased.copy()
vals['IndI'] = pop.infected.copy()
vals['Age'] = ages.copy()
vals['NoInf'] = pop.infection_counter.copy()
vals['bact_load'] = pop.bact_load.copy()
vals['T_latent'] = np.array([0, 0, 0, 0, 0, 1, 2, 0, 0, 1])
vals['T_ID'] = np.array([1, 0, 2, 0, 0, 0, 0, 0, 0, 0])
vals['T_D'] = np.array([0, 0, 0, 3, 1, 0, 0, 0, 0, 0])

# dtv = (~vals['IndI'].astype(bool) & vals['IndD'].astype(bool))
# pop.clock[dtv] = vals['T_D'][dtv]
# itv = (~vals['IndD'].astype(bool) & vals['IndI'].astype(bool))
# pop.clock[itv] = vals['T_latent'][itv]
# idtv = (vals['IndD'].astype(bool) & vals['IndI'].astype(bool))
# pop.clock[idtv] = vals['T_ID'][idtv]

# pop.infection_counter = vals['No_Inf'].copy()
# pop.bact_load = vals['bact_load'].copy()
# pop.ID_period_base = vals['Ind_ID_period_base'].copy()
# pop.D_period_base = vals['Ind_D_period_base'].copy()
# pop.latent_period_base = vals['Ind_latent'].copy()

np.random.seed(0)
pop.tick()

np.random.seed(0)
vals = stepF_fixed(vals=vals, params=params, demog=demog, bet=bet)

print(f"Children 1 - 9 infection: {pop.infection((52, 10 * 52))}")
print(f"Children 1 - 9 prevalence: {pop.prevalence((52, 10 * 52))}")
