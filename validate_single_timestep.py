import numpy as np
import pickle
import init
import parameters as p
from state import Population

from trachoma.trachoma_functions import stepF_fixed

from test_params import params, demog, bet

with open("out.p", "rb") as f:
    vals = pickle.load(f)

pop = Population(vals['Age'].copy())
# pop.diseased = (vals['T_D'] > 0).astype(bool)
# pop.ID = (vals['T_ID'] > 0).astype(bool)
# pop.latent = (vals['T_latent'] > 0).astype(bool)

pop.diseased = vals['IndD'].astype(bool)
pop.infected = vals['IndI'].astype(bool)
pop.latent = (vals['T_latent'] > 0).astype(bool)

# dtv = (~vals['IndI'].astype(bool) & vals['IndD'].astype(bool))
# pop.clock[dtv] = vals['T_D'][dtv]
# itv = (~vals['IndD'].astype(bool) & vals['IndI'].astype(bool))
# pop.clock[itv] = vals['T_latent'][itv]
# idtv = (vals['IndD'].astype(bool) & vals['IndI'].astype(bool))
# pop.clock[idtv] = vals['T_ID'][idtv]

pop.clock = vals['T_latent'] + vals['T_ID'] + vals['T_D']

pop.infection_counter = vals['No_Inf'].copy()
pop.bact_load = vals['bact_load'].copy()
pop.ID_period_base = vals['Ind_ID_period_base'].copy()
pop.D_period_base = vals['Ind_D_period_base'].copy()
pop.latent_period_base = vals['Ind_latent'].copy()

np.random.seed(0)
pop.tick()

np.random.seed(0)
vals = stepF_fixed(vals=vals, params=params, demog=demog, bet=bet)

print(f"Children 1 - 9 infection: {pop.infection((52, 10 * 52))}")
print(f"Children 1 - 9 prevalence: {pop.prevalence((52, 10 * 52))}")
