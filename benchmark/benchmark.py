from step import stepF_fixed
from load_parameters import params, demog
from init import Set_inits, Seed_infection

BET = 0.21

vals = Set_inits(params=params, demog=demog)  # set initial conditions
vals = Seed_infection(params=params, vals=vals)  # seed infection

for t in range(0, 10000):
    vals = stepF_fixed(vals, params, demog, BET)
