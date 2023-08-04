import setup

from trachoma.trachoma_functions import stepF_fixed

from test_params import params, demog, bet

params['N'] = 4096
vals = setup.vals.copy()
for i in range(1000):
    vals = stepF_fixed(vals=vals, params=params, demog=demog, bet=bet)
