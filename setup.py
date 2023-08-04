import init
import periods
from infection import get_load
import numpy as np
from parameters import AV_D_DURATION, AV_ID_DURATION

n = 4096
rng = np.random.default_rng()
ages = np.sort(init.ages(n, 60 * 52, 20 * 52, rng)).astype(np.int32)

clock = np.zeros(n, dtype=np.int32) - 1
bact_load = np.zeros(n, dtype=np.float64)
prob = np.zeros(n, dtype=np.float64)
infection_counter = np.zeros(n, dtype=np.int32)
latent = init.infected(n, 0.01, rng)
infected = latent.copy()
diseased = np.zeros(n, dtype=np.bool_)

latent_base = np.zeros(n, dtype=np.int32) + 1
ID_base = rng.poisson(lam=AV_ID_DURATION, size=n).astype(np.int32)
D_base = rng.poisson(lam=AV_D_DURATION, size=n).astype(np.int32)

clock[latent] = periods.latent_time(latent_base[latent])
bact_load[latent] = get_load(infection_counter[latent])
infection_counter[latent] = 1

vals = {}
vals['Age'] = ages
vals['IndI'] = infected.astype(int)
vals['IndD'] = diseased.astype(int)
vals['T_latent'] = np.zeros(n)
vals['T_D'] = np.zeros(n)
vals['T_ID'] = np.zeros(n)

vals['T_latent'][latent] = clock[latent].copy()
vals['No_Inf'] = infection_counter.copy()
vals['bact_load'] = bact_load.copy()
vals['Ind_ID_period_base'] = ID_base.copy()
vals['Ind_D_period_base'] = D_base.copy()
vals['Ind_latent'] = latent_base.copy()
