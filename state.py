import numpy as np

import periods
from infection import get_new_infections, get_load
import parameters as p
from init import infected

class Population:

    def __init__(self, ages):

        self.rng = np.random.default_rng()

        self.size = len(ages)
        self.ages = ages

        self.diseased = np.zeros(self.size, dtype=np.bool_)
        self.clock = np.zeros(self.size) - 1
        self.infection_counter = np.zeros(self.size)
        self.bact_load = np.zeros(self.size)


        # Baseline disease periods
        self.latent_period_base = np.zeros(p.POP_SIZE) + 1
        self.ID_period_base=self.rng.poisson(lam=p.AV_ID_DURATION, size=p.POP_SIZE)
        self.D_period_base=self.rng.poisson(lam=p.AV_D_DURATION, size=p.POP_SIZE)

        self.latent = infected(self.size, 0.01, self.rng)
        self.infected = self.latent.copy()
        self.clock[self.latent] = periods.latent_time(self.latent_period_base[self.latent])
        self.bact_load[self.latent] = get_load(self.infection_counter[self.latent])
        self.infection_counter[self.latent] = 1


    def tick(self):
        self.clock += -1
        transition = np.logical_not(self.clock.astype(np.bool_))
        new_s = self.diseased & ~self.infected & transition
        new_d = self.infected & ~self.latent & transition
        clearinf = self.latent & transition
        new_i = get_new_infections(~self.infected, self.ages, self.bact_load, self.rng)
        idx = np.array(range(0, 1000))
        import pdb; pdb.set_trace()

        self.clock[new_d] = periods.d_time(
            self.D_period_base[new_d], self.infection_counter[new_d], self.ages[new_d]
        )
        self.clock[clearinf] = periods.id_time(
            self.ID_period_base[clearinf], self.infection_counter[clearinf]
        )
        self.clock[new_i] = periods.latent_time(self.latent_period_base[new_i])

        self.diseased = self.diseased & ~new_s | clearinf
        self.infected = self.infected & ~new_d | new_i
        self.latent = self.latent & ~clearinf | new_i

        # bacterial load
        self.bact_load[new_d] = 0
        self.bact_load[clearinf] = get_load(self.infection_counter[clearinf])

        # housekeeping
        self.infection_counter[new_i] += 1
        self.ages += 1

        # Death
        dead = (np.random.uniform(0, 1, self.size) < p.bgrd_death_rate) | (self.ages > p.MAX_AGE)
        self.clock[dead] = -1
        self.diseased[dead] = False
        self.infected[dead] = False
        self.latent[dead] = False
        self.infection_counter[dead] = 0
        self.ages[dead] = 0

    def prevalence(self, age_bounds):
        i = (self.ages >= age_bounds[0]) | (self.ages < age_bounds[1])
        return np.count_nonzero(self.diseased[i]) / np.count_nonzero(i)

    def infection(self, age_bounds):
        i = (self.ages >= age_bounds[0]) | (self.ages < age_bounds[1])
        return np.count_nonzero(self.infected[i]) / np.count_nonzero(i)

    def group_size(self, age_bounds):
        return np.count_nonzero(
            (self.ages >= age_bounds[0]) | (self.ages < age_bounds[1])
        )

    def yearly_high_infection_count(threshold):
        """Return propotion of individuals with no of infection
        greater than threshold, per year group
        """
        i = self.infection_counter > threshold

        # Cast weights to integer to be able to count
        a, _ = np.histogram(self.age, bins=p.MAX_AGE, weights=i.astype(int))
        return a / self.size
