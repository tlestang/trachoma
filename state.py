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
        self.latent_period_base = np.zeros(p.POP_SIZE) + 2
        self.ID_period_base=self.rng.poisson(lam=p.AV_ID_DURATION, size=p.POP_SIZE)
        self.D_period_base=self.rng.poisson(lam=p.AV_D_DURATION, size=p.POP_SIZE)

        self.infected = infected(self.size, 0.01, self.rng)
        self.clock[self.infected] = periods.latent_time(self.latent_period_base[self.infected])
        self.bact_load[self.infected] = get_load(self.infection_counter[self.infected])
        self.infection_counter[self.infected] = 1


    def tick(self):
        transition = np.logical_not(self.clock.astype(np.bool_))
        new_s = self.diseased & ~self.infected & transition
        new_d = self.infected & self.diseased & transition
        new_id = self.infected & ~self.diseased & transition
        new_i = get_new_infections(~self.infected, self.ages, self.bact_load, self.rng)

        self.clock[new_i] = periods.latent_time(self.latent_period_base[new_i])
        self.clock[new_id] = periods.id_time(
            self.ID_period_base[new_id], self.infection_counter[new_id]
        )
        self.clock[new_d] = periods.d_time(
            self.D_period_base[new_d], self.infection_counter[new_d], self.ages[new_d]
        )

        self.diseased = self.diseased & ~new_s | new_d
        self.infected = self.infected & ~new_d | new_i

        # bacterial load
        self.bact_load[new_d] = 0
        self.bact_load[new_i] = get_load(self.infection_counter[new_i])

        # housekeeping
        self.infection_counter[new_i] += 1
        self.clock -= -1
        self.ages += 1

        # Death
        dead = (self.rng.uniform(0, 1, self.size) < p.bgrd_death_rate) | (self.ages > p.MAX_AGE)
        self.clock[dead] = 0
        self.diseased[dead] = 0
        self.infected[dead] = 0
        self.infection_counter[dead] = 0
        self.ages[dead] = 0

    def prevalence(age_bounds):
        i = self.age <= age_bounds[0] | self.age > age_bounds[1]
        return np.count_nonzero(self.diseased[i]) / np.count_nonzero(i)

    def infection(age_bounds):
        i = self.age <= age_bounds[0] | self.age > age_bounds[1]
        return np.count_nonzero(self.infected[i]) / np.count_nonzero(i)

    def group_size(age_bounds):
        return np.count_nonzero(
            i = self.age <= age_bounds[0] | self.age > age_bounds[1]
        )

    def yearly_high_infection_count(threshold):
        """Return propotion of individuals with no of infection
        greater than threshold, per year group
        """
        i = self.infection_counter > threshold

        # Cast weights to integer to be able to count
        a, _ = np.histogram(self.age, bins=p.MAX_AGE, weights=i.astype(int))
        return a / self.size
