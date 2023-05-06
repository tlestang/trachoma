import numpy as np


class Population:

    def __init__(self, ages):
        self.size = len(ages)
        self.ages = ages
        self.infected = np.zeros(self.size, dtype=np.bool_)
        self.diseased = np.zeros(self.size, dtype=np.bool_)
        self.clock = np.zeros(self.size) - 1
        self.infection_counter = np.zeros(self.size)
        self.bact_load = np.zeros(self.size)

        self.rng = np.random.default_rng()

        # Baseline disease periods
        self.latent_period_base = np.zeros(p.POP_SIZE) + 2
        self.ID_period_base=self.rng.poisson(lam=p.AV_ID_DURATION, size=p.POP_SIZE)
        self.D_period_base=self.rng.poisson(lam=p.AV_D_DURATION, size=p.POP_SIZE)


    def tick(self):
        transition = np.logical_not(clock.astype(np.bool_))
        new_s = self.diseased & ~self.infected & self.transition
        new_d = self.infected & self.diseased & self.transition
        new_id = self.infected & ~self.diseased & self.transition
        new_i = get_new_infections(~self.infected, self.ages, self.bact_load, rng)

        self.clock[new_i] = periods.latent_time(self.latent_period_base[new_i])
        self.clock[new_id] = periods.id_time(
            self.ID_period_base[new_id], self.infection_counter[new_id]
        )
        self.clock[new_d] = periods.d_time(
            self.D_period_base[new_d], self.infection_counter[new_d], ages[new_d]
        )

        self.diseased = self.diseased & ~new_s | new_d
        self.infected = self.infected & ~new_d | new_i

        # bacterial load
        self.bact_load[new_d] = 0
        self.bact_load[new_i] = get_load(infection_counter[new_i])

        # housekeeping
        self.infection_counter[new_i] += 1
        self.clock -= -1
        self.age += 1