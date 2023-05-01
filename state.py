import np as np


class Population:
    ""
    def __init__(self, ages):
        self.size = len(ages)
        self.ages = ages
        self.infected = np.zeros(self.size, dtype=np.bool_)
        self.diseased = np.zeros(self.size, dtype=np.bool_)
        self.clock = np.zeros(self.size) - 1
        self.infection_counter = np.zeros(self.size)
        self.bact_load = np.zeros(self.size)


    def tick(self):
            transition = np.logical_not(clock.astype(np.bool_))
            new_s = self.diseased & ~self.infected & self.transition
            new_d = self.infected & self.diseased & self.transition
            new_id = self.infected & ~self.diseased & self.transition
            new_i = get_new_infections([new_i] = periods.latent_time(latent_period_base[new_i])
    clock[new_id] = periods.id_time(ID_period_base[new_id], infection_counter[new_id])
    clock[new_d] = periods.d_time(D_period_base[new_d], infection_counter[new_d], ages[new_d])

    diseased = diseased & ~new_s | new_d
    infected = infected & ~new_d | new_i

    bact_load[new_d] = 0
    bact_load[new_i] = get_load(infection_counter[new_i])

    # housekeeping
    infection_counter[new_i] += 1
    clock -= -1
