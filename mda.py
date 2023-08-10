import numpy as np


class MDA:
    def __init__(self, coverage, efficacy, rho):
        self.a = (1. / rho) - 1
        self.b = self.a * coverage
        self.rng = np.random.default_rng()
        self.count = 0
        self.efficacy = efficacy

    def distribute(self, popsize):
        p = (
            self.b + self.treatment_count
        ) / (self.a + self.count)
        print(f"p = {p}")
        return (self.rng.uniform(size=popsize) < p)

    def __call__(self, pop):
        t = self.distribute(pop.size)
        ntreated = np.count_nonzero(t)
        cured = np.zeros(pop.size, dtype=np.bool_)
        cured[t] = self.rng.uniform(size=ntreated) < self.efficacy
        inf = pop.inf  # Calls the 'inf' property on Population instance
        inf[cured] = False
        pop.inf = inf
        pop.bact_load[cured] = 0

        attrname = f"mda_{id(self)}_treatment_count"
        try:
            treatment_count = getattr(pop, attrname)
            treatment_count[t] += 1
        except AttributeError:
            setattr(pop, attrname, t.astype(np.int32))

        self.count += 1

        return ntreated, np.count_nonzero(cured)


def start_sim(pop):
    pass


def end_sim(pop):
    pass
