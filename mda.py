import numpy as np


class MDA:
    def __init__(self, n, coverage, efficacy, rho):
        self.n = n
        self.treatment_count = np.zeros(n, dtype=np.int32)
        self.a = (1. / rho) - 1
        self.b = self.a * coverage
        self.rng = np.random.default_rng()
        self.count = 0
        self.efficacy = efficacy

    def distribute(self):
        p = (
            self.b + self.treatment_count
        ) / (self.a + self.count)
        print(f"p = {p}")
        return (self.rng.uniform(size=self.n) < p)

    def __call__(self, pop):
        t = self.distribute()
        ntreated = np.count_nonzero(t)
        cured = np.zeros(self.n, dtype=np.bool_)
        cured[t] = self.rng.uniform(size=ntreated) < self.efficacy
        inf = pop.inf  # Calls the 'inf' property on Population instance
        inf[cured] = False
        pop.inf = inf
        pop.bact_load[cured] = 0

        self.count += 1
        self.treatment_count[t] += 1

        return ntreated, np.count_nonzero(cured)
