import numpy as np
import init
import parameters as p
from state import Population

rng = np.random.default_rng()
ages = init.ages(p.POP_SIZE, p.MAX_AGE, p.MEAN_AGE, rng)
pop = Population(ages)

children = (52, 10 * 52)

nt = 100
data = np.zeros((nt, 2))
for i in range(nt):
    pop.tick()
    data[i, 0] = pop.infection(children)
    data[i, 1] = pop.prevalence(children)

with open("validate.csv", "w") as f:
    np.savetxt(f, data, delimiter=",", fmt="%.3f")
    
