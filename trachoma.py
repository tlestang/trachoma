import numpy as np
import init
import parameters as p
from state import Population


rng = np.random.default_rng()
ages = init.ages(p.POP_SIZE, p.MAX_AGE, p.MEAN_AGE, rng)

pop = Population(ages)

for t in range(0,10000):
    pop.tick()


    # def compose(f, g):
#     return lambda x: f(g(x))

# def identity(f)

# from functools import reduce

# reduce(compose, [f]* n)
