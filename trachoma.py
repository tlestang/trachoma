import numpy as np
import init
import parameters as p
from state import Population


rng = np.random.default_rng()
ages = init.ages(p.POP_SIZE, p.MAX_AGE, p.MEAN_AGE, rng)

pop = Population(ages)

children = (1, 9)

fchildren = open("data_children.csv", "w")
f = open("data.csv", "w")

for t in range(25 * 52):
    data_children = f"{2*t},{pop.prevalence(children)},{pop.infection(children)}\n"
    fchildren.write(data_children)
    data_all = f"{2*t},{pop.prevalence((0, p.MAX_AGE))},{pop.infection((0, p.MAX_AGE))}\n"
    f.write(data_all)
    for tt in range(2):
        pop.tick()

fchildren.close()
f.close()


    # def compose(f, g):
#     return lambda x: f(g(x))

# def identity(f)

# from functools import reduce

# reduce(compose, [f]* n)
