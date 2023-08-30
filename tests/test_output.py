from pathlib import Path
from ntdmc_trachoma.simulation import Simulation

SCENARIO_PATH = Path("tests/input/scenario.json")
PARAMETERS_PATH = Path("tests/input/parameters.json")
BETA_PATH = Path("tests/input/k_values.txt")

with BETA_PATH.open() as f:
    betavals = [
        float(line) for line in f
    ]
sim = Simulation(PARAMETERS_PATH)
sim.simulate(SCENARIO_PATH, betavals, record=False)
sim.output.write()
