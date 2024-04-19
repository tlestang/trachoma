from pathlib import Path
from ntdmc_trachoma.simulation import Simulation

param_path = Path("parameters.json")
sim = Simulation(param_path)

betavals = [0.2, 0.35, 0.4]
scenario_path = Path('scenario.json')
sim.simulate(scenario_path, betavals, record=True)

sim.output.write()
