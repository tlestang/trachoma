from pathlib import Path

from ntdmc_trachoma.simulation import Simulation

SCENARIO_PATH = Path("scenario.json")
BETA_PATH = Path("k_values.txt")
PARAMETERS_PATH = Path("parameters.json")


def main():

    sim = Simulation(PARAMETERS_PATH)

    with BETA_PATH.open() as f:
        betavals = [
            float(line) for line in f
        ]
    sim.simulate(SCENARIO_PATH, betavals, record=False)


if __name__ == "__main__":
    main()
