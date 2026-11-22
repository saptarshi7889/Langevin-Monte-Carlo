"""Run both models and resume a checkpoint without installing NumPy."""

from dataclasses import replace
from pathlib import Path
from tempfile import TemporaryDirectory

from langevin_monte_carlo import Parameters, commensurate, incommensurate

SPACE_SIZE = 3
TIME_SIZE = 4
BURN_IN_STEPS = 2
SAMPLE_INTERVAL = 2
CHECKPOINT_STEP = 7
FINAL_STEP = 14
COMMENSURATE_G = 0.2


def main():
    parameters = Parameters(L=SPACE_SIZE, Beta=TIME_SIZE, burn_in=BURN_IN_STEPS,
                            sample_every=SAMPLE_INTERVAL)
    initial = incommensurate.Simulation(parameters)
    initial.advance_to(CHECKPOINT_STEP)
    with TemporaryDirectory(prefix="langevin-python-example-") as directory:
        checkpoint = Path(directory) / "checkpoint.txt"
        initial.save_checkpoint(checkpoint)
        continued = incommensurate.Simulation.from_checkpoint(checkpoint)
        continued.advance_to(FINAL_STEP)
        state = continued.state
        print("Incommensurate field shape:", state.phi.shape)
        print("Cumulative mean order:", state.mean_order_parameter)

    simulation = commensurate.Simulation(replace(parameters, g=COMMENSURATE_G))
    while simulation.current_step < FINAL_STEP:
        sample = simulation.step()
        if sample is not None:
            print("Commensurate sample:", sample)


if __name__ == "__main__":
    main()
