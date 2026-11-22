"""Python-library regressions against the same native C++ simulations."""

from concurrent.futures import ThreadPoolExecutor
import csv
from dataclasses import FrozenInstanceError, replace
import gc
import math
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

from langevin_monte_carlo import Parameters, commensurate, incommensurate

ROOT = Path(__file__).resolve().parents[1]
BUILD = ROOT / os.environ.get("LANGEVIN_BUILD_DIR", "build")
MODELS = (("incommensurate", incommensurate), ("commensurate", commensurate))
PARAMETERS = Parameters(L=3, Beta=4, alpha=0.7, s=0.5, dt=0.01,
                        burn_in=2, sample_every=2, seed=73)
SPLIT_STEP = 7
FINAL_STEP = 14
EXPECTED_SAMPLES = 6
OBSERVATION_STEP = 4
THREAD_STEPS = 100


class PythonLibraryTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.directory = Path(self.temporary.name)

    @staticmethod
    def parameters(model):
        return replace(PARAMETERS, g=0.3) if model == "commensurate" else PARAMETERS

    def test_python_matches_cli_and_bidirectional_checkpoints(self):
        for name, model in MODELS:
            with self.subTest(model=name):
                p = self.parameters(name)
                simulation = model.Simulation(p)
                simulation.advance_to(FINAL_STEP)
                python_checkpoint = self.directory / f"{name}-python.txt"
                simulation.save_checkpoint(python_checkpoint)
                args = [f"{key}={value}" for key, value in p.__dict__.items()]
                output = self.directory / f"{name}-cli"
                subprocess.run([str(BUILD / f"initial-{name}"), *args,
                                f"steps={FINAL_STEP}", f"output={output}"],
                               capture_output=True, check=True, timeout=30)
                self.assertEqual(python_checkpoint.read_bytes(), (output / "checkpoint.txt").read_bytes())
                state = simulation.state
                field = [[float(value) for value in row.split()]
                         for row in (output / "configuration.dat").read_text().splitlines()]
                self.assertEqual(state.phi.tolist(), field)
                with (output / "roughness.csv").open() as handle:
                    for row in csv.DictReader(handle):
                        sums = state.rough_x if row["axis"] == "space" else state.rough_tau
                        self.assertEqual(sums[int(row["separation"]) - 1] / state.samples,
                                         float(row["mean_squared_difference"]))
                from_cli = model.Simulation.from_checkpoint(output / "checkpoint.txt")
                from_cli.advance_to(FINAL_STEP + 1)
                continued_output = self.directory / f"{name}-continued"
                subprocess.run([str(BUILD / f"continue-{name}"),
                                f"checkpoint={python_checkpoint}", f"steps={FINAL_STEP + 1}",
                                f"output={continued_output}"], capture_output=True, check=True, timeout=30)
                from_cli.save_checkpoint(self.directory / f"{name}-after.txt")
                self.assertEqual((self.directory / f"{name}-after.txt").read_bytes(),
                                 (continued_output / "checkpoint.txt").read_bytes())

    def test_exact_restart(self):
        for name, model in MODELS:
            with self.subTest(model=name):
                p = self.parameters(name)
                initial = model.Simulation(p)
                initial.advance_to(SPLIT_STEP)
                checkpoint = self.directory / f"{name}.txt"
                initial.save_checkpoint(checkpoint)
                continued = model.Simulation.from_checkpoint(checkpoint)
                continued.advance_to(FINAL_STEP)
                full = model.Simulation(p)
                full.advance_to(FINAL_STEP)
                continued.save_checkpoint(self.directory / f"{name}-continued.txt")
                full.save_checkpoint(self.directory / f"{name}-full.txt")
                self.assertEqual((self.directory / f"{name}-continued.txt").read_bytes(),
                                 (self.directory / f"{name}-full.txt").read_bytes())
                self.assertEqual(continued.samples, EXPECTED_SAMPLES)
                self.assertEqual(continued.state.parameters, p)

    def test_schedule_and_owned_snapshots(self):
        for name, model in MODELS:
            with self.subTest(model=name):
                simulation = model.Simulation(self.parameters(name))
                original = simulation.state
                self.assertTrue(original.phi.readonly)
                self.assertEqual(original.phi.shape, (PARAMETERS.Beta, PARAMETERS.L))
                self.assertTrue(all(value == 0 for row in original.phi.tolist() for value in row))
                self.assertIsNone(original.mean_order_parameter)
                with self.assertRaises(TypeError):
                    original.phi[0, 0] = 1
                self.assertIsNone(simulation.step())
                simulation.advance_to(OBSERVATION_STEP - 1)
                observation = simulation.step()
                self.assertEqual(observation.step, OBSERVATION_STEP)
                self.assertEqual(observation.time, OBSERVATION_STEP * PARAMETERS.dt)
                self.assertTrue(math.isfinite(observation.order_parameter))
                self.assertEqual(simulation.current_step, OBSERVATION_STEP)
                self.assertEqual(simulation.time, observation.time)
                self.assertEqual(simulation.samples, 1)
                self.assertEqual(simulation.state.mean_order_parameter, observation.order_parameter)
                current = simulation.state
                self.assertTrue(any(value != 0 for row in current.phi.tolist() for value in row))
                self.assertEqual(len(current.rough_x), PARAMETERS.L - 1)
                self.assertEqual(len(current.rough_tau), PARAMETERS.Beta - 1)
                self.assertTrue(all(value >= 0 for value in current.rough_x))
                self.assertTrue(simulation.coefficients.readonly)
                self.assertEqual(len(simulation.coefficients), PARAMETERS.Beta - 1)
                del simulation
                gc.collect()
                self.assertEqual(original.phi[0, 0], 0)
                self.assertTrue(math.isfinite(current.phi[0, 0]))

    def test_validation_and_exceptions(self):
        for name, model in MODELS:
            with self.subTest(model=name):
                with self.assertRaises(TypeError):
                    model.Simulation({"L": 3})
                with self.assertRaises(ValueError):
                    model.Simulation(replace(PARAMETERS, Beta=3))
                with self.assertRaises(ValueError):
                    model.Simulation(replace(PARAMETERS, dt=float("nan")))
                with self.assertRaises(OverflowError):
                    model.Simulation(replace(PARAMETERS, seed=-1))
                with self.assertRaises(OverflowError):
                    model.Simulation(replace(PARAMETERS, L=2**63))
                with self.assertRaises(TypeError):
                    model.Simulation(replace(PARAMETERS, L=3.5))
                simulation = model.Simulation(self.parameters(name))
                simulation.advance_to(SPLIT_STEP)
                with self.assertRaises(ValueError):
                    simulation.advance_to(SPLIT_STEP - 1)
                with self.assertRaises(TypeError):
                    simulation.advance_to(1.5)
                with self.assertRaises(OverflowError):
                    simulation.advance_to(2**63)
                simulation.advance_to(SPLIT_STEP)
                path = self.directory / f"{name}.txt"
                simulation.save_checkpoint(path)
                wrong_model = commensurate if name == "incommensurate" else incommensurate
                with self.assertRaises(RuntimeError):
                    wrong_model.Simulation.from_checkpoint(path)
                with self.assertRaises(RuntimeError):
                    model.Simulation.from_checkpoint(self.directory / "missing.txt")
                with self.assertRaises(RuntimeError):
                    simulation.save_checkpoint(self.directory / "missing" / "checkpoint.txt")
                path.write_text("truncated")
                with self.assertRaises(RuntimeError):
                    model.Simulation.from_checkpoint(path)
                with self.assertRaises(FrozenInstanceError):
                    PARAMETERS.g = 1
        with self.assertRaises(ValueError):
            incommensurate.Simulation(replace(PARAMETERS, g=0.3))

    def test_numpy_interoperability_when_available(self):
        try:
            import numpy as np
        except ImportError:
            self.skipTest("NumPy is optional")
        for _, model in MODELS:
            simulation = model.Simulation(PARAMETERS)
            simulation.advance_to(OBSERVATION_STEP)
            state = simulation.state
            field = np.asarray(state.phi)
            self.assertEqual(field.shape, (PARAMETERS.Beta, PARAMETERS.L))
            self.assertFalse(field.flags.writeable)
            self.assertEqual(field.tolist(), state.phi.tolist())
            self.assertTrue(np.isfinite(field).all())
            self.assertTrue(np.isfinite(np.asarray(simulation.coefficients)).all())
            snapshot = field.copy()
            simulation.advance_to(FINAL_STEP)
            np.testing.assert_array_equal(field, snapshot)

    def test_independent_instances_from_threads(self):
        def evolve(seed):
            simulation = commensurate.Simulation(replace(PARAMETERS, seed=seed, g=0.2))
            simulation.advance_to(THREAD_STEPS)
            return simulation.state.phi.tobytes()
        seeds = (73, 74)
        sequential = [evolve(seed) for seed in seeds]
        with ThreadPoolExecutor(max_workers=len(seeds)) as executor:
            parallel = list(executor.map(evolve, seeds))
        self.assertEqual(sequential, parallel)
        self.assertNotEqual(*parallel)


if __name__ == "__main__":
    unittest.main()
