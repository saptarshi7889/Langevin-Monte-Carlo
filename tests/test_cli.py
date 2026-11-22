"""Executable-level regression tests, using only the Python standard library."""

import csv
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
BUILD = ROOT / os.environ.get("LANGEVIN_BUILD_DIR", "build")
MODELS = ("incommensurate", "commensurate")
BASE_ARGUMENTS = ("L=3", "Beta=4", "alpha=0.7", "s=0.5", "dt=0.01", "burn_in=2", "sample_every=2", "seed=73")
FULL_STEPS = 14
SPLIT_STEPS = 7
EXPECTED_SAMPLES = 6


class SimulationTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.directory = Path(self.temporary.name)

    def run_program(self, executable, *arguments, success=True):
        result = subprocess.run([str(BUILD / executable), *arguments], capture_output=True, text=True, timeout=30)
        if success:
            self.assertEqual(result.returncode, 0, result.stderr)
        else:
            self.assertNotEqual(result.returncode, 0, result.stdout)
        return result

    def test_public_headers_and_example(self):
        self.run_program("test-headers", str(self.directory))
        self.run_program("library-example", str(self.directory / "example-checkpoint.txt"))

    def test_exact_restart(self):
        for model in MODELS:
            with self.subTest(model=model):
                full, split, resumed = (self.directory / f"{model}-{tag}" for tag in ("full", "split", "resumed"))
                physics = BASE_ARGUMENTS + (("g=0.3",) if model == "commensurate" else ())
                self.run_program(f"initial-{model}", *physics, f"steps={FULL_STEPS}", f"output={full}")
                self.run_program(f"initial-{model}", *physics, f"steps={SPLIT_STEPS}", f"output={split}")
                self.run_program(f"continue-{model}", f"checkpoint={split / 'checkpoint.txt'}", f"steps={FULL_STEPS}", f"output={resumed}")
                for filename in ("checkpoint.txt", "configuration.dat", "roughness.csv", "summary.txt", "kernel.csv"):
                    self.assertEqual((full / filename).read_bytes(), (resumed / filename).read_bytes(), filename)
                def measurements(folder):
                    with (folder / "measurements.csv").open() as handle:
                        return list(csv.DictReader(handle))
                self.assertEqual(measurements(full), measurements(split) + measurements(resumed))
                self.assertIn(f"samples {EXPECTED_SAMPLES}\n", (full / "summary.txt").read_text())
                wrong_model = "commensurate" if model == "incommensurate" else "incommensurate"
                self.run_program(f"continue-{wrong_model}", f"checkpoint={split / 'checkpoint.txt'}", f"output={self.directory / 'wrong-model'}", success=False)
                self.run_program(f"continue-{model}", f"checkpoint={split / 'checkpoint.txt'}", "s=1.5", f"output={self.directory / 'wrong-parameter'}", success=False)

    def test_cli_rejects_invalid_input(self):
        bad_arguments = ("Beta=3", "L=1", "s=0", "s=2", "dt=nan", "dt=-1", "alpha=-1", "sample_every=0", "seed=-1", "L=4x", "steps=0", "unknown=1", "g=0.2")
        for argument in bad_arguments:
            with self.subTest(argument=argument):
                base = tuple(value for value in BASE_ARGUMENTS if value.split("=", 1)[0] != argument.split("=", 1)[0])
                self.run_program("initial-incommensurate", *base, argument, f"output={self.directory / 'invalid'}", success=False)
        self.run_program("initial-incommensurate", "L=3", "L=4", success=False)
        for model in MODELS:
            self.run_program(f"continue-{model}", success=False)
            self.run_program(f"initial-{model}", "checkpoint=missing", success=False)
            self.run_program(f"initial-{model}", "--help")

    def test_no_overwrite_and_no_samples(self):
        output = self.directory / "short"
        self.run_program("initial-incommensurate", *BASE_ARGUMENTS, "steps=1", f"output={output}")
        before = (output / "checkpoint.txt").read_bytes()
        self.assertFalse((output / "roughness.csv").exists())
        self.assertIn("unavailable", (output / "summary.txt").read_text())
        self.run_program("initial-incommensurate", *BASE_ARGUMENTS, "steps=1", f"output={output}", success=False)
        self.assertEqual(before, (output / "checkpoint.txt").read_bytes())

    def test_truncated_checkpoint_rejected(self):
        checkpoint = self.directory / "bad.txt"
        checkpoint.write_text("LANGEVIN_CHECKPOINT 1\n0\n16 16")
        self.run_program("continue-incommensurate", f"checkpoint={checkpoint}", f"output={self.directory / 'bad-output'}", success=False)


if __name__ == "__main__":
    unittest.main()
