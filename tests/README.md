# Numerical verification

Run `make test EIGEN_INCLUDE=/path/to/eigen3` from the repository root. No test
framework or Python packages are required; CLI tests use `unittest` and temporary
output directories. `PYTHON=/path/to/python3` selects the test interpreter.

The numerical test executables include the public headers directly, without
including or renaming any executable's `main`. The CLI suite executes all four
thin command-line wrappers. A separate API executable links two translation units,
includes both models together, and verifies initial runs and checkpoint
continuation through the public classes. The standalone library example is also
built and executed by the suite.

- The Ohmic image sum is checked against its rational expression, including its
  prefactor and half-open cutoff; s=0.5 and s=1.5 are checked with an independent
  beta-function recurrence. Zero coupling must produce zero coefficients.
- Periodic Laplacians are checked against analytical Fourier eigenvalues.
- Rectangular-lattice roughness is checked against an analytical Fourier-mode
  result in both directions.
- Both interaction forces are checked against central finite differences of the
  lattice action. The commensurate physical-g normalization is checked separately.
- A constant field with nonzero phase checks the distinct order-parameter
  conventions and ensures measurements do not modify the evolving field.
- A seeded harmonic simulation checks the variance of a nonzero spatial Fourier
  mode against pi/16 (absolute tolerance 0.02, including sampling uncertainty and
  finite-step error). It does not test the diffusing zero mode or production
  parameter equilibration.
- CLI regressions compare complete checkpoints, final configurations, roughness,
  summaries, kernels and combined measurement rows after split versus uninterrupted
  trajectories. The commensurate restart case includes nonzero g, and both cases
  use rectangular lattices. Invalid parameters, duplicate arguments, truncated or
  wrong-model checkpoints, changed restart physics, and output overwrites fail.
- The public API is tested for zero-field initialization, kernel access, optional
  scheduled measurements, read-only state access, exact checkpoint continuation,
  rejection of backward integration and invalid parameters, and model separation.
  Repeated header inclusion and multi-translation-unit linking check header-only
  use without duplicate symbols.

For address and undefined-behavior checks, use a separate build directory:

```bash
make -j2 test BUILD=build-sanitize EIGEN_INCLUDE=/path/to/eigen3 \
  CXXFLAGS='-O1 -g -std=c++17 -Wall -Wextra -Wpedantic -fsanitize=address,undefined -fno-omit-frame-pointer'
```

Local verification used GCC 15.2 and an existing Eigen installation. Optimized and
sanitizer suites passed. Leak detection was disabled for the local sanitizer run
because LeakSanitizer cannot operate under the execution environment's process
tracing; no leak-check result is claimed. Clang and GitHub-hosted execution are
configured in CI but were not run locally.

These checks do not reproduce the published phase diagrams or establish
production equilibration, finite-size scaling, or time-step convergence.
