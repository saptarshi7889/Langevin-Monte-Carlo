# Langevin Monte Carlo for a dissipative XXZ spin chain

Research code associated with Saptarshi Majumdar's doctoral work on quantum spin chains coupled to local dissipative baths.

![Stylized three-dimensional schematic of an XXZ chain: six blue spins joined by exchange bonds, each coupled to a separate coral bath of harmonic oscillators at zero temperature.](assets/dissipative-xxz-schematic.png)

Each spin couples to its own independent oscillator bath at zero temperature; $\alpha$ measures the spin–bath coupling strength. The colored bonds distinguish two nearest-neighbor exchange terms:

- **Teal — transverse exchange:**

  $J_{xy}\left(S_i^x S_{i+1}^x + S_i^y S_{i+1}^y\right)$.


- **Purple — longitudinal exchange:**

  $J_z S_i^z S_{i+1}^z.$

Their different heights distinguish the two terms acting on the same neighboring sites, rather than separate spatial paths. Spin directions and oscillator arrangements are illustrative and do not specify a phase. AI-generated illustration inspired by Figure 1 of [Majumdar et al., Physical Review B 107, 165113 (2023)](https://doi.org/10.1103/PhysRevB.107.165113).

## Research context

We investigate the zero-temperature phase diagram of a one-dimensional XXZ spin chain in which each spin couples to an independent bath of harmonic oscillators: a many-body extension of the spin-boson model. Bosonization and a path-integral formulation map the quantum problem onto a classical field theory in one spatial dimension and imaginary time. Integrating out the baths produces a cosine interaction that is local in space and long-ranged along imaginary time.

The associated papers explore how dissipation changes ordering and transport:

- At finite, incommensurate magnetization, coupling to the baths can drive a transition from a Luttinger liquid to a gapless spin density wave. For subohmic baths, the ordered phase has vanishing DC conductivity, interpreted as bath-induced localization [2].
- At zero magnetization, the phase diagram contains a Luttinger liquid and a gapped antiferromagnetic phase. Subohmic dissipation can drive the ordering transition; for superohmic baths, interactions within the chain drive the transition [1].
- For Ohmic baths, the dissipative phase has vanishing spin stiffness and finite resistivity, distinct from the vanishing DC conductivity found for subohmic baths [3].

## Numerical method in the published work

Appendix C of [*Bath-induced phase transition in a Luttinger liquid*, Physical Review B 107, 165113 (2023)](https://doi.org/10.1103/PhysRevB.107.165113) describes Langevin sampling of the discretized bosonized action. The field $\phi_{ij}$ lives on an $L\times\beta$ lattice, with $i$ labeling space and $j$ labeling imaginary time, and periodic boundary conditions in both directions. Starting from a flat configuration, $\phi_{ij}(0)=0$, it evolves in an auxiliary sampling time $t$:

$$
\frac{d\phi_{ij}}{dt}=-\frac{\delta S}{\delta\phi_{ij}}+\eta_{ij}(t),
\qquad
\langle\eta_{ij}(t)\rangle=0,
$$

$$
\langle\eta_{ij}(t)\eta_{i'j'}(t')\rangle
=2\delta_{ii'}\delta_{jj'}\delta(t-t').
$$

After equilibration, configurations sample the weight $e^{-S}$. The sampling time $t$ is distinct from both imaginary time $\tau$ and physical real time.

For the Ohmic problem, the lattice equation is

$$
\begin{aligned}
\frac{d\phi_{ij}}{dt}
&=\frac{u}{\pi K}\left(\phi_{i+1,j}+\phi_{i-1,j}-2\phi_{ij}\right)\\
&\quad+\frac{1}{\pi uK}\left(\phi_{i,j+1}+\phi_{i,j-1}-2\phi_{ij}\right)\\
&\quad+\frac{\alpha}{\pi^2}\sum_{j'\ne j}D_\beta(|j-j'|)
\sin\left[2(\phi_{ij'}-\phi_{ij})\right]+\eta_{ij}(t).
\end{aligned}
$$

Here $u$ is the Luttinger velocity and $K$ the Luttinger parameter. The equal-time term is omitted because it produces no dissipative force. Appendix C discretizes the long-range kernel as $D_\infty(r)=1/(r^2-1/4)$ for nonzero integer separation $r$, and implements periodicity at finite $\beta$ through the image sum

$$
D_\beta(r)=\sum_{k=-\beta/2}^{\beta/2}
\frac{1}{(r+k\beta)^2-\tfrac14},
\qquad r=|j-j'|>0.
$$

The published simulations use a **second-order stochastic Runge–Kutta algorithm for white noise**, following R. L. Honeycutt, [*Stochastic Runge-Kutta algorithms. I. White noise*, Physical Review A **45**, 600 (1992)](https://doi.org/10.1103/PhysRevA.45.600), reference **34** of the paper. The Langevin step is $\Delta t=0.05$: drift increments scale as $\Delta t$, while noise increments scale as $\sqrt{\Delta t}$. Equilibration times are benchmarked against the analytically solvable harmonic approximation,

$$
\sin\left[2(\phi_{ij'}-\phi_{ij})\right]
\longrightarrow 2(\phi_{ij'}-\phi_{ij}).
$$

## Build and test

The simulation package requires a C++17 compiler, Make, and Eigen headers (Eigen 3.3 or later). The CLI regression tests use Python 3's standard library; NumPy and Matplotlib are only needed for the optional plotting script. MPI is not required. No dependencies are downloaded by the Makefile.

```bash
make -j2 EIGEN_INCLUDE=/path/to/eigen3
make test EIGEN_INCLUDE=/path/to/eigen3
```

`EIGEN_INCLUDE` must contain `Eigen/`; it defaults to `/usr/include/eigen3`. An optional [project Conda environment](environment.yml) specifies the build, test, and plotting dependencies. Creating that environment installs packages; review it before running `conda env create -f environment.yml`. After activation, build with `make EIGEN_INCLUDE="$CONDA_PREFIX/include/eigen3"`.

The [MIT license](LICENSE) covers this repository's code. [CITATION.cff](CITATION.cff) identifies the software; please also cite the relevant papers below. This package supports reproducible small examples and numerical tests, but does not contain a full reproduction of all published phase diagrams.

## Package layout

```text
src/
  incommensurate/  # Initial and continued finite-magnetization simulations
  commensurate/    # Initial and continued zero-magnetization simulations
include/langevin/
  incommensurate.hpp  # Header-only API: initial and continued incommensurate runs
  commensurate.hpp    # Header-only API: initial and continued commensurate runs
  simulation.hpp     # Shared state, integration, measurements and checkpoints
  numerics.hpp       # Shared lattice and interaction functions
  kernel.hpp         # Generalized dissipative kernel
  cli.hpp            # Optional command-line adapter
python/            # Native CPython bindings and model-specific Python modules
examples/          # Standalone C++ and Python library usage
tests/             # Numerical tests and executable-level restart tests
scripts/           # Optional plotting utilities
figures/           # Example kernel figure
```

| Model | Initial executable | Continuation executable | Order parameter |
| --- | --- | --- | --- |
| Incommensurate | `build/initial-incommensurate` | `build/continue-incommensurate` | $\langle\cos(\phi-\overline\phi)\rangle$ |
| Commensurate | `build/initial-commensurate` | `build/continue-commensurate` | $\langle\cos(2\phi)\rangle$ |

Initial runs start from $\phi=0$. Continuation restores a complete checkpoint, including the unshifted field, absolute step, parameters, accumulated measurements, random-generator state, and the Gaussian distribution's cached state. The original source filenames are preserved under their model directories, including the spelling `continuted` in the incommensurate continuation source.

## Using the header-only C++ libraries

The two related pairs of programs are available as **two C++17 header-only libraries**:

- `#include <langevin/incommensurate.hpp>` provides `langevin::incommensurate::Simulation`.
- `#include <langevin/commensurate.hpp>` provides `langevin::commensurate::Simulation`.

Each `Simulation` supports both zero-field initialization and checkpoint continuation. You can include both model headers in one program and in multiple translation units. All library definitions are inline or templates; there is no compiled library to link, no global namespace import, and no dependency on the four executable source files. Add this repository's `include/` and your Eigen header directory to the compiler's include search path. The four CLI executables now use these same libraries.

```cpp
#include <langevin/incommensurate.hpp>
#include <langevin/commensurate.hpp>

int main()
{
    constexpr int INITIAL_FINAL_STEP = 2000;
    constexpr int CONTINUED_FINAL_STEP = 4000;
    constexpr double COMMENSURATE_G = 0.2;

    langevin::Parameters parameters;
    langevin::incommensurate::Simulation initial(parameters);
    initial.advance_to(INITIAL_FINAL_STEP);
    initial.save_checkpoint("checkpoint.txt");

    auto continued =
        langevin::incommensurate::Simulation::from_checkpoint("checkpoint.txt");
    continued.advance_to(CONTINUED_FINAL_STEP);
    const auto& field = continued.state().phi;

    parameters.g = COMMENSURATE_G;
    langevin::commensurate::Simulation commensurate(parameters);
    if (const auto sample = commensurate.step()) {
        // Use sample->step, sample->time and sample->order_parameter here.
    }
    return field.allFinite() ? 0 : 1;
}
```

| API | Behavior |
| --- | --- |
| `Simulation(parameters)` | Validate parameters and create a zero-field initial simulation |
| `Simulation::from_checkpoint(path)` | Restore the same model, parameters, counters, accumulated measurements and RNG state |
| `step()` | Advance one step; return `std::optional<langevin::Measurement>` only when a measurement is scheduled |
| `advance_to(final_step)` | Advance to an absolute step, automatically accumulating scheduled measurements; an equal target is a no-op and a backward target throws |
| `state()` | Read-only `langevin::State`, including the field, parameters, absolute step, sample count and accumulated observable sums |
| `coefficients()` | Read-only interaction coefficients, including $\alpha/\pi^2$ |
| `save_checkpoint(path)` | Explicitly write or replace a complete checkpoint; the parent directory must already exist |

Initialization and stepping perform no file output. File writes happen only when you request a checkpoint; the optional CLI adapter handles its own output tables. Invalid inputs and I/O failures throw standard C++ exceptions rather than terminating the host application. The field has `Beta` rows and `L` columns. Mean order and roughness are obtained by dividing the corresponding sums in `state()` by `state().samples`, after checking that the sample count is nonzero. Separate simulation objects own their state and random generators; copies reproduce the same trajectory, so use distinct seeds for independent runs.

A [complete example](examples/library_usage.cpp) demonstrates both models and continuation:

```bash
make examples EIGEN_INCLUDE=/path/to/eigen3
build/library-example /tmp/langevin-example-checkpoint.txt
```

To use the same example from another project, compile it with `-std=c++17 -I/path/to/Langevin-Monte-Carlo/include -I/path/to/eigen3`. No source files from `src/` need to be compiled into that project.

## For Python users

The Python package is `langevin_monte_carlo` (distribution name: `langevin-monte-carlo-xxz`). It calls the same C++ header libraries through a native CPython extension; simulation loops are not reimplemented in Python. Use `incommensurate.Simulation` for finite, incommensurate magnetization and `commensurate.Simulation` for the zero-magnetization model with the two additional forces described below. Both classes support initial runs and checkpoint continuation.

### Build or install in your project environment

Requirements are **CPython 3.10 or later**, its development headers, a C++17 compiler, Eigen, setuptools 64 or later, and wheel for building. The core Python API needs no NumPy, pybind11, or other binding-framework runtime. NumPy and Matplotlib are useful for analysis and plotting. The provided [Conda environment](environment.yml) includes the compiler, Eigen, Python, packaging tools and plotting packages. Linux is tested locally; macOS is supported by the build configuration but has not been tested here. Windows/MSVC and free-threaded Python have not been validated. Native binaries must be built for the Python interpreter and platform that will load them.

After reviewing the environment specification, create and activate a dedicated environment if you do not already have one:

```bash
conda env create -f environment.yml
conda activate langevin-monte-carlo
export EIGEN_INCLUDE="$CONDA_PREFIX/include/eigen3"
```

Creating the environment installs packages. These commands are instructions for your own setup; no environment or dependency has been installed automatically. To install from this repository into that environment:

```bash
python -m pip install --no-build-isolation --no-deps .
```

This uses the existing environment's build tools and avoids downloading build dependencies. It is a local source installation, not a download from PyPI; this project has not been published there. There is no need to build the CLI executables for ordinary Python use.

Alternatively, **build and use Python bindings without installing a package**:

```bash
make python-bindings PYTHON=python EIGEN_INCLUDE="$EIGEN_INCLUDE"
PYTHONPATH="$PWD/build/python" python examples/python_usage.py
```

Choose `PYTHON` from the intended environment; its matching development headers must be available. The build writes only into this checkout. For Jupyter, use a kernel from the same environment and, if using the uninstalled build, add the absolute `build/python` directory to `sys.path` before importing. Restart the kernel after rebuilding an already imported native extension.

### Initialize, sample and continue

```python
from dataclasses import replace
from pathlib import Path
from langevin_monte_carlo import Parameters, commensurate, incommensurate

INITIAL_FINAL_STEP = 2000
CONTINUED_FINAL_STEP = 4000
COMMENSURATE_G = 0.2
CHECKPOINT = Path("python-checkpoint.txt")

parameters = Parameters(
    L=16, Beta=16, K=1.0, u=1.0, alpha=1.0, s=0.5,
    dt=0.05, burn_in=1000, sample_every=20, seed=73,
)
simulation = incommensurate.Simulation(parameters)  # Zero field at t=0.
simulation.advance_to(INITIAL_FINAL_STEP)
state = simulation.state
print(state.time, state.samples, state.mean_order_parameter)

simulation.save_checkpoint(CHECKPOINT)
continued = incommensurate.Simulation.from_checkpoint(CHECKPOINT)
continued.advance_to(CONTINUED_FINAL_STEP)

# Parameters are immutable: replace selected values for a new simulation.
afm = commensurate.Simulation(replace(parameters, g=COMMENSURATE_G, seed=74))
observations = []
while afm.current_step < INITIAL_FINAL_STEP:
    observation = afm.step()
    if observation is not None:
        observations.append(observation)
```

`advance_to(n)` targets the **absolute step** $n$, not an additional duration; $t=n\Delta t$. It accumulates scheduled measurements but does not return their individual rows. Use `step()` when you want a time series: it returns a frozen `Measurement(step, time, order_parameter)` on scheduled steps and `None` otherwise. Defaults, parameter ranges, burn-in, measurement intervals and the physical meaning of `g` are the same as in the C++/CLI descriptions. `Parameters` and `Measurement` are frozen dataclasses; `dataclasses.replace` creates altered parameters for a new run.

Checkpoints are **interchangeable between Python and the C++ executables** for the same model and numerical environment. They restore the field, sampling schedule, observable sums and RNG state; there is no parameters argument to `from_checkpoint` that could silently change the saved physics. The original pre-package matrices and filename-based restart files remain unsupported. Checkpoint writes explicitly replace their target; choose a new filename if you want to preserve an older checkpoint. Initialization and stepping otherwise do not write files.

### NumPy arrays and plotting

`simulation.state` returns an **owned, read-only snapshot**, rather than live mutable C++ memory. The field is a native-double `memoryview` with shape `(Beta, L)`: axis 0 is imaginary time and axis 1 is space. Without NumPy, use `state.phi[time_index, space_index]` or `state.phi.tolist()`. With NumPy:

```python
import numpy as np
import matplotlib.pyplot as plt

FIGURE_SIZE = (7, 4)
state = continued.state
field = np.asarray(state.phi)       # Read-only view of the snapshot.
editable_field = field.copy()      # Editing this never changes the simulation.
kernel = np.asarray(continued.coefficients)

fig, ax = plt.subplots(figsize=FIGURE_SIZE)
image = ax.imshow(field, aspect="auto", origin="lower")
ax.set_xlabel("Spatial site")
ax.set_ylabel("Imaginary-time site")
fig.colorbar(image, ax=ax, label=r"$\phi$")
fig.tight_layout()
plt.show()
```

A snapshot remains valid after advancing or deleting its simulation. Each `state` access copies the field and accumulated roughness arrays; use `simulation.current_step`, `simulation.samples` and `simulation.time` to check progress without copying the field. `simulation.coefficients` is a separate read-only snapshot of the kernel, including $\alpha/\pi^2$.

`state.order_sum`, `state.rough_x` and `state.rough_tau` are **cumulative sums**, not averages. `state.mean_order_parameter` is `None` before the first measurement. After checking `state.samples > 0`, obtain mean roughness as `np.asarray(state.rough_x) / state.samples` and `np.asarray(state.rough_tau) / state.samples`. The incommensurate order parameter uses the mean-subtracted field; the commensurate order parameter preserves the absolute phase, matching the C++ libraries.

### Errors, parallelism and verification

Invalid physical parameters raise `ValueError`; invalid Python argument types raise `TypeError`; integers outside the C++ size/step/seed range raise `OverflowError`. Invalid, missing or incompatible checkpoints and native file-operation failures raise `RuntimeError`. These errors propagate as Python exceptions rather than terminating the interpreter.

Long C++ advances release the Python GIL. Independent simulation objects may run in separate threads; give them distinct seeds for independent trajectories. Concurrent access to the **same** simulation raises `RuntimeError` while that object is busy. For process-based parallelism, construct the simulation inside each worker; native simulation handles are not picklable. Do not use the same checkpoint path from multiple workers. Statistical uncertainty still requires accounting for autocorrelation.

Run the Python suite with:

```bash
make test-python PYTHON=python EIGEN_INCLUDE="$EIGEN_INCLUDE"
```

It builds the extension and CLI binaries without installing them, then checks Python/CLI agreement, exact continuation in both directions, measurement timing, input errors, snapshot ownership, independent threaded trajectories, and NumPy conversion when NumPy is available. The Python API has been tested locally with CPython 3.14; GitHub CI is configured to test it separately from the C++ suite. A [complete Python example](examples/python_usage.py) runs with the standard library alone.

## Running a simulation

Arguments use `name=value`; `--help` lists the interface. The defaults are deliberately small examples, rather than the production lattice sizes and durations used in the papers.

| Parameter | Default | Meaning |
| --- | --- | --- |
| `L`, `Beta` | `16`, `16` | Space and imaginary-time extents; `L >= 2`, even `Beta >= 2` |
| `u`, `K` | `1`, `1` | Positive Luttinger velocity and parameter |
| `alpha`, `s` | `1`, `0.5` | Nonnegative bath coupling and exponent $0<s<2$ |
| `g` | `0` | Nonnegative commensurate action coupling; must be zero in the incommensurate model |
| `dt` | `0.05` | Positive Langevin time step |
| `steps` | `2000` | Final **absolute step**, also for continuation; $t=\text{step}\times\Delta t$ |
| `burn_in` | `1000` | Initial steps excluded from measurements |
| `sample_every` | `20` | Step interval between measurements after burn-in |
| `seed` | `1234` | Explicit nonnegative integer seed; use a different seed for each independent trajectory |
| `id` | `0` | Run label; selects `run-ID` when `output` is omitted, but does not change the seed |
| `output` | `run-0` | New output directory; an existing directory is never overwritten |
| `checkpoint` | none | Required input path for continuation executables |

The first measurement is at `burn_in + sample_every`. All scheduling uses integer steps, avoiding floating-point time-loop drift. At zero samples, the summary reports that the order parameter is unavailable; no roughness averages are written.

```bash
build/initial-incommensurate L=16 Beta=16 s=0.5 alpha=1 \
  steps=2000 burn_in=1000 sample_every=20 seed=73 output=runs/inc-73

build/continue-incommensurate checkpoint=runs/inc-73/checkpoint.txt \
  steps=4000 output=runs/inc-73-continued

build/initial-commensurate L=16 Beta=16 s=0.5 alpha=1 g=0.2 \
  steps=2000 burn_in=1000 sample_every=20 seed=74 output=runs/comm-74

build/continue-commensurate checkpoint=runs/comm-74/checkpoint.txt \
  steps=4000 output=runs/comm-74-continued
```

Continuation rejects changes to the saved physical parameters, seed, and sampling schedule, as well as checkpoints for the wrong model. `steps` must exceed the saved step. Split and uninterrupted trajectories agree exactly with the same executable and floating-point environment; cross-platform random-number or bitwise reproducibility is not promised. The final checkpoint is written after a successful run; periodic crash-recovery checkpoints are not implemented.

### Output files

| File | Contents |
| --- | --- |
| `checkpoint.txt` | Versioned complete restart state, parameters and cumulative measurement sums |
| `configuration.dat` | Final **unshifted** field, with `Beta` rows and `L` columns |
| `summary.txt` | Model, run label, final step/time, cumulative sample count and mean order parameter |
| `measurements.csv` | Step, time and order parameter for samples generated in this invocation |
| `roughness.csv` | Cumulative mean-squared field differences, separately for space and imaginary time |
| `kernel.csv` | Stored interaction coefficients, including $\alpha/\pi^2$ |

Continuation keeps cumulative statistics from the checkpoint but writes only its new measurement rows. Concatenate measurement tables with one header if a complete time series is needed. Both models now initialize, update, save and restore roughness accumulators, including on rectangular lattices. Incommensurate observables use the mean-subtracted field; commensurate observables keep the absolute field phase because arbitrary uniform shifts are not a symmetry of that action. Neither model shifts its saved field.

### Independent trajectories in parallel

Each executable is serial. Parallelize independent trajectories with distinct seeds and output directories; `id` alone does not make trajectories independent. For example:

```bash
build/initial-commensurate seed=1001 id=0 output=runs/comm-1001 &
build/initial-commensurate seed=1002 id=1 output=runs/comm-1002 &
wait
```

Use one process per trajectory in a job array on a cluster. The programs do not read MPI ranks, distribute a lattice, or aggregate trajectories. Do not use identical arguments under `mpirun`: identical seeds would duplicate trajectories, and identical output paths would conflict. Estimate uncertainty using independent trajectories or time-series methods that account for autocorrelation; increasing the number of saved samples alone does not establish independence.

## Discretization and kernel

All four programs retain the two-stage additive-noise stochastic predictor–corrector update. A single matrix of independent standard Gaussian values $\xi$ is reused in both stages:

$$
\begin{aligned}
\phi^{\mathrm{pred}}&=\phi+\Delta t\,F(\phi)+\sqrt{2\Delta t}\,\xi,\\
\phi^{\mathrm{new}}&=\phi+\frac{\Delta t}{2}
\left[F(\phi)+F(\phi^{\mathrm{pred}})\right]+\sqrt{2\Delta t}\,\xi.
\end{aligned}
$$

In the source matrices, rows label imaginary time and columns label space, the reverse of the paper's index convention used in the equations above. Spatial and imaginary-time Laplacians use periodic boundary conditions. The direct interaction evaluation costs $O(L\beta^2)$ per step; these programs do not use FFT acceleration.

For even $\beta$, `langevin::kernel()` stores

$$
\mathcal{K}_{\beta,s}(r)=\frac{\alpha}{\pi^2}
\sum_{k=-\beta/2}^{\beta/2-1}
\mathrm{B}\left(|r+k\beta|-\frac{s}{2},s+1\right),
\qquad r=1,\ldots,\beta-1,
$$

where $\mathrm B$ is the Euler beta special function, distinct from the extent $\beta$. For integer $n\ne0$ and $0<s<2$, the analytic continuation of $\mathrm B(n-s/2,s+1)$ is even in $n$. The implementation now reflects negative integer images before evaluating the function, then uses log-gamma values with positive arguments. This avoids relying on standard-library behavior for negative beta-function arguments. The restriction $0<s<2$ ensures that the reflected arguments are positive; larger exponents are not supported by this package.

The original half-open image range is retained. This finite image truncation approximates a periodic image sum, with a small reflection error; no exact finite-temperature kernel is claimed. The stored coefficients already contain $\alpha/\pi^2$. The equal-time term is excluded. For $s=1$,

$$
\mathrm{B}\left(n-\tfrac12,2\right)=\frac{1}{n^2-\tfrac14},
$$

recovering the discretized Ohmic kernel. The [existing $\beta=256$ figure](figures/dissipative-kernel-beta256.png) and [plotting script](scripts/plot_dissipative_kernel.py) show the **printed Appendix C** Ohmic sum, which includes the upper image endpoint and omits $\alpha/\pi^2$; they are distinct from `kernel.csv`.

## Commensurate action and the two additional forces

Equation (12) and the discussion immediately following it in [*Antiferromagnetic order enhanced by local dissipation*, Physical Review B 109, 205148 (2024)](https://doi.org/10.1103/PhysRevB.109.205148) give the zero-magnetization action $S=S_{\mathrm{LL}}+S_g+S_\alpha$:

$$
S_g=-\frac{g}{2\pi^2a\tau_c}\int dx\,d\tau\,\cos[4\phi(x,\tau)],
$$

$$
S_\alpha=-\frac{\alpha}{2\pi^2a\tau_c^{1-s}}
\int_{|\tau-\tau'|>\tau_c}dx\,d\tau\,d\tau'
\frac{\cos[2\phi(x,\tau)]\cos[2\phi(x,\tau')]}{|\tau-\tau'|^{1+s}}.
$$

Here $a$ and $\tau_c$ are the spatial and imaginary-time cutoffs. The identity

$$
\cos(2\phi)\cos(2\phi')=
\frac12\cos[2(\phi-\phi')]+\frac12\cos[2(\phi+\phi')]
$$

shows precisely what changes relative to the incommensurate action: the difference-field cosine remains, and **two additional terms** survive commensurability:

1. The local sine-Gordon (umklapp) term $S_g$, proportional to $-g\cos(4\phi)$.
2. A nonlocal **sum-field cosine**, proportional to $-\alpha\cos[2(\phi+\phi')]$, with the same long-range kernel as the difference-field cosine. Its prefactor in the action is $-\alpha/(4\pi^2a\tau_c^{1-s})$.

Taking the negative functional derivative adds two terms to the Langevin drift. In lattice units, using the stored kernel $\mathcal K_{\beta,s}$ defined above, `Diss_comm()` implements

$$
\begin{aligned}
F^{\mathrm{comm}}_{ij}
&=\frac{u}{\pi K}\Delta_x\phi_{ij}
+\frac{1}{\pi uK}\Delta_\tau\phi_{ij}\\
&\quad+\sum_{j'\ne j}\mathcal K_{\beta,s}(|j-j'|)
\sin[2(\phi_{ij'}-\phi_{ij})]\\
&\quad-\sum_{j'\ne j}\mathcal K_{\beta,s}(|j-j'|)
\sin[2(\phi_{ij'}+\phi_{ij})]
-g_{\mathrm{code}}\sin(4\phi_{ij}).
\end{aligned}
$$

The Langevin equation is $d\phi_{ij}/dt=F^{\mathrm{comm}}_{ij}+\eta_{ij}$, with the same noise normalization as above. Both additional forces are evaluated in each predictor–corrector stage. The command-line parameter `g` is the **paper's action coupling** in unit-cutoff lattice units; the runner converts it to $g_{\mathrm{code}}=2g/\pi^2$ before calling `Diss_comm()`. The default is $g=0$, leaving the local force inactive while the nonlocal sum-field force remains active.

## Verification and migration from the original programs

`make test` checks the public libraries directly for both models: the Ohmic kernel against its rational formula, generalized kernels against a beta-function recurrence, periodic Laplacians on Fourier modes, roughness on a rectangular lattice, interaction forces against finite-difference action derivatives, and a harmonic equilibrium variance benchmark. Executable-level tests compare split and uninterrupted runs, including checkpoint bytes and accumulated observables, and reject invalid inputs, changed restart parameters, wrong models, truncated checkpoints and output overwrites. Additional API tests verify both headers together, multi-translation-unit linkage, zero-field initialization, measurement scheduling, model validation, and exact checkpoint continuation. Tests of a small lattice do not establish equilibration, time-step convergence or finite-size scaling for production parameters.

GitHub Actions is configured to run the suite with GCC and Clang on pushes and pull requests. It installs Eigen only in the disposable CI runner. Production runs should additionally check time-step convergence, equilibration and autocorrelation for their chosen parameters.

The command-line and restart format replace the original hard-coded durations and filename-based restart scheme. Old matrices and auxiliary files are **not** accepted as complete checkpoints: they lack the random state and, in some cases, contain centered fields or invalid roughness data. Keep those data separately. The local force is still $-2g\sin(4\phi)/\pi^2$, but the CLI `g` now names the physical action coupling, rather than the old already-rescaled variable. Commensurate measurements no longer subtract the field mean. These corrections can change measured commensurate order parameters relative to the original programs.

## Related publications

1. Oscar Bouverot-Dupuis, Saptarshi Majumdar, Alberto Rosso, and Laura Foini, **Antiferromagnetic order enhanced by local dissipation**, *Physical Review B* **109**, 205148 (2024). [arXiv:2312.13095](https://arxiv.org/abs/2312.13095) · [Published article](https://doi.org/10.1103/PhysRevB.109.205148).
2. Saptarshi Majumdar, Laura Foini, Thierry Giamarchi, and Alberto Rosso, **Localization induced by spatially uncorrelated subohmic baths in one dimension**, *Physical Review B* **108**, 205138 (2023). [arXiv:2307.07989](https://arxiv.org/abs/2307.07989) · [Published article](https://doi.org/10.1103/PhysRevB.108.205138).
3. Saptarshi Majumdar, Laura Foini, Thierry Giamarchi, and Alberto Rosso, **Bath induced phase transition in a Luttinger liquid**, *Physical Review B* **107**, 165113 (2023). [arXiv:2210.01590](https://arxiv.org/abs/2210.01590) · [Published article](https://doi.org/10.1103/PhysRevB.107.165113).
