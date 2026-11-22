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
include/           # Shared kernel, command-line runner, checkpoints and measurements
tests/             # Numerical tests and executable-level restart tests
scripts/           # Optional plotting utilities
figures/           # Example kernel figure
```

| Model | Initial executable | Continuation executable | Order parameter |
| --- | --- | --- | --- |
| Incommensurate | `build/initial-incommensurate` | `build/continue-incommensurate` | $\langle\cos(\phi-\overline\phi)\rangle$ |
| Commensurate | `build/initial-commensurate` | `build/continue-commensurate` | $\langle\cos(2\phi)\rangle$ |

Initial runs start from $\phi=0$. Continuation restores a complete checkpoint, including the unshifted field, absolute step, parameters, accumulated measurements, random-generator state, and the Gaussian distribution's cached state. The original source filenames are preserved under their model directories, including the spelling `continuted` in the incommensurate continuation source.

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

For even $\beta$, `Kernel()` stores

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

`make test` checks both models: the Ohmic kernel against its rational formula, generalized kernels against a beta-function recurrence, periodic Laplacians on Fourier modes, roughness on a rectangular lattice, interaction forces against finite-difference action derivatives, and a harmonic equilibrium variance benchmark. Executable-level tests compare split and uninterrupted runs, including checkpoint bytes and accumulated observables, and reject invalid inputs, changed restart parameters, wrong models, truncated checkpoints and output overwrites. Tests of a small lattice do not establish equilibration, time-step convergence or finite-size scaling for production parameters.

GitHub Actions is configured to run the suite with GCC and Clang on pushes and pull requests. It installs Eigen only in the disposable CI runner. Production runs should additionally check time-step convergence, equilibration and autocorrelation for their chosen parameters.

The command-line and restart format replace the original hard-coded durations and filename-based restart scheme. Old matrices and auxiliary files are **not** accepted as complete checkpoints: they lack the random state and, in some cases, contain centered fields or invalid roughness data. Keep those data separately. The local force is still $-2g\sin(4\phi)/\pi^2$, but the CLI `g` now names the physical action coupling, rather than the old already-rescaled variable. Commensurate measurements no longer subtract the field mean. These corrections can change measured commensurate order parameters relative to the original programs.

## Related publications

1. Oscar Bouverot-Dupuis, Saptarshi Majumdar, Alberto Rosso, and Laura Foini, **Antiferromagnetic order enhanced by local dissipation**, *Physical Review B* **109**, 205148 (2024). [arXiv:2312.13095](https://arxiv.org/abs/2312.13095) · [Published article](https://doi.org/10.1103/PhysRevB.109.205148).
2. Saptarshi Majumdar, Laura Foini, Thierry Giamarchi, and Alberto Rosso, **Localization induced by spatially uncorrelated subohmic baths in one dimension**, *Physical Review B* **108**, 205138 (2023). [arXiv:2307.07989](https://arxiv.org/abs/2307.07989) · [Published article](https://doi.org/10.1103/PhysRevB.108.205138).
3. Saptarshi Majumdar, Laura Foini, Thierry Giamarchi, and Alberto Rosso, **Bath induced phase transition in a Luttinger liquid**, *Physical Review B* **107**, 165113 (2023). [arXiv:2210.01590](https://arxiv.org/abs/2210.01590) · [Published article](https://doi.org/10.1103/PhysRevB.107.165113).
