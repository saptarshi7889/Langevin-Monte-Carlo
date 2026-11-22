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

## Current programs

Both programs implement the nonlinear dissipative sine force, periodic spatial and imaginary-time Laplacians, and a two-stage stochastic predictor–corrector update. The same matrix of independent standard Gaussian values $\xi$ is reused in both stages:

$$
\begin{aligned}
\phi^{\mathrm{pred}}&=\phi+\Delta t\,F(\phi)+\sqrt{2\Delta t}\,\xi,\\
\phi^{\mathrm{new}}&=\phi+\frac{\Delta t}{2}
\left[F(\phi)+F(\phi^{\mathrm{pred}})\right]+\sqrt{2\Delta t}\,\xi.
\end{aligned}
$$

Here $F$ is the deterministic drift. This replaces the Euler update described for the earlier source. In the C++ matrices, rows label imaginary time and columns label space, the reverse of the paper's index convention above. Parameters other than the run identifier and continuation input folder are edited directly in `main` and require recompilation.

| Program | Role | Current defaults |
| --- | --- | --- |
| [`Langevin_simulation_initial_simulation.cpp`](Langevin_simulation_initial_simulation.cpp) | Starts at $t=0$ with every field value zero. | $u=1$, $K=1$, $s=1.5$, $\alpha=1$, $L=\beta=16$, $\Delta t=0.05$, $t_{\mathrm{fin}}=16$ |
| [`Langevin_simulation_continuted_simulation.cpp`](Langevin_simulation_continuted_simulation.cpp) | Loads a previously saved field, order-parameter average, and sample count, then continues sampling. The filename retains the spelling `continuted`. | $u=1$, $K=1$, $s=0.5$, $\alpha=5$, $L=\beta=256$, $\Delta t=0.05$, $t_{\mathrm{init}}=0$, $t_{\mathrm{fin}}=102400$ |

### Kernel implemented in both programs

The Ohmic formula in the published-method section describes Appendix C of the Ohmic paper. The current `Kernel(b, R, a, c)` functions instead use the Euler beta function, with $R=\beta$, $a=\alpha$, and $c=s$. For the even lattice sizes used here, the array stores

$$
\mathcal{K}_{\beta,s}(r)=\frac{\alpha}{\pi^2}
\sum_{k=-\beta/2}^{\beta/2-1}
\mathrm{B}\left(r+k\beta-\frac{s}{2},\,s+1\right),
\qquad r=1,\ldots,\beta-1.
$$

Here $\mathrm{B}(x,y)$ denotes the beta special function called as `beta(x,y)` in C++; it is distinct from the imaginary-time extent $\beta$. The implementation evaluates the signed argument $r+k\beta-s/2$ exactly as written, without an absolute value. Negative first arguments occur in the image sum, so their numerical evaluation should be checked when changing compiler or standard-library implementation. No independent validation of this generalized kernel has been performed here.

`Diss()` multiplies these stored coefficients directly by $\sin[2(\phi_{ij'}-\phi_{ij})]$ and sums over other imaginary-time sites at the same spatial site. Thus the stored kernel **already includes** $\alpha/\pi^2$; do not apply that prefactor again. The equal-time term is excluded. The Ohmic beta-function identity

$$
\mathrm{B}\left(x-\tfrac12,2\right)=\frac{1}{x^2-\tfrac14}
$$

relates the $s=1$ expression to the earlier Ohmic kernel. Unlike the printed Appendix C sum, the code excludes the upper image endpoint $k=\beta/2$. The existing [kernel plot at $\beta=256$](figures/dissipative-kernel-beta256.png) and [plotting script](scripts/plot_dissipative_kernel.py) illustrate the printed Ohmic Appendix C formula, without $\alpha/\pi^2$; they do not plot the generalized beta-function kernel above.

### Initial runs and saved observables

The initial program accepts `id=N` (default `0`). Output filenames use `N+1`, including `coreN+1` in configuration filenames. At each measurement, it subtracts the field's spatial–imaginary-time mean to obtain $\phi^\star$, saves that centered configuration, and accumulates spatial and imaginary-time roughness and the order parameter $\langle\cos\phi^\star\rangle$. It also writes a time log, running averages, a sample count, and a kernel file. Configuration files are plain whitespace-separated matrices with `Beta` rows and `L` columns.

**Before generating continuation input, increase the initial run duration:** measurements start at `t_eq1 = 4*t_eq2 = 64`, while the current `t_fin` is only `16`. These defaults therefore produce no saved configurations or measured observables. Measurements are scheduled every `int(t_eq2/delta_t)` steps, with `t_eq2 = 16`. Keep $L=\beta$ with the present roughness-array initialization and output loops.

### Continuing a run

The continuation program accepts `id=N` and `folder=/path/to/input/`. Include the trailing `/`: the code concatenates the folder and filenames directly. It reads the order-parameter file (second value on its data line), sample-count file, and configuration whose index equals that saved count. It reconstructs the accumulated order parameter from its previous average and count, then adds subsequent measurements. It does not restore the random-generator state, so this is a statistical continuation rather than a bit-for-bit restart.

Match the input's $L$, `Beta`, $K$, $s$, $\alpha$, time step, and run identifier. The supplied initial and continuation defaults **do not match**. Input filenames are constructed using `t_fin_prev`, currently assigned from the continuation's `t_fin`; set `t_fin_prev` to the duration tag in the existing filenames when extending to a different final time. The input configuration index comes from the saved count, not from `t_init` or a time log. Set `t_init` explicitly to the intended starting Langevin time; its current value is `0`. Changing it alone does not select a different input configuration.

Continuation measures every `int(t_eq/delta_t)` steps, with `t_eq = 256` and the condition `t > 0`. It evaluates the order parameter on the centered field, but saves the uncentered `Phi`. **Continuation roughness output is currently invalid:** loading, restoring, and updating `B_X` and `B_T` are commented out, but those uninitialized arrays are still written to a roughness file. Do not use that output for analysis. Existing order-parameter and sample-count files are replaced with updated totals; the time log is opened afresh.

### Running independent simulations in parallel

Neither source contains MPI, OpenMP, shared-memory work sharing, or a reduction across processes. Parallelism is at the level of **independent trajectories**: launch separate serial processes with distinct `id` values, and combine their statistics afterward. Each process creates its own random generator from `std::random_device`; `id` controls filenames, not the random seed. There is no seed argument or built-in ensemble aggregation.

Use a separate working directory for every trajectory, even with distinct IDs: the initial program's kernel filename has no ID, so concurrent runs in one directory would overwrite it. Outputs are always written in the current working directory. For continuation, use a new output directory and retain the original files in a separate input directory; the `folder` argument changes only the input location. Reuse the same ID for each trajectory's continuation. Do not launch the binaries with `mpirun -n N` and identical arguments: no MPI rank is read, so the processes would share filenames and overwrite outputs.

For example, after compiling the programs as `build/initial` and `build/continue`, launch two initial trajectories from the repository root:

```bash
mkdir -p runs/trajectory-0 runs/trajectory-1
(cd runs/trajectory-0 && ../../build/initial id=0) &
(cd runs/trajectory-1 && ../../build/initial id=1) &
wait
```

After adjusting and recompiling the continuation parameters to match a saved run, continue trajectory 0 in a new directory (replace the input path below with the absolute path to its saved files):

```bash
mkdir -p runs/trajectory-0-continuation
(cd runs/trajectory-0-continuation && ../../build/continue id=0 folder=/absolute/path/to/runs/trajectory-0/)
```

## Dependencies and validation

The programs require Eigen headers and a C++ compiler providing the C++17 beta special function. They also use variable-length stack arrays, a compiler extension rather than standard C++. MPI is not a dependency of these two sources. With an existing Eigen installation, example GNU C++ build commands are:

```bash
mkdir -p build
g++ -O2 -std=gnu++17 -I/path/to/eigen3 Langevin_simulation_initial_simulation.cpp -o build/initial
g++ -O2 -std=gnu++17 -I/path/to/eigen3 Langevin_simulation_continuted_simulation.cpp -o build/continue
```

Replace `/path/to/eigen3` with the directory containing `Eigen/`. No dependencies have been installed, and these two programs have not been compiled or numerically validated here. The descriptions above are based on source inspection; the kernel plot illustrates a formula, not a simulation benchmark. No separate commensurate pinning term for the zero-magnetization problem is implemented.

## Related publications

1. Oscar Bouverot-Dupuis, Saptarshi Majumdar, Alberto Rosso, and Laura Foini, **Antiferromagnetic order enhanced by local dissipation**, *Physical Review B* **109**, 205148 (2024). [arXiv:2312.13095](https://arxiv.org/abs/2312.13095) · [Published article](https://doi.org/10.1103/PhysRevB.109.205148).
2. Saptarshi Majumdar, Laura Foini, Thierry Giamarchi, and Alberto Rosso, **Localization induced by spatially uncorrelated subohmic baths in one dimension**, *Physical Review B* **108**, 205138 (2023). [arXiv:2307.07989](https://arxiv.org/abs/2307.07989) · [Published article](https://doi.org/10.1103/PhysRevB.108.205138).
3. Saptarshi Majumdar, Laura Foini, Thierry Giamarchi, and Alberto Rosso, **Bath induced phase transition in a Luttinger liquid**, *Physical Review B* **107**, 165113 (2023). [arXiv:2210.01590](https://arxiv.org/abs/2210.01590) · [Published article](https://doi.org/10.1103/PhysRevB.107.165113).
