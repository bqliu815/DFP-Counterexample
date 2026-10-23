# Polar DFP experiments

MATLAB implementation of the prescribed polar construction with the
feedback used in the manuscript and the finite-function comparisons in Sections
6.1–6.2. The geometry and finite-solver comparisons retain the `2:1` angle
shape and `beta=1/2`. Separate diagnostics vary the Hölder exponent or the
Wolfe parameter pair.

## Requirements

MATLAB with Optimization Toolbox and Statistics and Machine Learning Toolbox.
The exact Hessian certificates also require Symbolic Math Toolbox. The reported
experiments use MATLAB R2025a Update 1 with one computational thread.

## Reproduction

Add this directory to the MATLAB path and run:

```matlab
run_polar_experiments('smoke');
run_polar_experiments('all');
run_polar_experiments('scales'); % beta = 1/3, 1/2, 2/3
run_polar_experiments('wolfe'); % three Wolfe parameter pairs
```

The complete run generates the prescribed sequence, compares the DFP and BFGS
options of `fminunc`, checks the finite objective functions with exact rational
arithmetic, and exports three numerical figures. An optional second argument selects
the output directory. Generated data are written to `results/` by default.

For a single-threaded run from this directory:

```sh
matlab -singleCompThread -batch "run_polar_experiments('all')"
```

To repeat the certificates or figures from generated data:

```matlab
run_polar_experiments('certify');
run_polar_experiments('figures');
```

## Implementation

- `PolarDFP.m`: normalized matrix recurrence, finite interpolation, and solvers.
- `run_polar_experiments.m`: experiment settings and execution.
- `certify_finite.m`: exact separation and Hessian-bound checks for stored data.
- `make_polar_figures.m`: trajectory and solver figure export.
- `scaling_diagnostics.m`: long recurrence runs with bounded output storage.
- `make_scaling_figure.m`: two-panel exponent comparison.

The prescribed recurrence is evaluated directly, without a line search.
A common scale sets the final recorded radius to one. The dashed unit circle
marks this normalized terminal radius. The DFP and BFGS comparisons use fixed
finite interpolants and MATLAB's native line search and
update safeguards. Both solvers retain the prescribed initial search direction
through the same linear change of variables. Figure 1 pairs the prescribed
DFP sequence at J = 32^3 with the BFGS iterates on the finite interpolant at
J = 400^3; the two panels come from the respective experiments in Sections 6.1
and 6.2. Their titles are “Prescribed DFP recurrence” and “BFGS on a finite
interpolant”. Figure 2 compares the solvers for `J = 400^3`. Figure 3 plots
the radial and correction ratios for three exponents at `J = 32^3` with
10,000,000 cycles each. Every cycle is evaluated, and both ratios initially
move away from one before turning back toward it. The saved logarithmic
samples are plotted without smoothing; the cycle axis is logarithmic.
These output filenames differ from the manuscript figure
numbers; the manuscript also includes a geometric schematic.
The `figures` mode only redraws stored data and
does not rerun the recurrence, solvers, or certificates.

Native non-descent-direction errors are recorded with the last accepted
iterate; other exceptions are rethrown. Solver termination is reported
separately from the gradient tolerance.
