# Polar DFP experiments

MATLAB implementation of the prescribed cubic-root polar construction and the
finite-function comparisons in Sections 6.1–6.2.

## Requirements

MATLAB with Optimization Toolbox and Statistics and Machine Learning Toolbox.
The exact Hessian certificates also require Symbolic Math Toolbox. The reported
experiments use MATLAB R2025a Update 1 with one computational thread.

## Reproduction

Add this directory to the MATLAB path and run:

```matlab
run_polar_experiments('smoke');
run_polar_experiments('all');
```

The complete run generates the prescribed sequence, compares the DFP and BFGS
options of `fminunc`, checks the finite objective functions with exact rational
arithmetic, and exports both paper figures. An optional second argument selects
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
- `make_polar_figures.m`: figure export.

The prescribed recurrence uses fixed step lengths and does not perform a line
search. A common scale sets the final recorded radius to one. The DFP and BFGS
comparisons use fixed finite interpolants and MATLAB's native line search and
update safeguards. Both solvers retain the prescribed initial search direction
through the same linear change of variables. Figure 1 pairs the prescribed
DFP sequence with the BFGS iterates on the finite interpolant with J = 400^3;
the two panels come from the respective experiments in Sections 6.1 and 6.2.

Native non-descent-direction errors are recorded with the last accepted
iterate; other exceptions are rethrown. Solver termination is reported
separately from the gradient tolerance.
