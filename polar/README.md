# Polar DFP implementation

The revised manuscript uses one prescribed orbit to build the finite objective for both numerical figures. `PolarDFP.orbit(100000, 32^3, 1.03)` returns the 100,000 two-step cycles. `PolarDFP.finite` builds the interpolant at all 200,001 sites, and `PolarDFP.solve` runs MATLAB’s DFP or BFGS update from the prescribed starting data.

From the repository root:

```matlab
outDir = fullfile(pwd, 'results', 'section6');
run_experiments('paper', outDir);
```

`figures` redraws the trajectories and gradient norms from the stored run. `certify` verifies support separation and the global Hessian bounds on the stored interpolation data. See the [reproduction guide](../docs/REPRODUCIBILITY.md) for inputs and outputs and the [recorded Section 6 dataset](https://github.com/bqliu815/DFP-Counterexample/releases/download/section6-matched-data/DFP-Section6-data.zip) for the reference figures and traces.

The `legacy`, `scales`, and `wolfe` modes preserve experiments from earlier revisions. The legacy MATLAB routine `certify_finite.m` uses a separation test tailored to those smaller interpolation sets.
