# Reproducing the numerical experiments

## Environment and commands

The reference experiments use MATLAB R2025a Update 1 on an Apple M5 processor
with one computational thread and IEEE 754 double precision. Optimization
Toolbox and Statistics and Machine Learning Toolbox are required. The exact
certificates also require Symbolic Math Toolbox.

From the repository root:

```matlab
run_experiments('test');
run_experiments('smoke');
run_experiments('all');
```

For the single-threaded setup used in the paper:

```sh
matlab -singleCompThread -batch "run_experiments('all')"
```

| Mode | Purpose | Required input |
| --- | --- | --- |
| `test` | Run the MATLAB test suite | None |
| `smoke` | Run 200 prescribed cycles and a small BFGS comparison | None |
| `all` | Generate both experiments, certify finite functions, and export figures | None |
| `certify` | Repeat the exact finite-function checks | Stored finite MAT files from `all` |
| `figures` | Export vector PDF and PNG figures | Stored trajectory MAT files from `all` |

The default mode is `smoke`. Its outputs go to `results/smoke/`; the complete
run and its follow-up stages use `results/paper/`. A second argument selects
an output directory. Results are generated locally and excluded from version
control. Figure export requires MATLAB with the JVM.

## Prescribed polar recurrence

`PolarDFP.orbit(100000, 32^3, 1.03)` generates Figure 1(a) and Table 1.
Each cycle contains two prescribed steps. No line search is performed.
All radii are multiplied by one common factor so that the final recorded
radius is one. This changes neither the matrix recurrence nor the secant
and line-search ratios. The dashed curve is the unit circle.

Figure 1(a) retains all recorded points. The inset magnifies the boxed region
using the same coordinates and equal axis scaling; the arrow identifies its
location. No points are removed and no smoothing is applied.

Table 1 uses medians over the last 10,000 cycles for the first four rows.
The secant and update residuals are maxima over all steps. Armijo and strong
Wolfe ratios are computed from the quadratic endpoint values and analytic
gradients. The failure counts use `(c1,c2) = (0.25,0.75)` and tolerance `1e-12`.
These finite-run diagnostics are distinct from the asymptotic limits.

## Fixed finite objective functions

The two solver comparisons use `J = 200^3` and `J = 400^3`, each with 4,002
cycles and 8,005 interpolation points. The terminal radius is normalized to
one. The quadratic term is `0.5*norm(x)^2`, the corrections are
`(zeta_k - 1)*x_k`, and each support radius is 0.24 times the nearest-neighbor
distance. The radial cutoff is one on `[0,1/3]`, zero on `[1,infinity)`, and
`1 - 10*z^3 + 15*z^4 - 6*z^5` between them, where `z = (3*t - 1)/2`.
The objective function is fixed before either solver starts.

Both methods use `fminunc` with `Algorithm='quasi-newton'`, an analytic
gradient, and `HessUpdate='dfp'` or `'bfgs'`. The change of variables
`x = x0 + L*z`, where `L*L' = H0`, preserves the prescribed initial search
direction. MATLAB retains its native line search and update safeguards.
The common budget is 5,000 iterations and 100,000 solver evaluations.
An output function stops when the original gradient norm reaches `1e-10`.
`TolFun` and `TolX` are zero. Native termination messages are also recorded;
a positive solver exit flag alone is not treated as reaching the gradient target.

Figure 1(b) and Figure 2 use the finite interpolant with `J = 400^3`.
Thus the two panels of Figure 1 display the respective experiments in
Sections 6.1 and 6.2. Figure 2 compares the solvers on the same fixed function.

`certify_finite` interprets the stored binary64 coordinates, corrections,
and radii as exact dyadic rationals. It checks support separation through
disjoint projections on the second coordinate and certifies global Hessian
bounds using rational norm comparisons and a rational upper bound for `sqrt(3)`.
These certificates apply to the fixed finite interpolants.

## Output files

| Output | Contents |
| --- | --- |
| `geometry.mat`, `geometry.json` | Prescribed recurrence and Table 1 diagnostics |
| `finite.mat`, `finite_200.mat` | Fixed finite objective functions and prescribed data |
| `dfp.mat`, `bfgs.mat` and `_200` counterparts | Solver traces and summaries |
| `dfp.csv`, `bfgs.csv` and `_200` counterparts | Iterations, values, gradient norms, and coordinates |
| `dfp.json`, `bfgs.json` and `_200` counterparts | Termination, iterations, final gradients, and timings |
| `certificate.json`, `finite_200_certificate.json` | Exact finite-function certificates |
| `environment.json` | MATLAB release, toolboxes, and execution settings |
| `Fig1.pdf`, `Fig2.pdf` | Vector figures, with PNG copies |
| `complete.json` | Successful completion of the full run |

To evaluate a stored objective function:

```matlab
addpath('polar');
S = load('results/paper/finite.mat', 'objective', 'finiteOrbit');
S.objective.tree = KDTreeSearcher(S.objective.points);
[f, g] = PolarDFP.valueGrad(S.objective, S.finiteOrbit.x(1, :)');
```

The floating-point estimate `objective.band` is separate from the exact
bounds in the certificate. Iteration counts and termination can vary with
MATLAB release and floating-point behavior. The GitHub Actions workflow
runs the current MATLAB tests on R2025a.
