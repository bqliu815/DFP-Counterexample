# Reproducing the numerical experiments

## Environment and commands

The reference runs use MATLAB R2023a on an Intel Xeon Gold 6326 (2.90 GHz),
Rocky Linux 8.6, one computational thread, and 4 GiB allocated memory.
The stored R2023a trajectories were checked and plotted with MATLAB R2025a
on macOS; the per-stage environment records identify these releases.
Optimization Toolbox and Statistics and Machine Learning Toolbox are required
for the experiments and tests. Symbolic Math Toolbox is required for the
exact finite-function certificates.

```matlab
run_experiments('test');
run_experiments('all');
```

`all` runs `full`, `certify`, `verify`, and `figures` in order. A second
argument selects the output directory. Otherwise the complete protocol writes
to `results/paper/`; `test`, `smoke`, and `probe` use separate subdirectories.
Generated results are excluded from version control.

| Mode | Purpose | Required input |
| --- | --- | --- |
| `probe` | Record the MATLAB environment | — |
| `test` | Run the test suite | — |
| `smoke` | Run 100 prescribed cycles and a finite-function BFGS comparison | — |
| `full` | Generate the paper trajectories and finite functions | — |
| `certify` | Check exact Hessian bounds and support separation | `full` |
| `verify` | Reevaluate recorded iterates and export tables | `full`, `certify` |
| `figures` | Export vector PDF and PNG figures | `full` |

For a single-threaded headless run, use:

```sh
matlab -nodisplay -singleCompThread -batch "run_experiments('all')"
```

Figure export requires the JVM. The numerical stage can also be run separately with `-nojvm`. The caller's working directory and
MATLAB search path are preserved. The GitHub Actions workflow runs the tests
on MATLAB R2023a with Optimization Toolbox and Statistics and Machine
Learning Toolbox.

## Prescribed recurrence

`DFPExperiment.oracle(100000, 0.03)` generates the 100,000 two-step cycles
in Figure 1(a). Step lengths are prescribed; this calculation does not call
an optimization solver or perform a line search. Its initial data are

```text
r0 = epsilon0^2
p0 = 2 + (198/5)*epsilon0^3 - (9/5)*epsilon0^4
h0 = 1 + 8*epsilon0^3
H0 = diag([h0*p0*r0^2, h0])
g0 = x0 = [1; p0*r0]
```

The two legs use `(mu, tau) = (epsilon, 2/3)` and `(-2*epsilon, 1/3)`.
The eigenframe is recomputed before each leg. The explicit DFP update in
`DFPExperiment.m` is used only for this prescribed recurrence.

The dashed circle has center `C_(2N)`, with `N = 100000`, and radius
`G_N*exp(-13*epsilon_N/3)`. Table 1 reports medians over the last 10,000 cycles
and maximum algebraic residuals over all steps. The Armijo diagnostic uses
surrogate endpoint values `0.5*norm(x_k - C_(2N))^2`; its counts and the
strong-curvature counts use `(c1,c2) = (0.25,0.75)` with tolerance `1e-12`.
The infinite construction's Wolfe conditions are proved in Section 4.3.

## Finite objective functions

`DFPExperiment.finite(epsilon0)` forms a fixed finite objective function.
With `B = max(100, ceil(0.5*epsilon0^(-1.5)))`, the interpolation prefix has
`2*B + 4` steps and includes the initial point. The quadratic term is centered
at the last reference center. Support radii are 0.24 times the nearest-neighbor
distances. The cutoff is one for `t <= 1/3`, zero for `t >= 1`, and
`1 - 10*z^3 + 15*z^4 - 6*z^5` in between, where `z = (3*t - 1)/2`.
The function and its analytic gradient are evaluated by `valueGrad`.

Figure 1(b) uses the 205-point interpolant at `epsilon0 = 0.03`, for which
no global-convexity certificate is claimed. Figure 2 uses `epsilon0 = 0.0025`
and 8,005 interpolation points. Table 2 uses `epsilon0 = 0.001, 0.002, 0.0025`.
The latter three finite functions have certified global Hessian bounds.

## fminunc comparisons

All finite-function optimization runs call `fminunc` with
`Algorithm='quasi-newton'`, an analytic gradient, and `HessUpdate='dfp'` or
`'bfgs'`. Both methods have the same budget of 5,000 iterations and 100,000
solver function evaluations. An output function stops when the gradient
norm in the original coordinates is at most `1e-10`, or after 5,000 iterations.
`TolFun` and `TolX` are zero so that the common stopping test is not replaced
by a tolerance in transformed coordinates. Native line-search termination
can still stop a run earlier; its exit flag and message are saved separately.
A positive exit flag alone is not counted as reaching the gradient target.

Two initializations are recorded for each Table 2 parameter:

- `prescribed`: minimize `f(x0 + L*z)` from `z = 0`, where `L*L' = H0`.
  The gradient passed to MATLAB is `L'*grad_f`. Its initial direction maps
  to `-H0*grad_f` in the original coordinates. This is the main comparison.
- `identity`: use `L = I`, giving identity initialization in the original
  coordinates. These controls are exported in `tables/Identity.csv`.

The same transform is used for DFP and BFGS. MATLAB's built-in line search,
first-update scalar rescaling, and curvature safeguards are left unchanged.
In R2023a, the internal line-search parameters are `rho = 0.01` and
`sigma = 0.9`. The adapter does not supply a custom line search, prescribe
accepted steps, or modify MATLAB's implementation. Thus these runs measure
the built-in solver, while Figure 1(a) measures the prescribed recurrence.

All plotted gradient norms and stated Hessian bounds refer to the original
coordinates. The finite functions have Lipschitz Hessians; these experiments
are not numerical proofs of nonconvergence of the infinite construction.

## Certificates, checks, and outputs

`certify_bounds` interprets the stored binary64 coefficients as exact dyadic
rationals. It verifies support separation and global Hessian bounds for the
three small-parameter functions, checking both stored corrections and exact
differences of stored centers. The certificates concern the finite functions,
not floating-point trajectories or the infinite construction.

`verify_results` independently reevaluates all recorded iterates from the
13 optimization runs, checks final gradients and stopping classifications,
and checks the recurrence diagnostics and three exact certificates. It does
not require one method to outperform another or turn a native stopping flag
into a convergence claim.

| Output | Contents |
| --- | --- |
| `raw/geometry.mat`, `.json` | Prescribed recurrence and Table 1 diagnostics |
| `raw/finite_<parameter>.mat` | Fixed finite function and prescribed data |
| `raw/<run>.mat`, `.csv`, `.json` | Recorded iterates, native output, options, and summary |
| `raw/certificate_<parameter>.json` | Exact finite-function bounds |
| `raw/verification.json` | Independent result checks |
| `raw/environment_<mode>.json` | MATLAB release, toolboxes, and execution settings |
| `tables/Table1.csv`, `Table2.csv`, `Identity.csv` | Recurrence and solver comparisons |
| `figures/Fig1.pdf`, `Fig2.pdf` | Vector figures, with PNG copies |

Each trajectory row records the iteration, point, function value, original
gradient, line-search step, cumulative solver evaluations, and step norm.
Armijo and curvature ratios are recorded as diagnostics, with `NaN` when the
computed step has no positive descent denominator. `function_evaluations`
counts MATLAB's solver calls. Evaluations used to monitor the original
gradient are recorded separately; setup and final checks are not solver calls.
Native `output.iterations` can include a failed final attempt, so the number
of recorded steps is also saved.

To evaluate a stored objective, rebuild its nearest-neighbor tree:

```matlab
addpath('src');
S = load('results/paper/raw/finite_0p0025.mat', 'obj', 'o');
obj = S.obj;
obj.tree = KDTreeSearcher(obj.points);
[f, g] = DFPExperiment.valueGrad(obj, S.o.x(1, :)');
```

`obj.band` is a floating-point estimate; the exact bounds are given in the
certificate. Iteration counts and stopping reasons can vary with MATLAB
release and floating-point behavior.
