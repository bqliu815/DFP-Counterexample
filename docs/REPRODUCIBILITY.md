# Reproducing the numerical experiments

## Environment and commands

The reference runs used MATLAB R2023a on an Intel Xeon Gold 6326 (2.90 GHz),
Rocky Linux 8.6, one computational thread, and 4 GiB allocated memory.
All computations use double precision and analytic gradients. DFP and BFGS
are implemented directly, with symmetrization after each update and no
damping, restart, or matrix reset.

The release was validated on 19 September 2026. MATLAB R2023a on Linux passed
all 14 tests, the smoke run, and the complete protocol. All 28 trajectory CSV
files were byte-identical to the original MATLAB results. The eight exact
certificates, all three CSV tables, and both PDF/PNG figures were also checked.
MATLAB R2025a on macOS also passed the 14 tests, smoke run, and complete
protocol. The [GitHub Actions workflow](https://github.com/bqliu815/DFP-Counterexample/actions/workflows/matlab.yml)
runs the 14 tests on MATLAB R2023a for pushes and pull requests.

The complete headless Linux run used:

```sh
env -u DISPLAY LC_ALL=C matlab -softwareopengl -nodisplay -singleCompThread \
  -batch "run_experiments('all')"
```

Some preliminary runs on this Linux installation stalled during MATLAB
startup or shutdown. The completed full run exited normally; the final test
run appended `quit(0,'force')` after the checks. The entry point initializes
the symbolic engine before the numerical stage in `all` mode.

Run commands from the repository root. The `all` mode performs `full`,
`certify`, `verify`, and `figures` in that order.

| Command | Purpose | Required input |
| --- | --- | --- |
| `run_experiments('probe')` | Record the MATLAB environment | — |
| `run_experiments('test')` | Run the test suite | — |
| `run_experiments('smoke')` | Run a short example | — |
| `run_experiments('full')` | Generate the paper data | — |
| `run_experiments('certify')` | Check exact finite-function bounds | `full` output |
| `run_experiments('verify')` | Check recorded steps and export tables | `full` and `certify` outputs |
| `run_experiments('figures')` | Export PDF and PNG figures | `full` output |

The Statistics and Machine Learning Toolbox is required by `test`, `smoke`,
and `full`; `certify` requires the Symbolic Math Toolbox. Figure export needs
the JVM, so omit `-nojvm`. To choose an output directory, use, for example,
`run_experiments('all', fullfile(tempdir, 'dfp-reproduction'))`. All stages
must use the same directory; repeated runs overwrite the corresponding output.
The caller's working directory and MATLAB search path are preserved.
The `smoke` mode runs 100 prescribed cycles and one finite-function BFGS
comparison. Dedicated tests for the updates, line searches, interpolation,
and entry point are in `tests/test_core.m` and run through the `test` mode.

## Section 6.1: geometry and asymptotic behavior

`DFPExperiment.oracle(100000, 0.03)` generates 100,000 two-step cycles.
The initial coordinates are

```text
r0 = epsilon0^2
p0 = 2 + (198/5)*epsilon0^3 - (9/5)*epsilon0^4
h0 = 1 + 8*epsilon0^3
H0 = diag([h0*p0*r0^2, h0])
g0 = x0 = [1; p0*r0]
```

The two legs use `(mu, tau) = (epsilon, 2/3)` and `(-2*epsilon, 1/3)`.
The eigenframe is recomputed before each leg. These steps follow the
prescribed recurrence; no line search is performed in this calculation.
The computed steps are checked against the Wolfe inequalities.

Figure 1(a) shows all 200,000 steps. Its dashed circle uses the last reference
center `C_k = x_k - g_k` and the radius estimate `G_N*exp(-13*epsilon_N/3)`.
Table 1 reports medians over the last 10,000 cycles and maximum residuals
over all steps. Figure 1(b) shows BFGS on a finite interpolant with 205 points
at `epsilon0 = 0.03`. It shares the initial point and matrix of Figure 1(a);
the two panels do not compare algorithms on the same finite function.
This larger-parameter interpolant is used for geometry and is outside the
uniform-convexity certificates below.

## Section 6.2: DFP and BFGS on the same objective function

`DFPExperiment.finite` forms a fixed finite objective function from the
prescribed endpoints. With `B = max(100, ceil(0.5*epsilon0^(-1.5)))`, the
interpolation prefix has `2*B + 4` steps and includes the initial point.
The quadratic term is centered at the final reference center; support radii
are 0.24 times the nearest-neighbor distances. The cutoff is one for
`t <= 1/3`, zero for `t >= 1`, and
`1 - (10*z^3 - 15*z^4 + 6*z^5)` in between, where `z = (t - 1/3)/(2/3)`.

Figure 2 uses `epsilon0 = 0.0025` and 8,005 interpolation points. Table 2
uses `epsilon0 = 0.001, 0.002, 0.0025`. Both methods start from the same
point and matrix on the same finite function at each parameter. DFP has a
5,000-iteration limit and BFGS a 1,000-iteration limit; the gradient tolerance
is `1e-10`. Invalid curvature, a non-descent direction, or line-search failure
terminates the run with a failure status.

All line searches use `c1 = 0.25` and `c2 = 0.75`:

| Policy | Search | First trial |
| --- | --- | --- |
| `weak_unit` | Weak Wolfe bracketing and bisection | 1 |
| `zoom_unit` | Strong Wolfe bracketing and cubic zoom | 1 |
| `mt_unit` | Strong Wolfe, More–Thuente | 1 |
| `zoom_history` | Strong Wolfe bracketing and cubic zoom | History-based |
| `mt_history` | Strong Wolfe, More–Thuente | History-based |

The main comparison uses `zoom_unit`. The five-policy DFP diagnostic uses
`epsilon0 = 0.001, 0.002`; weak/strong BFGS agreement is checked at all three
Table 2 parameters. History initialization starts at 1, then uses
`min(1, 1.01*2*(f_k - f_(k-1))/(g_k'*d_k))`, reverting to 1 for a negative
value. Strong Wolfe trial steps are capped at 64; the weak search uses 64
in its bracketing rule. The weak search allows 100 evaluations, cubic zoom
allows 40 outer and 11 zoom iterations, and More–Thuente allows 99 iterations.

The departure diagnostic runs DFP with `weak_unit` on the grid
`0.0005, 0.00065, 0.0008, 0.001, 0.00125, 0.0015, 0.002, 0.0025`.
At step `k`, `target_ratio` is the distance to the prescribed endpoint divided
by its support radius. The first exceedances of `1/3` and `1` give the plateau
and support departure indices. Runs stop at support departure or the budget
`B`; the plateau indices are fitted by log–log least squares.

These finite experiments illustrate the geometry and long DFP transients.
The infinite nonconvergence result is proved in the paper. The uniformly
convex finite interpolants have Lipschitz Hessians, so their long transients
are consistent with the planar convergence theorem.

## Certificates and recorded-step checks

`certify_bounds` treats the stored binary64 points, radii, and corrections
as exact dyadic rationals using `sym(value, 'f')`. It checks global Hessian
perturbation bounds and support separation for the eight small-parameter
finite functions. Both stored corrections and exact differences of stored
centers are checked. These certificates concern the fixed finite functions,
not floating-point trajectories or the infinite construction.

`verify_results` checks all 28 full-protocol trajectory files, the agreement
records, and the eight certificates. The checks include finite values,
curvature, positive definiteness, secant residuals, and Wolfe ratios. Recorded
Wolfe ratios use a scaled tolerance of `1e-11`; the line searches use their
direct computed inequalities when accepting steps. Successful verification
also exports the CSV tables.

## Output files

Paths below are relative to `results/paper/` or the chosen output directory.
The `test`, `smoke`, and `probe` modes default to their own directories.

| Path | Contents |
| --- | --- |
| `raw/geometry.mat`, `raw/geometry.json` | Prescribed recurrence and Table 1 quantities |
| `raw/finite_<parameter>.mat` | Finite objective coefficients and prescribed data |
| `raw/<run>.mat`, `.csv`, `.json` | Full run, accepted-step table, and summary |
| `raw/departure_fit.json`, `raw/agreements.json` | Departure fit and trajectory comparisons |
| `raw/certificate_<parameter>.json` | Exact finite-function bounds |
| `raw/verification.json`, `raw/environment_<mode>.json` | Recorded-step checks and software environment |
| `tables/Table1.csv`, `tables/Table2.csv`, `tables/Departure.csv` | Paper tables and departure indices |
| `figures/Fig1.pdf`, `figures/Fig2.pdf` | Vector figures; PNG copies are also exported |

Parameter tags replace the decimal point by `p`, as in `0p0025`.
Each `result.trace` row records an accepted step: iterate, gradient norm,
step length, support distances, Wolfe ratios, curvature, secant residual,
and matrix eigenvalue bounds. `result.trials{k}` retains all trials, including
rejections, as `[alpha, function_value, directional_derivative]`. The
`evaluations` column counts line-search trials; it excludes the repeated
evaluation of an accepted endpoint when recording the update.

To evaluate a saved finite function, rebuild its nearest-neighbor tree:

```matlab
addpath('src');
S = load('results/paper/raw/finite_0p0025.mat', 'obj', 'o');
obj = S.obj;
obj.tree = KDTreeSearcher(obj.points);
[f, g] = DFPExperiment.valueGrad(obj, S.o.x(1, :)');
```

`obj.band` is a floating-point estimate; exact bounds are in the certificate
file. [Reference summaries](../reference/expected_results.json) and
[Figure 1](../reference/figures/Fig1.pdf) / [Figure 2](../reference/figures/Fig2.pdf)
are supplied for comparison, not as computational inputs. Complete
trajectories are regenerated and excluded from version control. Finite-step
departure indices can vary with MATLAB release and floating-point behavior.
