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
| `all` | Generate all diagnostics and solver comparisons, certify finite functions, and export figures | None |
| `scales` | Compare `beta=1/3,1/2,2/3`, with 10,000,000 cycles each | None |
| `wolfe` | Check 10,000 cycles for each of `(0.25,0.75)`, `(0.2,0.3)`, and `(0.8,0.9)` | None |
| `certify` | Repeat the exact finite-function checks | Stored finite MAT files from `all` |
| `figures` | Export vector PDF and PNG figures | Stored trajectory MAT files from `all` |

The default mode is `smoke`. Its outputs go to `results/smoke/`; the complete
run and its follow-up stages use `results/paper/`. A second argument selects
an output directory. Results are generated locally and excluded from version
control. Figure export requires MATLAB with the JVM.

## Prescribed polar recurrence

All current runs use `phi=p*gamma_star/(gamma*nu(gamma)*(1/a+1/b))`, where
`p=1/(1+beta)` and `nu(gamma)` is the local increasing inverse of
`gamma(nu)=(b-a+2*nu)/(nu*(a-nu)*(b+nu))`. The same cycle-start feedback is
used in both stages. The code follows this inverse by Newton iteration and
checks the branch and residual.
The default shape keeps `a=4*pi/9`, `b=2*pi/9`, `nu_star=1.03*(2*pi/9)`,
and `beta=1/2` for the geometry and finite-solver experiments.

`PolarDFP.orbit(100000, 32^3, 1.03)` generates `Fig1.pdf` panel (a) and
the recurrence-verification table.
Each cycle contains two prescribed steps. No line search is performed.
All radii are multiplied by one common factor so that the final recorded
radius is one. This changes neither the matrix recurrence nor the secant
and line-search ratios. In the paper's notation, the recorded radii are
`rho_k^(N) = rho_k / rho_(2N)`. The dashed unit circle marks this normalized terminal
radius. This finite normalization does not compute the infinite sequence's
accumulation radius.

Figure 1(a) retains all recorded points. The inset magnifies the boxed region
using the same coordinates and equal axis scaling; the arrow identifies its
location. No points are removed and no smoothing is applied.

Table 1 uses medians over the last 10,000 cycles for the first four rows.
The secant and update residuals are maxima over all steps. Armijo and strong
Wolfe ratios are computed from the quadratic endpoint values and analytic
gradients. The failure counts use `(c1,c2) = (0.25,0.75)` and tolerance `1e-12`.
The first four rows compare finite-run records with asymptotic limits.
The zero reference values for the secant and update residuals instead come
from exact identities. The remaining rows diagnose the Wolfe endpoint
conditions.

The `scales` stage fixes the representative angle coefficients and Wolfe pair
`(0.25,0.75)` and runs `beta=1/3,1/2,2/3` with `p=1/(1+beta)`,
`J=32^3`, and 10,000,000 cycles each. The feedback uses the corresponding `p`.
`Fig3.pdf` plots the stored arrays `ratio_radial` and `ratio_correction`:

- `log(r_(2j)*r_(2j+1))/((r_E+r_O)*tau_j^2)`;
- `(J+j)*(1-zeta_(2j+1))/phi(gamma_star)`.

`scaling_diagnostics` executes every recurrence cycle and stores logarithmically
spaced samples, the first 101 cycles, and the final cycle. Matrix validity,
Wolfe inequalities, and the two ratio extrema are checked over all steps.
The more expensive secant and DFP update residuals are checked only at the
stored cycles and are labelled `sampled` in the JSON summaries. Small blocks
of scale values avoid allocating full trajectory arrays for these long runs.
The test suite compares the sampled diagnostics and extrema with `PolarDFP.orbit`.

The implementation uses `log1p` and `expm1` to avoid cancellation. Both ratios
have theoretical limit one. In these runs they initially move away from one,
reach a maximum, and then decrease toward it. The initial state lies on the
limiting matrix family, but finite-angle terms produce a transient displacement
before the restoring feedback reduces it. Larger Hölder exponents give slower
decay, consistent with the bound `O((J+j)^(-(1-beta)/(1+beta)))`.
The figure connects all saved samples without smoothing and uses a logarithmic
cycle axis; cycle zero is omitted from the plot. These finite diagnostics
illustrate the recurrence scales, not global Hölder regularity or an infinite limit.

The `wolfe` stage uses `J=400^3`, `beta=1/2`, and 10,000 cycles for each
of the following settings:

| Wolfe pair | `(a,b,nu_star)` |
| --- | --- |
| `(0.25,0.75)` | `(4*pi/9,2*pi/9,1.03*(2*pi/9))` |
| `(0.2,0.3)` | `(7,13,3)` |
| `(0.8,0.9)` | `(7/12,1/6,1/4)` |

The latter two shapes are selected by `PolarDFP.wolfeParameters`. The table
in Section 6.1 reports the minimum Armijo and maximum strong Wolfe curvature
ratios over all recorded steps. No finite-function solver comparison is
performed for these diagnostic shapes. The complete `all` mode includes
both `scales` and `wolfe`.

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
The certified Hessian bounds refer to the original coordinates. At the
initial point, the objective Hessian in the solver coordinates has condition
number approximately `1.716e9` for `J = 200^3` and `2.745e10` for `J = 400^3`.

The common budget is 5,000 iterations and 100,000 solver evaluations.
An output function stops when the original gradient norm reaches `1e-10`.
`TolFun` and `TolX` are zero. Native termination messages are also recorded;
a positive solver exit flag alone is not treated as reaching the gradient target.

Figure 1(b) and Figure 2 use the finite interpolant with `J = 400^3`.
The Figure 1 titles are “(a) Prescribed DFP recurrence” and
“(b) BFGS on a finite interpolant”. The two panels use different values of the
paper's starting index `j_0` (the code parameter `J`): `32^3` in panel (a) and
`400^3` in panel (b). They display the respective experiments in Sections 6.1
and 6.2. Panel (a) shows the prescribed geometric trajectory; panel (b) shows
the BFGS trajectory generated by MATLAB on the fixed finite interpolant.
Figure 2 compares the solvers for `J = 400^3`; the finite-comparison table
reports both methods for each of the two finite interpolants.

To redraw the figures without running either numerical experiment, place the
stored `geometry.mat`, `dfp.mat`, and `bfgs.mat` files in one output directory
and run:

```matlab
run_experiments('figures', '/absolute/path/to/stored/data');
```

This mode redraws the stored trajectories without rerunning the recurrence,
solvers, or certificates. If the three `exponent_case_*.mat` files are also
present, it redraws `Fig3.pdf` from them.

`certify_finite` interprets the stored binary64 coordinates, corrections,
and radii as exact dyadic rationals. It checks support separation through
disjoint projections on the second coordinate and certifies global Hessian
bounds using rational norm comparisons and a rational upper bound for `sqrt(3)`.
These certificates apply to the fixed finite interpolants.

## Output files

| Output | Contents |
| --- | --- |
| `geometry.mat`, `geometry.json` | Prescribed recurrence and Table 1 diagnostics |
| `exponent_case_1.mat/.json` through `exponent_case_3.mat/.json` | Three-exponent samples, exact ratio extrema, and run summaries |
| `scales_environment.json` | MATLAB release, architecture, thread count, and precision for the exponent runs |
| `wolfe_reference.mat/.json`, `wolfe_case_1.mat/.json`, `wolfe_case_2.mat/.json` | Three Wolfe-parameter diagnostics |
| `finite.mat`, `finite_200.mat` | Fixed finite objective functions and prescribed data |
| `finite.json`, `finite_200.json` | Finite-orbit summaries and floating-point Hessian estimates |
| `finite_fixed_data.csv`, `finite_200_fixed_data.csv` | Stored coordinates, corrections, and support radii |
| `dfp.mat`, `bfgs.mat` and `_200` counterparts | Solver traces and summaries |
| `dfp.csv`, `bfgs.csv` and `_200` counterparts | Iterations, values, gradient norms, and coordinates |
| `dfp.json`, `bfgs.json` and `_200` counterparts | Termination, iterations, final gradients, and timings |
| `certificate.json`, `finite_200_certificate.json` | Exact finite-function certificates |
| `environment.json` | MATLAB release, toolboxes, and execution settings |
| `Fig1.pdf`, `Fig2.pdf`, `Fig3.pdf` | Trajectories, solver gradient norms, and exponent diagnostics; PNG copies included |
| `figures_complete.json` | Figure-export completion and plotting metadata |
| `complete.json` | Successful completion of the full run |

The manuscript includes an additional geometric schematic. Its Figures 2, 3,
and 4 correspond to the code outputs `Fig1.pdf`, `Fig3.pdf`, and `Fig2.pdf`,
respectively. The LaTeX asset filenames are `Fig2.pdf`, `Fig4.pdf`, and
`Fig3.pdf`; labels determine their printed numbers.

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
