# Reproducing Section 6

## Use the recorded data

Extract the ZIP into a directory of your choice. It contains the prescribed orbit, the interpolation data, the two solver traces, the Hessian certificate, and the figures. [Download the Section 6 data](https://github.com/bqliu815/DFP-Counterexample/releases/download/section6-matched-data/DFP-Section6-data.zip)

To redraw the figures, start MATLAB in the repository root and point it to the directory that contains `geometry.mat`, `bfgs.mat`, and `dfp.mat`. For example, if you extracted the archive into `results/recorded`, run:

```matlab
dataDir = fullfile(pwd, 'results', 'recorded');
run_experiments('figures', dataDir);
```

This writes `Fig2.pdf` (trajectories) and `Fig3.pdf` (gradient norms), with PNG copies, in the selected directory. Figure 1 in the manuscript is a geometric schematic. To inspect the traces directly, open `bfgs.csv` and `dfp.csv`. Each row gives the iteration, objective value, original-coordinate gradient norm, and two coordinates.

## Repeat the experiment from scratch

The reference environment is MATLAB R2025a in double precision with a single computational thread. Running the solvers requires Optimization Toolbox and Statistics and Machine Learning Toolbox. The exact certificate uses Python 3 with NumPy and h5py. From the repository root:

```sh
python3 -m pip install numpy h5py
matlab -singleCompThread -batch "run_experiments('paper', fullfile(pwd,'results','section6'))"
```

The `paper` and `all` modes each generate the data, certificate, and both figures. The default output directory is `results/section6/`. Use a fresh output directory for an independent run. The same entry point also accepts `test` for the MATLAB tests and `smoke` for a short example.

### Figure 2. Prescribed DFP and MATLAB BFGS

The numerical example uses the family of Section 4. The Hessian Hölder exponent is $\beta=1/2$, the angular coefficients are $\bar\delta_{\mathrm E}=4\pi/9$ and $\bar\delta_{\mathrm O}=2\pi/9$, and $\nu_*=1.03\bar\delta_{\mathrm O}$. Starting at $j_0=32^3$, the code uses the integrated angular increments in (82) over 100,000 two-step cycles and reconstructs the preceding radii from their ratios with the terminal radius fixed to one. The resulting 200,000 prescribed steps pass the strong Wolfe checks with parameters $(1/4,3/4)$ and tolerance $10^{-12}$. The concrete analytic example in Proposition 4 uses $(1/4,1/2)$.

The objective interpolates all 200,001 sites before either solver runs. Its baseline is $\|x\|^2/2$. At site $x_k$, the gradient correction is $(\zeta_k-1)x_k$, supported in a disk whose radius is 0.24 times the nearest-neighbor distance. Figure 2(a) shows the prescribed orbit and Figure 2(b) shows the accepted BFGS iterates on that objective. The dashed circle in each panel has radius one.

### Figure 3. MATLAB gradient norms

Both solvers start from the prescribed $x_0$ and $H_0$. In the coordinates $x=x_0+Lz$ with $LL^T=H_0$, identity initialization at $z=0$ realizes this choice of $H_0$. Each call to `fminunc` uses the quasi-Newton algorithm, analytic gradients, its native line search and safeguards, and `HessUpdate` set to `dfp` or `bfgs`. The iteration limit is 5,000, the evaluation budget is 100,000, and the gradient stopping tolerance is $10^{-10}$ in the original coordinates. Figure 3 reuses the BFGS trace from Figure 2(b).

| Update | Accepted steps in the reference run | Final original-coordinate gradient norm |
| --- | ---: | ---: |
| BFGS | 23 | $5.39\times10^{-12}$ |
| DFP | 5,000 | $1.10818$ |

Solver histories may vary with MATLAB release and floating-point behavior. The recorded output files allow comparison at each accepted iterate.

### Hessian certificate

`matched_finite.mat` stores the two site coordinates, two gradient-correction components, and the support radius in each row of `fixedData`. The certificate treats these binary64 entries as exact dyadic rationals. It checks disjointness of the closed support disks using integer squared distances and combines the correction-to-radius bound with an analytic cutoff bound. The resulting global Hessian eigenvalue interval in the original coordinates is $[1937/12500,23063/12500]=[0.15496,1.84504]$.

The `hessian_band` in the solver summaries is the coarser estimate from objective construction. The paper reports the bound in `certificate.json`.

To recompute this certificate for the recorded objective, start MATLAB at the repository root, with the Python dependencies installed, and use the same extracted data directory as above:

```matlab
run_experiments('certify', dataDir);
```

To run the same check directly in a shell from the repository root:

```sh
python3 polar/certify_section6.py results/recorded/matched_finite.mat
```

The script writes `certificate.json` beside the MAT file. Use `--output PATH` for a different destination. To choose the Python interpreter used by the MATLAB entry point, set the `PYTHON` environment variable to its executable path.

## Earlier experiments

`run_experiments('legacy', outDir)` retains the earlier batch, including solver comparisons using 8,005 sites at starting indices $200^3$ and $400^3$. Use a separate output directory for those runs and their figure numbering. The optional `scales` and `wolfe` modes reproduce the earlier exponent and Wolfe-parameter diagnostics. The legacy symbolic certificate uses Symbolic Math Toolbox.
