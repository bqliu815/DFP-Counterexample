# DFP-Counterexample

This repository accompanies the paper

*A counterexample to global convergence of classical DFP under the standard
strong Wolfe conditions*.

Uniform convexity alone does not guarantee global convergence of classical
DFP under the standard strong Wolfe conditions. The paper gives a
two-dimensional $C^2$ counterexample whose gradient norms converge to a
positive constant and whose iterates accumulate on a circle. It also proves
convergence under the standard weak Wolfe conditions for strongly convex
$C^2$ objective functions in two dimensions when the Hessian is locally
Lipschitz continuous near the initial level set.

The MATLAB code reproduces the numerical experiments in Sections 6.1–6.2:
the prescribed two-step recurrence, finite interpolation, and DFP/BFGS
comparisons with Wolfe line searches.

[中文说明](README_zh.md)

## Repository contents

| Paper part | Implementation | Output |
| --- | --- | --- |
| Section 6.1: geometry and asymptotics | `DFPExperiment.oracle`, `geometrySummary` | Figure 1(a) and Table 1 |
| Section 6.1: finite-function BFGS trajectory | `DFPExperiment.finite`, `run` | Figure 1(b) |
| Section 6.2: finite-function comparisons | `DFPExperiment.finite`, `run`, `wolfe_search` | Figure 2, Table 2, departure and line-search diagnostics |
| Finite-function Hessian bounds | `experiments/certify_bounds.m` | Eight exact rational certificates |

The two files in `src/` implement the recurrence, interpolation, updates,
and line searches. `run_experiments.m` calls the fixed paper protocol in
`experiments/`, which also exports the tables and plots. Tests are in `tests/`;
compact reference results and manuscript figures are in `reference/`.

## Installation

The reference experiments use MATLAB R2023a with the Statistics and Machine
Learning Toolbox. The Symbolic Math Toolbox is also required for the exact
certificates (`certify` and `all`). The optimization routines are implemented
directly; neither Python nor the Optimization Toolbox is required.

Download the repository or clone it:

```sh
git clone https://github.com/bqliu815/DFP-Counterexample.git
```

Open the repository root folder in MATLAB. Run the tests and a short example with:

```matlab
run_experiments('test');
run_experiments('smoke');
```

To run the DFP/BFGS comparison directly:

```matlab
addpath('src');
[objective, orbit] = DFPExperiment.finite(0.0025);
dfp = DFPExperiment.run(objective, orbit, 'dfp', 'zoom_unit', 5000);
bfgs = DFPExperiment.run(objective, orbit, 'bfgs', 'zoom_unit', 1000);
dfp.summary
bfgs.summary
```

## Reproducing the numerical experiments

The complete protocol generates the data, checks the finite-function bounds
and recorded steps, then exports the tables and figures:

```matlab
run_experiments('all');
```

For the single-threaded setup used in the paper, run from the repository root:

```sh
matlab -singleCompThread -batch "run_experiments('all')"
```

Results are written to `results/paper/`, with `raw/`, `tables/`, and `figures/`
subdirectories. A second argument selects a different output directory.
See [docs/REPRODUCIBILITY.md](docs/REPRODUCIBILITY.md) for individual stages,
experimental parameters, output files, and the current verification status.
The 14 tests and complete protocol have been run in MATLAB R2023a, including
the exact certificates and table and figure exports. The 28 trajectory CSV
files agree byte for byte with the reference MATLAB results.

## Citation

If you use this code, please cite the accompanying manuscript:

```bibtex
@unpublished{LiuWangWenYuanZhang2026DFP,
  author = {Benqi Liu and Zichen Wang and Zaiwen Wen and
            Yaxiang Yuan and Liwei Zhang},
  title  = {A counterexample to global convergence of classical {DFP}
            under the standard strong {Wolfe} conditions},
  year   = {2026},
  note   = {Manuscript}
}
```

Machine-readable citation metadata is provided in [CITATION.cff](CITATION.cff).

## Authors and contact

- Benqi Liu: [bqliu@pku.edu.cn](mailto:bqliu@pku.edu.cn)
- Zichen Wang: [zichenwang25@stu.pku.edu.cn](mailto:zichenwang25@stu.pku.edu.cn)
- Zaiwen Wen: [wenzw@pku.edu.cn](mailto:wenzw@pku.edu.cn)
- Yaxiang Yuan: [yyx@lsec.cc.ac.cn](mailto:yyx@lsec.cc.ac.cn)
- Liwei Zhang: [zhanglw@mail.neu.edu.cn](mailto:zhanglw@mail.neu.edu.cn)

## License

Original code is released under the [MIT License](LICENSE).
The SciPy adaptations in `src/wolfe_search.m` retain the BSD 3-Clause license;
see [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md) for attribution.

## Questions and contributions

For questions or bug reports, open an issue with the MATLAB release, command,
and full error message. Pull requests are welcome; describe the change and
run `run_experiments('test')` and `run_experiments('smoke')` before submitting.
