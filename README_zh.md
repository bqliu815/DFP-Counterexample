# DFP-Counterexample

本仓库配套论文

*A counterexample to global convergence of classical DFP under the standard
strong Wolfe conditions*。

在标准强 Wolfe 条件下，一致凸性并不能保证经典 DFP 方法的全局收敛。
论文给出一个二维 $C^2$ 反例，其梯度范数趋于正常数，迭代点的聚点集为一个圆。
另一方面，对于二维强凸 $C^2$ 目标函数，若 Hessian 在初始水平集的邻域内
局部 Lipschitz 连续，则标准弱 Wolfe 条件下的经典 DFP 方法收敛。

这里的 MATLAB 程序对应第 6.1–6.2 节，包括预定两步递推、有限插值、
Wolfe 线搜索，以及 DFP 与 BFGS 的数值比较。

[English](README.md)

## 仓库内容

| 论文内容 | 实现 | 输出 |
| --- | --- | --- |
| 第 6.1 节：几何与渐近行为 | `DFPExperiment.oracle`、`geometrySummary` | 图 1(a)、表 1 |
| 第 6.1 节：有限函数上的 BFGS 轨道 | `DFPExperiment.finite`、`run` | 图 1(b) |
| 第 6.2 节：有限函数上的比较 | `DFPExperiment.finite`、`run`、`wolfe_search` | 图 2、表 2、离轨与线搜索诊断 |
| 有限函数的 Hessian 界 | `experiments/certify_bounds.m` | 八份精确有理数证书 |

`src/` 保存递推与优化程序，`experiments/` 保存实验流程、表格与作图程序，
`tests/` 保存测试，`reference/` 保存少量参考数值及论文图件。
数值核心集中在 `src/` 的两个文件中，论文复现统一由 `run_experiments.m` 运行。

## 安装

参考实验使用 MATLAB R2023a 和 Statistics and Machine Learning Toolbox。
精确证书阶段（`certify`、`all`）还需 Symbolic Math Toolbox。优化方法直接
实现，无需 Python 或 Optimization Toolbox。

下载仓库，或使用以下命令克隆：

```sh
git clone https://github.com/bqliu815/DFP-Counterexample.git
```

将 MATLAB 当前文件夹切换到仓库根目录，执行测试和小规模示例：

```matlab
run_experiments('test');
run_experiments('smoke');
```

也可以直接调用核心程序，比较同一有限目标函数上的 DFP 与 BFGS：

```matlab
addpath('src');
[objective, orbit] = DFPExperiment.finite(0.0025);
dfp = DFPExperiment.run(objective, orbit, 'dfp', 'zoom_unit', 5000);
bfgs = DFPExperiment.run(objective, orbit, 'bfgs', 'zoom_unit', 1000);
dfp.summary
bfgs.summary
```

## 复现实验

完整流程生成数据、检查有限函数的 Hessian 界及记录的迭代步，再导出图表：

```matlab
run_experiments('all');
```

与论文一致的单计算线程运行方式为：

```sh
matlab -singleCompThread -batch "run_experiments('all')"
```

结果位于 `results/paper/`，其中 `raw/`、`tables/`、`figures/` 分别保存原始
数据、表格和图件。第二个参数可指定其他输出目录。各阶段命令、实验参数、
输出文件及验证状态见 [复现说明](docs/REPRODUCIBILITY.md)。14 项测试和完整
复现流程均已在 MATLAB R2023a 中运行，包括精确证书及图表导出；新生成的
28 份轨迹 CSV 文件与参考 MATLAB 结果逐字节一致。

## 引用与作者

论文作者为 Benqi Liu、Zichen Wang、Zaiwen Wen、Yaxiang Yuan 和 Liwei Zhang。
引用格式和联系方式见 [英文首页](README.md#citation)，机器可读信息见
[CITATION.cff](CITATION.cff)。

## 许可

原创代码采用 [MIT 许可](LICENSE)。`src/wolfe_search.m` 中转写自 SciPy 的代码
保留 BSD 3-Clause 许可与来源，详见 [第三方说明](THIRD_PARTY_NOTICES.md)。

## 问题反馈与贡献

欢迎提交 issue 或 pull request。报告问题时请附 MATLAB 版本、运行命令和
完整错误信息；提交代码修改前请运行 `test` 与 `smoke` 两个检查。
