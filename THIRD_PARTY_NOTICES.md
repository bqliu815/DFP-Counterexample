# Third-party code

`src/wolfe_search.m` contains MATLAB adaptations of the strong Wolfe searches
in SciPy 1.17.0:

- [scipy/optimize/_linesearch.py](https://github.com/scipy/scipy/blob/v1.17.0/scipy/optimize/_linesearch.py): bracketing and cubic zoom.
- [scipy/optimize/_dcsrch.py](https://github.com/scipy/scipy/blob/v1.17.0/scipy/optimize/_dcsrch.py): More–Thuente search and safeguarded step updates.

The upstream copyright notice and BSD 3-Clause license are retained verbatim
in [licenses/SCIPY_LICENSE.txt](licenses/SCIPY_LICENSE.txt). The MATLAB file
is distributed under that license. The top-level MIT License applies to
original project code and does not replace the SciPy notice.

The upstream `_dcsrch.py` identifies the MINPACK-1 and MINPACK-2 origins:
Jorge J. Moré and David J. Thuente; Brett M. Averick, Richard G. Carter, and
Jorge J. Moré; Argonne National Laboratory and the University of Minnesota.
This provenance is retained here.

The MATLAB adaptation uses the paper's fixed Wolfe parameters and iteration
limits, records all trials, and offers unit-first and history-based first
trials. No SciPy or Python installation is required to run it.

MATLAB and its toolboxes are external dependencies and are not redistributed.
