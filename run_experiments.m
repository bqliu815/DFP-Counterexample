function outDir = run_experiments(mode, outDir)
% RUN_EXPERIMENTS Reproduce the polar DFP counterexample experiments.
%   Default: 'smoke'. Use 'all' for the complete paper protocol.
%   Other modes: 'test', 'wolfe', 'scales', 'certify', 'figures'.
%   'scales' varies beta in the prescribed recurrence without finite solvers.
% SPDX-License-Identifier: MIT
if nargin < 1, mode = 'smoke'; end
mode = validatestring(mode, {'smoke', 'all', 'test', 'wolfe', 'scales', 'certify', 'figures'});
repoDir = fileparts(mfilename('fullpath'));
if nargin < 2 || isempty(outDir)
    if strcmp(mode, 'smoke')
        outDir = fullfile(repoDir, 'results', 'smoke');
    else
        outDir = fullfile(repoDir, 'results', 'paper');
    end
end
oldPath = path;
restorePath = onCleanup(@()path(oldPath));
addpath(fullfile(repoDir, 'polar'));
if strcmp(mode, 'test')
    results = runtests(fullfile(repoDir, 'tests'));
    assertSuccess(results);
else
    run_polar_experiments(mode, outDir);
end
end
