function outDir = run_experiments(mode, outDir)
% RUN_EXPERIMENTS Reproduce the polar DFP counterexample experiments.
%   'paper' or 'all': reproduce Section 5, Figures 2 and 3, and the certificate.
%   'figures' and 'certify' reuse the Section 5 data. Default: 'smoke'.
%   Earlier experiments: 'legacy', 'legacy-figures', 'legacy-certify', 'scales', 'wolfe'.
% SPDX-License-Identifier: MIT
if nargin < 1, mode = 'smoke'; end
mode = validatestring(mode, {'smoke', 'paper', 'all', 'test', 'certify', 'figures', ...
    'legacy', 'legacy-figures', 'legacy-certify', 'wolfe', 'scales'});
repoDir = fileparts(mfilename('fullpath'));
if nargin < 2 || isempty(outDir)
    if strcmp(mode, 'smoke')
        outDir = fullfile(repoDir, 'results', 'smoke');
    elseif ismember(mode, {'paper', 'all', 'figures', 'certify'})
        outDir = fullfile(repoDir, 'results', 'section6');
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
