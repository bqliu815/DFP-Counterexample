function run_suite(mode, outDir)
% RUN_SUITE Reproduce the prescribed recurrence and fminunc comparisons.
% SPDX-License-Identifier: MIT
rawDir = fullfile(outDir, 'raw');
if ~isfolder(rawDir), mkdir(rawDir); end
started = tic;
switch mode
    case 'smoke'
        o = DFPExperiment.oracle(100, .03);
        s = DFPExperiment.geometrySummary(o);
        assert(s.max_secant_residual < 1e-12 && s.min_metric_eigenvalue > 0);
        [obj, o] = DFPExperiment.finite(.0025);
        b = run_fminunc(obj, o, 'bfgs');
        DFPExperiment.saveRun(fullfile(rawDir, 'smoke_bfgs'), b);
        assert(b.summary.final_gradient_norm < b.summary.initial_gradient_norm);
        DFPExperiment.writeJSON(fullfile(rawDir, 'smoke.json'), ...
            struct('passed', true, 'geometry', s, 'bfgs', b.summary));
    case 'full'
        o = DFPExperiment.oracle(100000, .03);
        save(fullfile(rawDir, 'geometry.mat'), 'o', '-v7');
        DFPExperiment.writeJSON(fullfile(rawDir, 'geometry.json'), ...
            DFPExperiment.geometrySummary(o));
        [obj, o] = DFPExperiment.finite(.03);
        b = run_fminunc(obj, o, 'bfgs');
        DFPExperiment.saveRun(fullfile(rawDir, 'geometry_bfgs'), b);
        saveObjective(fullfile(rawDir, 'finite_0p03.mat'), obj, o);
        summaries = {};
        for e = [.001, .002, .0025]
            tag = strrep(sprintf('%.8g', e), '.', 'p');
            [obj, o] = DFPExperiment.finite(e);
            saveObjective(fullfile(rawDir, ['finite_', tag, '.mat']), obj, o);
            for initialization = {'prescribed', 'identity'}
                for method = {'dfp', 'bfgs'}
                    result = run_fminunc(obj, o, method{1}, 5000, initialization{1});
                    prefix = '';
                    if strcmp(initialization{1}, 'identity'), prefix = 'identity_'; end
                    DFPExperiment.saveRun(fullfile(rawDir, ...
                        [prefix, method{1}, '_', tag]), result);
                    summaries{end+1} = result.summary;
                end
            end
        end
        DFPExperiment.writeJSON(fullfile(rawDir, 'full_complete.json'), ...
            struct('complete', true, 'elapsed_seconds', toc(started), 'runs', {summaries}));
    case 'certify'
        certify_bounds(outDir);
    case 'figures'
        make_figures(outDir);
    case 'verify'
        verify_results(outDir);
    otherwise
        error('Unknown suite mode');
end
fprintf('MODE %s FINISHED in %.3f seconds\n', mode, toc(started));
end

function saveObjective(path, obj, o)
obj = rmfield(obj, 'tree');
save(path, 'obj', 'o', '-v7');
end
