function run_suite(mode, outDir)
    % RUN_SUITE Fixed experimental protocol for Section 6 of the paper.
    % SPDX-License-Identifier: MIT
    rawDir = fullfile(outDir, 'raw');
    if ~isfolder(rawDir)
        mkdir(rawDir);
    end
    started = tic;
    switch mode
        case 'smoke'
            o = DFPExperiment.oracle(100, .03);
            s = DFPExperiment.geometrySummary(o);
            assert(s.armijo_failures == 0 && s.strong_failures == 0);
            assert(s.max_secant_residual < 1e-12 && s.min_metric_eigenvalue > 0);
            [obj, o] = DFPExperiment.finite(.0025);
            r = DFPExperiment.run(obj, o, 'bfgs', 'zoom_unit', 1000);
            DFPExperiment.saveRun(fullfile(rawDir, 'smoke_bfgs'), r);
            assert(strcmp(r.summary.status, 'gradient_tolerance'));
            assert(r.summary.armijo_failures == 0 && r.summary.strong_failures == 0);
            DFPExperiment.writeJSON(fullfile(rawDir, 'smoke.json'), ...
                struct('passed', true, 'geometry', s, ...
                       'bfgs', r.summary, 'elapsed_seconds', toc(started)));
        case 'full'
            o = DFPExperiment.oracle(100000, .03);
            save(fullfile(rawDir, 'geometry.mat'), 'o', '-v7');
            s = DFPExperiment.geometrySummary(o);
            DFPExperiment.writeJSON(fullfile(rawDir, 'geometry.json'), s);
            fprintf('GEOMETRY %s\n', jsonencode(s));
            [obj, o] = DFPExperiment.finite(.03);
            r = DFPExperiment.run(obj, o, 'bfgs', 'zoom_unit', 1000);
            DFPExperiment.saveRun(fullfile(rawDir, 'geometry_bfgs'), r);
            saveObjective(fullfile(rawDir, 'finite_0p03.mat'), obj, o);
            grid = [.0005, .00065, .0008, .001, .00125, .0015, .002, .0025];
            departures = zeros(numel(grid), 3);
            summaries = {};
            diagnostics = {};
            agreements = {};
            for i = 1:numel(grid)
                e = grid(i);
                tag = strrep(sprintf('%.8g', e), '.', 'p');
                [obj, o] = DFPExperiment.finite(e);
                saveObjective(fullfile(rawDir, ['finite_', tag, '.mat']), obj, o);
                r = DFPExperiment.run(obj, o, 'dfp', 'weak_unit', obj.budget, true);
                DFPExperiment.saveRun(fullfile(rawDir, ['departure_', tag]), r);
                departures(i, :) = [e, r.summary.plateau_exit, r.summary.support_exit];
                summaries{end + 1} = r.summary;
                if any(e == [.001, .002, .0025])
                    r = DFPExperiment.run(obj, o, 'dfp', 'zoom_unit', 5000);
                    DFPExperiment.saveRun(fullfile(rawDir, ['dfp_', tag]), r);
                    rb = DFPExperiment.run(obj, o, 'bfgs', 'zoom_unit', 1000);
                    DFPExperiment.saveRun(fullfile(rawDir, ['bfgs_', tag]), rb);
                    rw = DFPExperiment.run(obj, o, 'bfgs', 'weak_unit', 1000);
                    DFPExperiment.saveRun(fullfile(rawDir, ['bfgs_weak_', tag]), rw);
                    same = isequal(rb.trace(:, {'alpha', 'x1', 'x2'}), rw.trace(:, {'alpha', 'x1', 'x2'}));
                    agreements{end + 1} = struct('epsilon0', e, 'bfgs_weak_strong_identical', same);
                end
                if any(e == [.001, .002])
                    n = min(obj.steps - 2, max(100, ceil(.8 * e^(-1.5))));
                    policies = {'weak_unit', 'zoom_unit', 'mt_unit', 'zoom_history', 'mt_history'};
                    reference = [];
                    for j = 1:numel(policies)
                        rr = DFPExperiment.run(obj, o, 'dfp', policies{j}, n);
                        DFPExperiment.saveRun(fullfile(rawDir, ['diagnostic_', tag, '_', policies{j}]), rr);
                        if j == 1
                            reference = rr.trace(:, {'alpha', 'x1', 'x2'});
                        elseif j <= 3
                            agreements{end + 1} = struct('epsilon0', e, 'policy', policies{j}, ...
                                'unit_trajectory_identical', ...
                                isequal(reference, rr.trace(:, {'alpha', 'x1', 'x2'})));
                        end
                        diagnostics{end + 1} = rr.summary;
                    end
                end
            end
            assert(all(isfinite(departures(:, 2))), 'Missing departure index');
            coef = polyfit(log(departures(:, 1)), log(departures(:, 2)), 1);
            observed = log(departures(:, 2));
            fitted = polyval(coef, log(departures(:, 1)));
            fit = struct('exponent', coef(1), 'log_prefactor', coef(2), ...
                         'r_squared', 1 - sum((observed - fitted).^2) / sum((observed - mean(observed)).^2), ...
                         'grid', grid, 'departure_indices', departures(:, 2)', 'support_indices', departures(:, 3)');
            DFPExperiment.writeJSON(fullfile(rawDir, 'departure_fit.json'), fit);
            DFPExperiment.writeJSON(fullfile(rawDir, 'agreements.json'), agreements);
            DFPExperiment.writeJSON(fullfile(rawDir, 'full_complete.json'), ...
                struct('complete', true, 'elapsed_seconds', toc(started), ...
                       'departure_runs', {summaries}, 'diagnostics', {diagnostics}, 'fit', fit));
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
