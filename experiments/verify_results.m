function verify_results(outDir)
% VERIFY_RESULTS Independently reevaluate the recorded fminunc iterates.
% SPDX-License-Identifier: MIT
rawDir = fullfile(outDir, 'raw');
names = {'geometry_bfgs'};
for e = [.001, .002, .0025]
    tag = strrep(sprintf('%.8g', e), '.', 'p');
    for prefix = {'dfp_', 'bfgs_', 'identity_dfp_', 'identity_bfgs_'}
        names{end+1} = [prefix{1}, tag];
    end
end
audits = cell(numel(names), 1);
for i = 1:numel(names)
    S = load(fullfile(rawDir, [names{i}, '.mat']), 'result');
    r = S.result;
    t = r.trace;
    tag = strrep(sprintf('%.8g', r.summary.epsilon0), '.', 'p');
    S = load(fullfile(rawDir, ['finite_', tag, '.mat']), 'obj');
    obj = S.obj;
    obj.tree = KDTreeSearcher(obj.points);
    assert(strcmp(r.summary.solver, 'fminunc'));
    assert(r.summary.function_evaluations == r.summary.counted_evaluations);
    assert(all(isfinite(t{:, 1:10}), 'all'));
    assert(height(t) == r.summary.recorded_steps);
    assert(all(diff(t.iteration) > 0));
    for k = 1:height(t)
        [f, g] = DFPExperiment.valueGrad(obj, [t.x1(k); t.x2(k)]);
        assert(abs(f - t.function_value(k)) < 1e-13 * max(1, abs(f)));
        assert(norm(g - [t.g1(k); t.g2(k)]) < 1e-13);
        assert(abs(norm(g) - t.gradient_norm(k)) < 1e-13);
    end
    [~, finalG] = DFPExperiment.valueGrad(obj, r.x);
    assert(abs(norm(finalG) - r.summary.final_gradient_norm) < 1e-13);
    assert(strcmp(r.summary.status, 'gradient_tolerance') == ...
           (norm(finalG) <= r.summary.gradient_tolerance));
    % Native stopping flags are recorded separately from the gradient test.
    audits{i} = struct('run', names{i}, 'recorded_steps', height(t), ...
        'exitflag', r.exitflag, 'status', r.summary.status, 'passed', true);
end
g = jsondecode(fileread(fullfile(rawDir, 'geometry.json')));
assert(g.armijo_failures == 0 && g.strong_failures == 0);
assert(g.max_secant_residual < 1e-12 && g.min_metric_eigenvalue > 0);
cert = jsondecode(fileread(fullfile(rawDir, 'certificates_complete.json')));
assert(numel(cert) == 3);
for i = 1:numel(cert)
    assert(cert(i).within_half_three_halves && cert(i).all_projected_supports_disjoint);
end
write_tables(outDir);
DFPExperiment.writeJSON(fullfile(rawDir, 'verification.json'), ...
    struct('passed', true, 'runs', {audits}, 'certificates_checked', 3, ...
    'claim_boundary', 'Stored finite objective functions and recorded fminunc iterates.'));
end
