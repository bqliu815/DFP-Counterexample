function verify_results(outDir)
    % VERIFY_RESULTS Check accepted steps, agreement records, and certificates.
    % SPDX-License-Identifier: MIT
    rawDir = fullfile(outDir, 'raw');
    names = expectedRuns();
    audits = {};
    for i = 1:numel(names)
        path = fullfile(rawDir, names{i});
        assert(isfile(path), 'dfp:MissingInput', 'Missing trajectory: %s', names{i});
        S = load(path, 'result');
        r = S.result;
        t = r.trace;
        s = r.summary;
        assert(height(t) == s.iterations && height(t) > 0);
        assert(all(isfinite(t{:, :}), 'all'));
        assert(all(t.curvature > 0) && all(t.metric_min > 0));
        assert(all(t.armijo_ratio >= .25 - 1e-11));
        assert(all(t.weak_margin >= -1e-11));
        if ~startsWith(s.policy, 'weak')
            assert(all(t.strong_ratio <= .75 + 1e-11));
        end
        assert(all(t.secant_residual < 1e-8));
        assert(~contains(s.status, 'failure') && ~contains(s.status, 'nonpositive') ...
               && ~strcmp(s.status, 'non_descent'));
        if strcmp(s.method, 'bfgs')
            assert(strcmp(s.status, 'gradient_tolerance'));
        end
        audits{end + 1} = struct('file', names{i}, 'steps', height(t), 'passed', true, ...
                                 'strong_wolfe_failures', sum(t.strong_ratio > .75 + 1e-11));
    end
    g = jsondecode(fileread(fullfile(rawDir, 'geometry.json')));
    assert(g.armijo_failures == 0 && g.strong_failures == 0 && g.min_metric_eigenvalue > 0);
    assert(g.max_secant_residual < 1e-12);
    agreement = jsondecode(fileread(fullfile(rawDir, 'agreements.json')));
    for i = 1:numel(agreement)
        if iscell(agreement)
            a = agreement{i};
        else
            a = agreement(i);
        end
        if isfield(a, 'bfgs_weak_strong_identical')
            assert(a.bfgs_weak_strong_identical);
        end
        if isfield(a, 'unit_trajectory_identical')
            assert(a.unit_trajectory_identical);
        end
    end
    cert = jsondecode(fileread(fullfile(rawDir, 'certificates_complete.json')));
    assert(numel(cert) == 8);
    for i = 1:numel(cert)
        assert(cert(i).within_half_three_halves && cert(i).all_projected_supports_disjoint);
    end
    write_tables(outDir);
    DFPExperiment.writeJSON(fullfile(rawDir, 'verification.json'), ...
        struct('passed', true, 'run_audits', {audits}, ...
               'geometry_passed', true, 'eight_exact_certificates_passed', true, ...
               'claim_boundary', 'Finite computed trajectories and fixed stored objective functions only.'));
end

function names = expectedRuns()
    names = {'geometry_bfgs.mat'};
    grid = [.0005, .00065, .0008, .001, .00125, .0015, .002, .0025];
    for e = grid
        tag = strrep(sprintf('%.8g', e), '.', 'p');
        names{end + 1} = ['departure_', tag, '.mat'];
        if any(e == [.001, .002, .0025])
            for prefix = {'dfp_', 'bfgs_', 'bfgs_weak_'}
                names{end + 1} = [prefix{1}, tag, '.mat'];
            end
        end
        if any(e == [.001, .002])
            policies = {'weak_unit', 'zoom_unit', 'mt_unit', 'zoom_history', 'mt_history'};
            for j = 1:numel(policies)
                names{end + 1} = ['diagnostic_', tag, '_', policies{j}, '.mat'];
            end
        end
    end
end
