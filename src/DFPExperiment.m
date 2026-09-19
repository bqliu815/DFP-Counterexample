classdef DFPExperiment
    % DFP recurrence, finite objective functions, and classical DFP/BFGS in 2D.
    % SPDX-License-Identifier: MIT
    methods (Static)

        function H = update(H, s, y, method)
            % UPDATE Apply one classical inverse DFP or BFGS update.
            % H is 2-by-2; s and y are 2-by-1 with s'*y > 0.
            % Symmetrization removes floating-point asymmetry only.
            sy = s' * y;
            assert(isfinite(sy) && sy > 0, 'Nonpositive curvature');
            if strcmp(method, 'dfp')
                Hy = H * y;
                yHy = y' * Hy;
                assert(yHy > 0, 'Nonpositive DFP denominator');
                H = H - (Hy * Hy') / yHy + (s * s') / sy;
            elseif strcmp(method, 'bfgs')
                V = eye(2) - (s * y') / sy;
                H = V * H * V' + (s * s') / sy;
            else
                error('dfp:UnknownMethod', 'Method must be dfp or bfgs.');
            end
            H = (H + H') / 2;
        end

        function [epsilon, G, Q, h, p] = state(H, g)
            % STATE Recover the positive eigenframe coordinates of a cycle.
            % The first column of Q is the low-eigenvalue direction oriented toward g.
            [V, D] = eig(H, 'vector');
            [lambda, idx] = sort(D);
            low = V(:, idx(1));
            if low' * g < 0
                low = -low;
            end
            high = [-low(2); low(1)];
            Q = [low, high];
            gc = Q' * g;
            G = gc(1);
            u = gc(2) / G;
            h = high' * H * high;
            r = lambda(1) / (h * u);
            p = u / r;
            assert(min([G, h, r, p]) > 0, 'Left positive coordinate chart');
            epsilon = sqrt(r);
        end

        function o = oracle(cycles, e0)
            % ORACLE Generate the prescribed two-step DFP recurrence.
            % Step lengths are prescribed, without a line search.
            % The initial p and h are the fourth-order truncations used in the paper.
            r = e0^2;
            p = 2 + (198 / 5) * e0^3 - (9 / 5) * e0^4;
            h = 1 + 8 * e0^3;
            H = diag([h * p * r * r, h]);
            g = [1; p * r];
            x = g;
            o.h0 = H;
            n = 2 * cycles;
            o.x = zeros(n + 1, 2);
            o.g = o.x;
            o.centers = o.x;
            o.x(1, :) = x';
            o.g(1, :) = g';
            o.alpha = zeros(n, 1);
            o.q = o.alpha;
            o.tau = o.alpha;
            o.line_ratio = o.alpha;
            o.secant_residual = o.alpha;
            o.steps = zeros(n, 2);
            o.strong_ratio = o.alpha;
            o.cycle_epsilon = zeros(cycles + 1, 1);
            o.cycle_amplitude = o.cycle_epsilon;
            angle = o.cycle_epsilon;
            o.min_metric_eigenvalue = Inf;
            k = 0;
            for j = 1:cycles
                [e, G, Q] = DFPExperiment.state(H, g);
                o.cycle_epsilon(j) = e;
                o.cycle_amplitude(j) = G;
                angle(j) = atan2(Q(2, 1), Q(1, 1));
                for leg = 1:2
                    [~, ~, Q] = DFPExperiment.state(H, g);
                    if leg == 1
                        mu = e;
                        tau = 2 / 3;
                    else
                        mu = -2 * e;
                        tau = 1 / 3;
                    end
                    A = Q * [1, mu; mu, 1] * Q';
                    Hg = H * g;
                    alpha = tau * (g' * Hg) / (Hg' * A * Hg);
                    s = -alpha * Hg;
                    y = A * s;
                    q = -g' * s;
                    gp = g + y;
                    xp = x + s;
                    Hp = DFPExperiment.update(H, s, y, 'dfp');
                    k = k + 1;
                    o.alpha(k) = alpha;
                    o.q(k) = q;
                    o.tau(k) = tau;
                    o.line_ratio(k) = (s' * y) / q;
                    o.steps(k, :) = s';
                    o.secant_residual(k) = norm(Hp * y - s) / norm(s);
                    o.strong_ratio(k) = abs(gp' * s) / abs(g' * s);
                    o.min_metric_eigenvalue = min(o.min_metric_eigenvalue, min(eig(Hp)));
                    H = Hp;
                    g = gp;
                    x = xp;
                    o.x(k + 1, :) = x';
                    o.g(k + 1, :) = g';
                    o.centers(k + 1, :) = (x - g)';
                end
            end
            [e, G, Q] = DFPExperiment.state(H, g);
            o.cycle_epsilon(end) = e;
            o.cycle_amplitude(end) = G;
            angle(end) = atan2(Q(2, 1), Q(1, 1));
            o.cycle_angle = unwrap(angle);
            o.hfinal = H;
            o.epsilon0 = e0;
            o.cycles = cycles;
        end

        function [obj, o] = finite(e0, prefixFactor)
            % FINITE Realize a finite prefix by disjoint local gradient corrections.
            % obj.band is a floating-point bound; certify_bounds checks exact bounds.
            if nargin < 2
                prefixFactor = 2;
            end
            budget = max(100, ceil(.5 * e0^(-1.5)));
            steps = ceil(prefixFactor * budget) + 4;
            steps = steps + mod(steps, 2);
            o = DFPExperiment.oracle(steps / 2, e0);
            obj.points = o.x;
            obj.centers = o.centers;
            obj.center_ref = o.centers(end, :)';
            obj.corrections = obj.center_ref' - obj.centers;
            obj.tree = KDTreeSearcher(obj.points);
            [~, dist] = knnsearch(obj.tree, obj.points, 'K', 2);
            obj.radii = .24 * dist(:, 2);
            assert(all(obj.radii > 0));
            obj.ratio = max(vecnorm(obj.corrections, 2, 2) ./ obj.radii);
            obj.perturbation = (135 / 16 + 15 * sqrt(3) / 2) * obj.ratio;
            obj.band = [1 - obj.perturbation, 1 + obj.perturbation];
            obj.epsilon0 = e0;
            obj.steps = steps;
            obj.budget = budget;
            obj.prefixFactor = prefixFactor;
        end

        function [f, g, nearest, dist] = valueGrad(obj, x)
            % VALUEGRAD Evaluate the finite objective function and its analytic gradient.
            % The radial cutoff is one up to rho/3 and zero from rho onward.
            d = x - obj.center_ref;
            f = .5 * (d' * d);
            g = d;
            [nearest, dist] = knnsearch(obj.tree, x');
            rho = obj.radii(nearest);
            if dist >= rho
                return
            end
            v = x - obj.points(nearest, :)';
            a = obj.corrections(nearest, :)';
            t = dist / rho;
            if t <= 1 / 3
                phi = 1;
                dp = 0;
            else
                width = 1 - 1 / 3;
                z = (t - 1 / 3) / width;
                phi = 1 - (10 * z^3 - 15 * z^4 + 6 * z^5);
                dp = -30 * z^2 * (1 - z)^2 / width;
            end
            av = a' * v;
            f = f + phi * av;
            g = g + phi * a;
            if dist > 0 && dp ~= 0
                g = g + av * dp * v / (rho * dist);
            end
        end

        function result = run(obj, o, method, policy, maxit, stopAtSupport)
            % RUN Apply classical DFP/BFGS to a fixed finite objective function.
            % Initial x and H come from o. Accepted steps and every line-search trial
            % are retained. The gradient tolerance is 1e-10, as in the paper.
            if nargin < 6
                stopAtSupport = false;
            end
            x = o.x(1, :)';
            H = o.h0;
            [f, g] = DFPExperiment.valueGrad(obj, x);
            prev = NaN;
            columns = {'iteration', 'alpha', 'function_value', 'gradient_norm', 'target_ratio', ...
                       'nearest_ratio', 'armijo_ratio', 'strong_ratio', 'weak_margin', 'curvature', ...
                       'secant_residual', 'metric_min', 'metric_max', 'evaluations', 'x1', 'x2', 'first_trial'};
            rows = zeros(maxit, numel(columns));
            events = cell(maxit, 1);
            status = 'iteration_limit';
            kdone = 0;
            for k = 1:maxit
                if norm(g) <= 1e-10
                    status = 'gradient_tolerance';
                    break
                end
                d = -H * g;
                dg = g' * d;
                if ~isfinite(dg) || dg >= 0
                    status = 'non_descent';
                    break
                end
                evaluator = @(z)DFPExperiment.valueGrad(obj, z);
                try
                    [alpha, ~, ~, trials] = wolfe_search(evaluator, x, f, g, d, policy, prev);
                catch err
                    status = ['line_search_failure: ', err.message];
                    break
                end
                s = alpha * d;
                xp = x + s;
                [fp, gp, idx, dist] = DFPExperiment.valueGrad(obj, xp);
                y = gp - g;
                if s' * y <= 0
                    status = 'nonpositive_curvature';
                    break
                end
                Hp = DFPExperiment.update(H, s, y, method);
                ev = eig(Hp);
                if min(ev) <= 0 || any(~isfinite(ev))
                    status = 'nonpositive_metric';
                    break
                end
                target = min(k + 1, size(o.x, 1));
                er = norm(xp - o.x(target, :)') / obj.radii(target);
                rows(k, :) = [k, alpha, fp, norm(gp), er, dist / obj.radii(idx), ...
                              (f - fp) / (-alpha * dg), abs(gp' * d) / abs(dg), (gp' * d - .75 * dg) / abs(dg), ...
                              s' * y, norm(Hp * y - s) / norm(s), min(ev), max(ev), size(trials, 1), xp', trials(1, 1)];
                events{k} = trials;
                kdone = k;
                prev = f;
                x = xp;
                f = fp;
                g = gp;
                H = Hp;
                if norm(g) <= 1e-10
                    status = 'gradient_tolerance';
                    break
                end
                if stopAtSupport && er > 1
                    status = 'support_exit';
                    break
                end
            end
            rows = rows(1:kdone, :);
            result.trace = array2table(rows, 'VariableNames', columns);
            result.trials = events(1:kdone);
            result.x0 = o.x(1, :);
            result.h0 = o.h0;
            result.hfinal = H;
            summary = struct('epsilon0', obj.epsilon0, 'method', method, 'policy', policy, ...
                             'status', status, 'iterations', kdone, ...
                             'prefix_steps', obj.steps, 'hessian_band', obj.band, ...
                             'initial_gradient_norm', norm(o.g(1, :)), 'final_gradient_norm', norm(g));
            if kdone > 0
                summary.plateau_exit = DFPExperiment.first(rows(:, 5) > 1 / 3);
                summary.support_exit = DFPExperiment.first(rows(:, 5) > 1);
                summary.first_nonunit = DFPExperiment.first(abs(rows(:, 2) - 1) > 1e-13);
                summary.armijo_failures = sum(rows(:, 7) < .25 - 1e-11);
                summary.strong_failures = sum(rows(:, 8) > .75 + 1e-11);
                summary.weak_failures = sum(rows(:, 9) < -1e-11);
                summary.min_armijo_ratio = min(rows(:, 7));
                summary.max_strong_ratio = max(rows(:, 8));
                summary.max_secant_residual = max(rows(:, 11));
                summary.min_metric_eigenvalue = min(rows(:, 12));
                summary.min_curvature = min(rows(:, 10));
                summary.alpha_range = [min(rows(:, 2)), max(rows(:, 2))];
                summary.third_alpha = NaN;
                if kdone >= 3
                    summary.third_alpha = rows(3, 2);
                end
            end
            result.summary = summary;
        end

        function j = first(mask)
            j = find(mask, 1);
            if isempty(j)
                j = NaN;
            end
        end

        function s = geometrySummary(o)
            % GEOMETRYSUMMARY Report recurrence residuals and asymptotic diagnostics.
            % The limiting center is estimated by the last computed reference center.
            e = o.cycle_epsilon;
            G = o.cycle_amplitude;
            theta = o.cycle_angle;
            tail = max(1, numel(e) - 10000):numel(e) - 1;
            vals = [diff(e) ./ e(1:end - 1).^4, (G(2:end) ./ G(1:end - 1) - 1) ./ e(1:end - 1).^4, ...
                    diff(theta) ./ e(1:end - 1).^2];
            fvals = .5 * sum((o.x - o.centers(end, :)).^2, 2);
            armijo = (fvals(1:end - 1) - fvals(2:end)) ./ o.q;
            s = struct('epsilon0', o.epsilon0, 'cycles', o.cycles, ...
                       'turns', abs(theta(end) - theta(1)) / (2 * pi), ...
                       'initial_gradient_norm', norm(o.g(1, :)), 'final_gradient_norm', norm(o.g(end, :)), ...
                       'limiting_radius_estimate', G(end) * exp(-13 * e(end) / 3), ...
                       'normalized_medians', median(vals(tail, :), 1), ...
                       'max_line_ratio_residual', max(abs(o.line_ratio - o.tau)), ...
                       'max_secant_residual', max(o.secant_residual), ...
                       'strong_ratio_medians', [median(o.strong_ratio(1:2:end)), median(o.strong_ratio(2:2:end))], ...
                       'min_armijo_ratio', min(armijo), 'armijo_failures', sum(armijo < .25 - 1e-12), ...
                       'strong_failures', sum(o.strong_ratio > .75 + 1e-12), ...
                       'min_metric_eigenvalue', o.min_metric_eigenvalue);
        end

        function writeJSON(path, value)
            fid = fopen(path, 'w');
            assert(fid >= 0);
            clean = onCleanup(@()fclose(fid));
            fprintf(fid, '%s\n', jsonencode(value, 'PrettyPrint', true));
        end

        function saveRun(stem, result)
            save([stem, '.mat'], 'result', '-v7.3');
            writetable(result.trace, [stem, '.csv']);
            DFPExperiment.writeJSON([stem, '.json'], result.summary);
            fprintf('%s\n', jsonencode(result.summary));
        end

    end
end
