classdef DFPExperiment
    % Prescribed DFP recurrence and finite objective functions in two dimensions.
    % SPDX-License-Identifier: MIT
    methods (Static)

        function H = update(H, s, y)
            % UPDATE Classical inverse DFP update for the prescribed recurrence.
            sy = s' * y;
            Hy = H * y;
            yHy = y' * Hy;
            assert(isfinite(sy) && sy > 0 && yHy > 0, 'Nonpositive curvature');
            H = H - (Hy * Hy') / yHy + (s * s') / sy;
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
                    Hp = DFPExperiment.update(H, s, y);
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
