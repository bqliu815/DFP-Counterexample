function [alpha, fp, gp, trials] = wolfe_search(evaluate, x, f0, g0, d, policy, previous)
    % WOLFE_SEARCH Weak Wolfe, cubic-zoom, and More--Thuente line searches.
    % Policies: weak_unit, zoom_unit, mt_unit, zoom_history, mt_history.
    % evaluate(x) returns [f,g]; previous is the preceding function value
    % for history initialization, or NaN for the first iteration.
    % Strong cubic zoom and More-Thuente follow SciPy 1.17.0 _linesearch.py and
    % _dcsrch.py. See ../THIRD_PARTY_NOTICES.md and ../licenses/SCIPY_LICENSE.txt.
    % SPDX-License-Identifier: BSD-3-Clause
    % Trials record [alpha, function value, directional derivative].
    policies = {'weak_unit', 'zoom_unit', 'mt_unit', 'zoom_history', 'mt_history'};
    if ~any(strcmp(policy, policies))
        error('dfp:UnknownPolicy', 'Unknown line-search policy: %s.', policy);
    end
    c1 = .25;
    c2 = .75;
    amax = 64;
    dg0 = g0' * d;
    assert(dg0 < 0);
    alpha = 1;
    trials = zeros(0, 3);
    fp = NaN;
    gp = NaN(2, 1);
    if contains(policy, 'history') && isfinite(previous)
        alpha = min(1, 1.01 * 2 * (f0 - previous) / dg0);
        if alpha < 0
            alpha = 1;
        end
    end
    if startsWith(policy, 'weak')
        lo = 0;
        flo = f0;
        hi = NaN;
        for it = 1:100
            [fp, gp, dg] = trial(alpha);
            if fp > f0 + c1 * alpha * dg0
                hi = alpha;
                break
            end
            if dg >= c2 * dg0
                return
            end
            lo = alpha;
            flo = fp;
            alpha = 2 * alpha;
            if alpha > amax
                hi = alpha;
                break
            end
        end
        assert(isfinite(hi), 'Failed to bracket weak Wolfe step');
        while size(trials, 1) < 100
            alpha = (lo + hi) / 2;
            [fp, gp, dg] = trial(alpha);
            if fp > f0 + c1 * alpha * dg0 || fp >= flo
                hi = alpha;
                continue
            end
            if dg >= c2 * dg0
                return
            end
            lo = alpha;
            flo = fp;
        end
        error('Weak Wolfe evaluation limit');
    elseif startsWith(policy, 'zoom')
        a0 = 0;
        fa0 = f0;
        dga0 = dg0;
        [fp, gp, dg] = trial(alpha);
        for it = 1:40
            assert(alpha > 0 && a0 <= amax, 'Invalid strong Wolfe trial');
            if fp > f0 + c1 * alpha * dg0 || (it > 1 && fp >= fa0)
                [alpha, fp, gp] = zoom(a0, alpha, fa0, fp, dga0);
                return
            end
            if abs(dg) <= -c2 * dg0
                return
            end
            if dg >= 0
                [alpha, fp, gp] = zoom(alpha, a0, fp, fa0, dg);
                return
            end
            a0 = alpha;
            fa0 = fp;
            dga0 = dg;
            alpha = min(2 * alpha, amax);
            [fp, gp, dg] = trial(alpha);
        end
        error('Strong Wolfe outer iteration limit');
    elseif startsWith(policy, 'mt')
        amin = 1e-12;
        xtol = 1e-14;
        assert(alpha >= amin && alpha <= amax);
        bracket = false;
        stage = 1;
        gtest = c1 * dg0;
        width = amax - amin;
        width1 = 2 * width;
        ax = 0;
        fx = f0;
        gx = dg0;
        ay = 0;
        fy = f0;
        gy = dg0;
        stmin = 0;
        stmax = 5 * alpha;
        for it = 1:99
            [fp, gp, dg] = trial(alpha);
            ftest = f0 + alpha * gtest;
            if stage == 1 && fp <= ftest && dg >= 0
                stage = 2;
            end
            if fp <= ftest && abs(dg) <= -c2 * dg0
                return
            end
            bad = (bracket && (alpha <= stmin || alpha >= stmax)) || ...
                (bracket && stmax - stmin <= xtol * stmax) || ...
                (alpha == amax && fp <= ftest && dg <= gtest) || ...
                (alpha == amin && (fp > ftest || dg >= gtest));
            assert(~bad, 'More-Thuente rounding or step-bound warning');
            if stage == 1 && fp <= fx && fp > ftest
                fm = fp - alpha * gtest;
                fxm = fx - ax * gtest;
                fym = fy - ay * gtest;
                [ax, fxm, gxm, ay, fym, gym, alpha, bracket] = dcstep( ...
                    ax, fxm, gx - gtest, ay, fym, gy - gtest, ...
                    alpha, fm, dg - gtest, bracket, stmin, stmax);
                fx = fxm + ax * gtest;
                fy = fym + ay * gtest;
                gx = gxm + gtest;
                gy = gym + gtest;
            else
                [ax, fx, gx, ay, fy, gy, alpha, bracket] = dcstep( ...
                    ax, fx, gx, ay, fy, gy, alpha, fp, dg, bracket, stmin, stmax);
            end
            assert(isreal(alpha) && isfinite(alpha), 'More-Thuente nonfinite trial');
            if bracket
                if abs(ay - ax) >= .66 * width1
                    alpha = ax + .5 * (ay - ax);
                end
                width1 = width;
                width = abs(ay - ax);
                stmin = min(ax, ay);
                stmax = max(ax, ay);
            else
                stmin = alpha + 1.1 * (alpha - ax);
                stmax = alpha + 4 * (alpha - ax);
            end
            alpha = min(amax, max(amin, alpha));
            if bracket && (alpha <= stmin || alpha >= stmax || stmax - stmin <= xtol * stmax)
                alpha = ax;
            end
        end
        error('More-Thuente iteration limit');
    else
        error('Unknown line search');
    end

    function [v, gv, dgv] = trial(a)
        [v, gv] = evaluate(x + a * d);
        dgv = gv' * d;
        trials(end + 1, :) = [a, v, dgv];
    end

    function [a, v, gv] = zoom(lo, hi, flo, fhi, dlo)
        arec = 0;
        frec = f0;
        for zi = 0:10
            da = hi - lo;
            lower = min(lo, hi);
            upper = max(lo, hi);
            a = NaN;
            if zi > 0
                a = cubicmin(lo, flo, dlo, hi, fhi, arec, frec);
            end
            if ~isfinite(a) || a > upper - .2 * da || a < lower + .2 * da
                db = hi - lo;
                B = (fhi - flo - dlo * db) / (db * db);
                a = lo - dlo / (2 * B);
                if ~isfinite(a) || a > upper - .1 * da || a < lower + .1 * da
                    a = lo + .5 * da;
                end
            end
            [v, gv, dv] = trial(a);
            if v > f0 + c1 * a * dg0 || v >= flo
                arec = hi;
                frec = fhi;
                hi = a;
                fhi = v;
            else
                if abs(dv) <= -c2 * dg0
                    return
                end
                if dv * (hi - lo) >= 0
                    arec = hi;
                    frec = fhi;
                    hi = lo;
                    fhi = flo;
                else
                    arec = lo;
                    frec = flo;
                end
                lo = a;
                flo = v;
                dlo = dv;
            end
        end
        error('Strong Wolfe zoom iteration limit');
    end

end

function a = cubicmin(x, fx, dx, y, fy, z, fz)
    db = y - x;
    dc = z - x;
    den = (db * dc)^2 * (db - dc);
    ab = [dc^2, -db^2; -dc^3, db^3] * [fy - fx - dx * db; fz - fx - dx * dc] / den;
    A = ab(1);
    B = ab(2);
    rad = B * B - 3 * A * dx;
    if rad < 0
        a = NaN;
    else
        a = x + (-B + sqrt(rad)) / (3 * A);
    end
    if ~isreal(a) || ~isfinite(a)
        a = NaN;
    end
end

function [ax, fx, dx, ay, fy, dy, ap, bracket] = dcstep(ax, fx, dx, ay, fy, dy, ap, fp, dp, bracket, amin, amax)
    sgnd = sign(dp) * sign(dx);
    if fp > fx
        theta = 3 * (fx - fp) / (ap - ax) + dx + dp;
        sc = max(abs([theta, dx, dp]));
        gamma = sc * sqrt((theta / sc)^2 - (dx / sc) * (dp / sc));
        if ap < ax
            gamma = -gamma;
        end
        r = ((gamma - dx) + theta) / (((gamma - dx) + gamma) + dp);
        ac = ax + r * (ap - ax);
        aq = ax + ((dx / ((fx - fp) / (ap - ax) + dx)) / 2) * (ap - ax);
        if abs(ac - ax) <= abs(aq - ax)
            af = ac;
        else
            af = ac + (aq - ac) / 2;
        end
        bracket = true;
    elseif sgnd < 0
        theta = 3 * (fx - fp) / (ap - ax) + dx + dp;
        sc = max(abs([theta, dx, dp]));
        gamma = sc * sqrt((theta / sc)^2 - (dx / sc) * (dp / sc));
        if ap > ax
            gamma = -gamma;
        end
        r = ((gamma - dp) + theta) / (((gamma - dp) + gamma) + dx);
        ac = ap + r * (ax - ap);
        aq = ap + (dp / (dp - dx)) * (ax - ap);
        if abs(ac - ap) > abs(aq - ap)
            af = ac;
        else
            af = aq;
        end
        bracket = true;
    elseif abs(dp) < abs(dx)
        theta = 3 * (fx - fp) / (ap - ax) + dx + dp;
        sc = max(abs([theta, dx, dp]));
        gamma = sc * sqrt(max(0, (theta / sc)^2 - (dx / sc) * (dp / sc)));
        if ap > ax
            gamma = -gamma;
        end
        r = ((gamma - dp) + theta) / ((gamma + (dx - dp)) + gamma);
        if r < 0 && gamma ~= 0
            ac = ap + r * (ax - ap);
        elseif ap > ax
            ac = amax;
        else
            ac = amin;
        end
        aq = ap + (dp / (dp - dx)) * (ax - ap);
        if bracket
            if abs(ac - ap) < abs(aq - ap)
                af = ac;
            else
                af = aq;
            end
            if ap > ax
                af = min(ap + .66 * (ay - ap), af);
            else
                af = max(ap + .66 * (ay - ap), af);
            end
        else
            if abs(ac - ap) > abs(aq - ap)
                af = ac;
            else
                af = aq;
            end
            af = min(amax, max(amin, af));
        end
    else
        if bracket
            theta = 3 * (fp - fy) / (ay - ap) + dy + dp;
            sc = max(abs([theta, dy, dp]));
            gamma = sc * sqrt((theta / sc)^2 - (dy / sc) * (dp / sc));
            if ap > ay
                gamma = -gamma;
            end
            r = ((gamma - dp) + theta) / (((gamma - dp) + gamma) + dy);
            af = ap + r * (ay - ap);
        elseif ap > ax
            af = amax;
        else
            af = amin;
        end
    end
    if fp > fx
        ay = ap;
        fy = fp;
        dy = dp;
    else
        if sgnd < 0
            ay = ax;
            fy = fx;
            dy = dx;
        end
        ax = ap;
        fx = fp;
        dx = dp;
    end
    ap = af;
end
