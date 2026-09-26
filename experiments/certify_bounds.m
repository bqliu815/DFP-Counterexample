function certify_bounds(outDir)
    % CERTIFY_BOUNDS Exact rational checks for the stored finite functions.
    % sym(x,'f') interprets each binary64 coefficient as its exact dyadic value.
    % SPDX-License-Identifier: MIT
    rawDir = fullfile(outDir, 'raw');
    assert(isequal(sym(.1, 'f'), sym(3602879701896397) / sym(36028797018963968)));
    grid = [.001, .002, .0025];
    den = sym(10)^40;
    sqrt3up = (floor(sqrt(sym(3)) * den) + 1) / den;
    assert(isAlways(sqrt3up^2 > 3));
    K = sym(135) / 16 + sym(15) / 2 * sqrt3up;
    allcert = cell(numel(grid), 1);
    for i = 1:numel(grid)
        started = tic;
        e = grid(i);
        tag = strrep(sprintf('%.8g', e), '.', 'p');
        S = load(fullfile(rawDir, ['finite_', tag, '.mat']), 'obj');
        obj = S.obj;
        if e == .001
            ds = '0.17860';
        elseif e == .002
            ds = '0.35732';
        elseif e == .0025
            ds = '0.4467372';
        else
            ds = '0.5';
        end
        delta = sym(ds);
        c = sym(obj.center_ref', 'f');
        n = size(obj.points, 1);
        [~, order] = sort(obj.points(:, 2));
        chunks = ceil(n / 1000);
        last = -sym(Inf);
        for b = 1:chunks
            ix = (b - 1) * 1000 + 1:min(n, b * 1000);
            rr = sym(obj.radii(ix), 'f');
            aa = sym(obj.corrections(ix, :), 'f');
            cc = c - sym(obj.centers(ix, :), 'f');
            assert(all(isAlways(rr > 0)));
            assert(all(isAlways(K^2 * sum(aa.^2, 2) <= delta^2 * rr.^2)), ...
                   'Stored correction bound failed');
            assert(all(isAlways(K^2 * sum(cc.^2, 2) <= delta^2 * rr.^2)), ...
                   'Exact difference of centers bound failed');
            si = order(ix);
            yy = sym(obj.points(si, 2), 'f');
            sr = sym(obj.radii(si), 'f');
            left = yy - sr;
            right = yy + sr;
            assert(isAlways(left(1) > last));
            assert(all(isAlways(left(2:end) > right(1:end - 1))));
            last = right(end);
            if mod(b, 20) == 0 || b == chunks
                fprintf('CERT %.8g %d/%d\n', e, b, chunks);
            end
        end
        cert = struct('epsilon0', e, 'endpoints', n, 'delta_decimal', ds, ...
                      'certified_lower', char(1 - delta), 'certified_upper', char(1 + delta), ...
                      'lower_decimal', double(1 - delta), 'upper_decimal', double(1 + delta), ...
                      'max_correction_over_radius_float', obj.ratio, 'floating_hessian_band', obj.band, ...
                      'sqrt3_upper_exact', char(sqrt3up), 'hessian_constant_upper_exact', char(K), ...
                      'arithmetic', 'Symbolic Math Toolbox exact rational arithmetic', ...
                      'stored_binary64_coefficients_exact', true, 'all_projected_supports_disjoint', true, ...
                      'stored_corrections_checked', true, 'exact_center_differences_checked', true, ...
                      'within_half_three_halves', isAlways(delta <= sym(1) / 2), 'elapsed_seconds', toc(started));
        DFPExperiment.writeJSON(fullfile(rawDir, ['certificate_', tag, '.json']), cert);
        allcert{i} = cert;
        fprintf('CERTIFIED %s\n', jsonencode(cert));
    end
    DFPExperiment.writeJSON(fullfile(rawDir, 'certificates_complete.json'), allcert);
end
