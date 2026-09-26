function write_tables(outDir)
    % WRITE_TABLES Export paper tables as CSV files.
    % Exact bounds are read from the completed finite-function certificates.
    % SPDX-License-Identifier: MIT

    rawDir = fullfile(outDir, 'raw');
    tableDir = fullfile(outDir, 'tables');
    if ~isfolder(tableDir)
        mkdir(tableDir);
    end
    geometry = jsondecode(fileread(fullfile(rawDir, 'geometry.json')));
    quantity = {'epsilon_increment'; 'amplitude_increment'; 'rotation_increment'; ...
                'strong_wolfe_ratio_first'; 'strong_wolfe_ratio_second'; ...
                'line_ratio_residual'; 'secant_residual'; ...
                'armijo_failures'; 'strong_wolfe_failures'};
    computed = [geometry.normalized_medians(:); geometry.strong_ratio_medians(:); ...
                geometry.max_line_ratio_residual; geometry.max_secant_residual; ...
                geometry.armijo_failures; geometry.strong_failures];
    limit = [-1.5; -6.5; -3; 1/3; 2/3; 0; 0; 0; 0];
    writetable(table(quantity, computed, limit), fullfile(tableDir, 'Table1.csv'));

    epsilon0 = [.001; .002; .0025];
    rows = zeros(3, 7);
    for k = 1:numel(epsilon0)
        tag = strrep(sprintf('%.8g', epsilon0(k)), '.', 'p');
        d = jsondecode(fileread(fullfile(rawDir, ['dfp_', tag, '.json'])));
        b = jsondecode(fileread(fullfile(rawDir, ['bfgs_', tag, '.json'])));
        certificate = jsondecode(fileread(fullfile(rawDir, ['certificate_', tag, '.json'])));
        rows(k, :) = [epsilon0(k), certificate.lower_decimal, certificate.upper_decimal, ...
                      d.iterations, d.final_gradient_norm, b.iterations, b.final_gradient_norm];
    end
    names = {'epsilon0', 'hessian_lower', 'hessian_upper', 'dfp_iterations', ...
             'dfp_final_gradient_norm', 'bfgs_iterations', 'bfgs_final_gradient_norm'};
    writetable(array2table(rows, 'VariableNames', names), fullfile(tableDir, 'Table2.csv'));

    rows = zeros(6, 5);
    k = 0;
    for e = epsilon0'
        tag = strrep(sprintf('%.8g', e), '.', 'p');
        for method = {'dfp', 'bfgs'}
            k = k + 1;
            s = jsondecode(fileread(fullfile(rawDir, ['identity_', method{1}, '_', tag, '.json'])));
            rows(k, :) = [e, strcmp(method{1}, 'bfgs'), s.iterations, ...
                          s.final_gradient_norm, s.exitflag];
        end
    end
    writetable(array2table(rows, 'VariableNames', ...
        {'epsilon0', 'is_bfgs', 'iterations', 'final_gradient_norm', 'exitflag'}), ...
        fullfile(tableDir, 'Identity.csv'));
end
