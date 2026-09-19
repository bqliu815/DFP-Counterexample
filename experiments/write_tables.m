function write_tables(outDir)
    % WRITE_TABLES Export paper tables and the departure diagnostic as CSV files.
    % Exact bounds are read from the completed finite-function certificates.
    % SPDX-License-Identifier: MIT

    rawDir = fullfile(outDir, 'raw');
    tableDir = fullfile(outDir, 'tables');
    if ~isfolder(tableDir)
        mkdir(tableDir);
    end
    geometry = jsondecode(fileread(fullfile(rawDir, 'geometry.json')));
    quantity = {'epsilon_increment'; 'amplitude_increment'; 'rotation_increment'; ...
                'line_ratio_residual'; 'secant_residual'};
    computed = [geometry.normalized_medians(:); geometry.max_line_ratio_residual; ...
                geometry.max_secant_residual];
    limit = [-1.5; -6.5; -3; 0; 0];
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

    fit = jsondecode(fileread(fullfile(rawDir, 'departure_fit.json')));
    departures = table(fit.grid(:), fit.departure_indices(:), fit.support_indices(:), ...
                       'VariableNames', {'epsilon0', 'plateau_exit', 'support_exit'});
    writetable(departures, fullfile(tableDir, 'Departure.csv'));
end
