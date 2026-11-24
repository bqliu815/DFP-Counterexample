function run_section6(mode)
% RUN_SECTION6 Generate or reuse the Section 6 data in the working directory.
% SPDX-License-Identifier: MIT
fullRun = ismember(mode, {'paper', 'all'});
drawFigures = fullRun || strcmp(mode, 'figures');
checkCertificate = fullRun || strcmp(mode, 'certify');
if fullRun && isfile('complete.json'), delete('complete.json'); end
if drawFigures && isfile('figures_complete.json'), delete('figures_complete.json'); end
if checkCertificate && isfile('certificate.json'), delete('certificate.json'); end

if fullRun
    geometry = PolarDFP.orbit(100000, 32^3, 1.03);
    assert(geometry.summary.armijo_failures == 0 && geometry.summary.strong_failures == 0, ...
        'The prescribed steps failed the recorded strong Wolfe checks.');
    save('geometry.mat', 'geometry', '-v7.3');
    PolarDFP.writeJSON('geometry.json', geometry.summary);

    objective = PolarDFP.finite(geometry);
    fixedData = [objective.points, objective.corrections, objective.radii];
    save('matched_finite.mat', 'fixedData', '-v7.3');
    PolarDFP.writeJSON('setup.json', struct('j0', geometry.J, 'cycles', geometry.cycles, ...
        'sites', size(objective.points,1), 'x0', geometry.x(1,:), 'H0', geometry.h0, ...
        'max_correction_ratio', objective.ratio, 'terminal_radius', geometry.r(end), ...
        'environment', version, 'threads', maxNumCompThreads));

    for method = {'bfgs', 'dfp'}
        name = method{1};
        run = PolarDFP.solve(objective, geometry, name, 5000);
        record = struct(name, run);
        save([name '.mat'], '-struct', 'record');
        writetable(run.trace, [name '.csv']);
        PolarDFP.writeJSON([name '.json'], run.summary);
        fprintf('%s: %s, %d steps, gradient norm %.12g\n', name, ...
            run.summary.status, run.summary.iterations, run.summary.final_gradient_norm);
    end
end
if drawFigures
    make_polar_figures('paper');
end
if checkCertificate
    python = getenv('PYTHON');
    if isempty(python)
        if ispc, python = 'python'; else, python = 'python3'; end
    end
    checker = fullfile(fileparts(mfilename('fullpath')), 'certify_section6.py');
    args = java.util.ArrayList;
    for arg = {python, checker, 'matched_finite.mat', '--output', 'certificate.json'}
        args.add(java.lang.String(arg{1}));
    end
    builder = java.lang.ProcessBuilder(args);
    builder.directory(java.io.File(pwd));
    builder.inheritIO();
    process = builder.start();
    assert(process.waitFor() == 0, ...
        'Certificate check failed. Set PYTHON to an interpreter with h5py and NumPy.');
end
if fullRun
    PolarDFP.writeJSON('complete.json', struct('status', 'completed', 'matlab', version));
end
end
