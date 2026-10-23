function tests = test_core
    % TEST_CORE Check the polar recurrence, interpolation, solver, and entry point.
    % SPDX-License-Identifier: MIT
    tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    root = fileparts(fileparts(mfilename('fullpath')));
    testCase.TestData.oldPath = path;
    testCase.TestData.root = root;
    addpath(root, fullfile(root, 'polar'));
    testCase.TestData.orbit = PolarDFP.orbit(200, 400^3, 1.03);
    testCase.TestData.objective = PolarDFP.finite(testCase.TestData.orbit);
end

function teardownOnce(testCase)
    path(testCase.TestData.oldPath);
end

function testPrescribedRecurrence(testCase)
    orbit = testCase.TestData.orbit;
    summary = orbit.summary;
    verifySize(testCase, orbit.x, [401, 2]);
    verifyEqual(testCase, summary.armijo_failures, 0);
    verifyEqual(testCase, summary.strong_failures, 0);
    verifyLessThan(testCase, summary.max_secant_residual, 1e-10);
    verifyLessThan(testCase, summary.max_update_residual, 1e-10);
    verifyGreaterThan(testCase, summary.min_normalized_eigenvalue, 0);
    verifyGreaterThan(testCase, summary.min_sigma, 0);
    verifyGreaterThan(testCase, min(eig(orbit.h0)), 0);
    verifyEqual(testCase, orbit.h0, orbit.h0', 'AbsTol', 1e-15);
    verifyEqual(testCase, orbit.r(end), 1, 'AbsTol', 1e-14);
    verifyTrue(testCase, all(diff(orbit.r) < 0));

    step = (orbit.x(2, :) - orbit.x(1, :))';
    direction = -orbit.h0 * orbit.g(1, :)';
    verifyLessThan(testCase, norm(step / norm(step) - direction / norm(direction)), 1e-8);
end

function testInterpolationAtNodes(testCase)
    orbit = testCase.TestData.orbit;
    objective = testCase.TestData.objective;
    verifyGreaterThan(testCase, objective.band(1), 0);
    for k = 1:size(objective.points, 1)
        x = objective.points(k, :)';
        [f, g] = PolarDFP.valueGrad(objective, x);
        verifyEqual(testCase, f, .5 * (x' * x), 'AbsTol', 1e-14);
        verifyEqual(testCase, g, orbit.g(k, :)', 'AbsTol', 1e-13);
    end
end

function testNewFeedbackAndGeneralWolfe(testCase)
    for pair=[.2,.8;.3,.9]
        pars=PolarDFP.wolfeParameters(pair(1),pair(2),.5);
        orbit=PolarDFP.orbit(40,400^3,pars);
        expected=pars.p*pars.gamma./(orbit.A(1:end-1).*orbit.nu*(1/pars.a+1/pars.b));
        verifyEqual(testCase,orbit.phi,expected,'RelTol',1e-14);
        verifyEqual(testCase,orbit.summary.armijo_failures,0);
        verifyEqual(testCase,orbit.summary.strong_failures,0);
        verifyLessThan(testCase,orbit.summary.max_update_residual,1e-10);
        inverse=(pars.b-pars.a+2*orbit.nu)./(orbit.nu.*(pars.a-orbit.nu).*(pars.b+orbit.nu));
        verifyEqual(testCase,inverse,orbit.A(1:end-1),'RelTol',2e-13);
    end
end

function testBetaRecurrenceAndScalingRatios(testCase)
    % Three exponents share a single shape, Wolfe pair, and cycle index range.
    % Direct coordinate differences independently check the stable diagnostics.
    N=80;J=32^3;base=2*pi/9;indices=J+(0:N-1)';
    for beta=[1/3,1/2,2/3]
        pars=PolarDFP.parameters(2*base,base,1.03*base,beta,.25,.75);
        orbit=PolarDFP.orbit(N,J,pars);
        verifyEqual(testCase,orbit.parameters.beta,beta);
        verifyEqual(testCase,orbit.parameters.p,1/(1+beta),'RelTol',2*eps);
        verifyEqual(testCase,orbit.summary.nu_over_b,1.03,'RelTol',2*eps);
        verifyEqual(testCase,orbit.tau,(J+(0:N)').^(-pars.p),'RelTol',2*eps);
        expectedAngles=((indices+1).^(1-pars.p)-indices.^(1-pars.p))/(1-pars.p);
        verifyEqual(testCase,orbit.delta(1:2:end)/pars.a,expectedAngles,'RelTol',1e-10);
        verifyEqual(testCase,orbit.delta(2:2:end)/pars.b,expectedAngles,'RelTol',1e-10);
        verifyEqual(testCase,orbit.summary.armijo_failures,0);
        verifyEqual(testCase,orbit.summary.strong_failures,0);
        verifyLessThan(testCase,orbit.summary.max_secant_residual,1e-10);
        verifyLessThan(testCase,orbit.summary.max_update_residual,1e-10);
        verifySize(testCase,orbit.ratio_radial,[N,1]);
        verifySize(testCase,orbit.ratio_correction,[N,1]);
        verifyTrue(testCase,all(isfinite(orbit.ratio_radial)));
        verifyTrue(testCase,all(isfinite(orbit.ratio_correction)));
        verifyGreaterThan(testCase,min(orbit.ratio_radial),0);
        verifyGreaterThan(testCase,min(orbit.ratio_correction),0);

        denominator=(orbit.targets.D0+orbit.targets.D1)*orbit.tau(1:end-1).^2;
        directRadial=log(orbit.r(1:2:end-2)./orbit.r(3:2:end))./denominator;
        directCorrection=indices.*(1-orbit.zeta(2:2:end))/orbit.targets.L;
        % The direct forms lose O(eps/tau^2) and O(n_j*eps/L), respectively.
        verifyEqual(testCase,orbit.ratio_radial,directRadial, ...
            'AbsTol',32*eps/min(denominator));
        verifyEqual(testCase,orbit.ratio_correction,directCorrection, ...
            'AbsTol',4*eps*max(indices)/orbit.targets.L);
        verifyEqual(testCase,orbit.targets.L,pars.p/(pars.nu*(1/pars.a+1/pars.b)), ...
            'RelTol',2*eps);
    end
end

function testStreamedScalingMatchesFullOrbit(testCase)
    % The storage reduction must preserve the recurrence and exact extrema.
    N=1000;J=32^3;base=2*pi/9;
    for beta=[1/3,1/2,2/3]
        pars=PolarDFP.parameters(2*base,base,1.03*base,beta,.25,.75);
        full=PolarDFP.orbit(N,J,pars);compact=scaling_diagnostics(N,J,pars);
        indices=compact.cycle+1;
        verifyEqual(testCase,compact.ratio_radial,full.ratio_radial(indices),'AbsTol',1e-13);
        verifyEqual(testCase,compact.ratio_correction,full.ratio_correction(indices),'AbsTol',1e-13);
        verifyEqual(testCase,compact.gamma,full.A(indices),'AbsTol',1e-13);
        verifyEqual(testCase,compact.phi,full.phi(indices),'AbsTol',1e-13);
        [value,index]=max(full.ratio_radial);
        verifyEqual(testCase,compact.summary.max_ratio_radial,value,'AbsTol',1e-13);
        verifyEqual(testCase,compact.summary.peak_radial_cycle,index-1);
        [value,index]=max(full.ratio_correction);
        verifyEqual(testCase,compact.summary.max_ratio_correction,value,'AbsTol',1e-13);
        verifyEqual(testCase,compact.summary.peak_correction_cycle,index-1);
        verifyEqual(testCase,compact.summary.final_gamma,full.A(end),'AbsTol',1e-13);
        verifyEqual(testCase,compact.summary.armijo_failures,full.summary.armijo_failures);
        verifyEqual(testCase,compact.summary.strong_failures,full.summary.strong_failures);
        verifyEqual(testCase,compact.summary.min_armijo,full.summary.min_armijo,'AbsTol',1e-13);
        verifyEqual(testCase,compact.summary.max_strong,full.summary.max_strong,'AbsTol',1e-13);
        verifyLessThan(testCase,compact.summary.sampled_max_secant_residual,1e-10);
        verifyLessThan(testCase,compact.summary.sampled_max_update_residual,1e-10);
        verifyEqual(testCase,compact.cycle([1,end]),[0;N-1]);
    end
end

function testScalingAcrossBlockBoundary(testCase)
    % Carrying the matrix into the next scale block must preserve the path.
    N=100005;J=32^3;base=2*pi/9;
    pars=PolarDFP.parameters(2*base,base,1.03*base,.5,.25,.75);
    full=PolarDFP.orbit(N,J,pars);compact=scaling_diagnostics(N,J,pars);
    indices=compact.cycle+1;
    verifyEqual(testCase,compact.ratio_radial,full.ratio_radial(indices),'AbsTol',1e-13);
    verifyEqual(testCase,compact.ratio_correction,full.ratio_correction(indices),'AbsTol',1e-13);
    verifyEqual(testCase,compact.gamma,full.A(indices),'AbsTol',1e-13);
    verifyEqual(testCase,compact.phi,full.phi(indices),'AbsTol',1e-13);
    verifyEqual(testCase,compact.summary.final_gamma,full.A(end),'AbsTol',1e-13);
    verifyEqual(testCase,compact.summary.final_ratio_radial,full.ratio_radial(end),'AbsTol',1e-13);
    verifyEqual(testCase,compact.summary.final_ratio_correction,full.ratio_correction(end),'AbsTol',1e-13);
end

function testGradientInsideTransition(testCase)
    % The larger scale keeps the cutoff correction above subtraction noise.
    orbit = PolarDFP.orbit(20, 32^3, 1.03);
    objective = PolarDFP.finite(orbit);
    k = 10;
    verifyGreaterThan(testCase, norm(objective.corrections(k, :)), 1e-8);
    x = objective.points(k, :)' + .6 * objective.radii(k) * [.6; .8];
    [~, gradient] = PolarDFP.valueGrad(objective, x);
    h = 1e-3 * objective.radii(k);
    finiteDifference = zeros(2, 1);
    for j = 1:2
        delta = zeros(2, 1);
        delta(j) = h;
        finiteDifference(j) = (PolarDFP.valueGrad(objective, x + delta) ...
                               - PolarDFP.valueGrad(objective, x - delta)) / (2 * h);
    end
    verifyLessThan(testCase, norm(finiteDifference - gradient), 1e-7);
end

function testQuadraticOutsideSupports(testCase)
    x = [10; -10];
    [f, g] = PolarDFP.valueGrad(testCase.TestData.objective, x);
    verifyEqual(testCase, g, x);
    verifyEqual(testCase, f, .5 * (x' * x), 'AbsTol', 1e-12);
end

function testBFGSFinalGradientAndTrace(testCase)
    objective = testCase.TestData.objective;
    orbit = testCase.TestData.orbit;
    result = PolarDFP.solve(objective, orbit, 'bfgs', 100);
    [~, finalGradient] = PolarDFP.valueGrad(objective, result.x);
    verifyEqual(testCase, result.summary.status, 'gradient_tolerance');
    verifyLessThanOrEqual(testCase, norm(finalGradient), 1e-10);
    verifyEqual(testCase, result.summary.final_gradient_norm, norm(finalGradient), 'AbsTol', 1e-14);
    verifyGreaterThan(testCase, height(result.trace), 1);
    verifyEqual(testCase, sum(result.trace.iteration == 0), 1);
    verifyTrue(testCase, all(diff(result.trace.iteration) > 0));
    verifyEqual(testCase, [result.trace.x1(end); result.trace.x2(end)], result.x, 'AbsTol', 1e-14);
    verifyEqual(testCase, [result.trace.x1(1), result.trace.x2(1)], orbit.x(1, :), 'AbsTol', 1e-14);
    for k = 1:height(result.trace)
        x = [result.trace.x1(k); result.trace.x2(k)];
        [f, g] = PolarDFP.valueGrad(objective, x);
        verifyEqual(testCase, result.trace.function_value(k), f, 'AbsTol', 1e-14);
        verifyEqual(testCase, result.trace.gradient_norm(k), norm(g), 'AbsTol', 1e-14);
    end
end

function testDFPShortRun(testCase)
    objective = testCase.TestData.objective;
    orbit = testCase.TestData.orbit;
    result = PolarDFP.solve(objective, orbit, 'dfp', 3);
    [f, g] = PolarDFP.valueGrad(objective, result.x);
    verifyEqual(testCase, result.summary.method, 'dfp');
    verifyEqual(testCase, result.summary.status, 'iteration_limit');
    verifyEqual(testCase, result.summary.iterations, 3);
    verifyEqual(testCase, result.trace.iteration, (0:3)');
    verifyEqual(testCase, [result.trace.x1(end); result.trace.x2(end)], result.x, 'AbsTol', 1e-14);
    verifyEqual(testCase, result.trace.function_value(end), f, 'AbsTol', 1e-14);
    verifyEqual(testCase, result.trace.gradient_norm(end), norm(g), 'AbsTol', 1e-14);
    verifyEqual(testCase, result.summary.final_gradient_norm, norm(g), 'AbsTol', 1e-14);
end

function testFailedCertificateRemovesOldResult(testCase)
    folder = tempname(tempdir);
    mkdir(folder);
    cleanup = onCleanup(@()rmdir(folder, 's'));
    stem = fullfile(folder, 'missing');
    certificatePath = [stem, '_certificate.json'];
    PolarDFP.writeJSON(certificatePath, struct('status', 'exact_dyadic_checks_passed'));
    verifyError(testCase, @()certify_finite(stem), 'MATLAB:load:couldNotReadFile');
    verifyFalse(testCase, isfile(certificatePath));
end

function testEntryPointPreservesSession(testCase)
    folder = tempname(tempdir);
    mkdir(folder);
    oldDirectory = pwd;
    oldPath = path;
    cleanup = onCleanup(@()restoreSession(oldDirectory, oldPath, folder));
    rmpath(fullfile(testCase.TestData.root, 'polar'));
    cd(folder);
    expectedDirectory = pwd;
    expectedPath = path;
    outputDir = fullfile(folder, 'output');
    run_experiments('smoke', outputDir);
    verifyEqual(testCase, pwd, expectedDirectory);
    verifyEqual(testCase, path, expectedPath);
    smoke = jsondecode(fileread(fullfile(outputDir, 'smoke.json')));
    verifyEqual(testCase, smoke.recurrence.armijo_failures, 0);
    verifyEqual(testCase, smoke.recurrence.strong_failures, 0);
    verifyLessThanOrEqual(testCase, smoke.solver.final_gradient_norm, 1e-10);
end

function restoreSession(oldDirectory, oldPath, folder)
    cd(oldDirectory);
    path(oldPath);
    if isfolder(folder)
        rmdir(folder, 's');
    end
end
