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
