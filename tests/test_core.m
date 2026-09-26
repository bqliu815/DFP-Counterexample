function tests = test_core
    % TEST_CORE Check the recurrence, interpolation, fminunc adapter, and entry point.
    % SPDX-License-Identifier: MIT
    tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    root = fileparts(fileparts(mfilename('fullpath')));
    testCase.TestData.oldPath = path;
    addpath(root, fullfile(root, 'src'), fullfile(root, 'experiments'));
end

function teardownOnce(testCase)
    path(testCase.TestData.oldPath);
end

function testDFPUpdate(testCase)
    H = [2, .2; .2, 1];
    s = [.3; -.2];
    y = [2, .1; .1, 3] * s;
    updated = DFPExperiment.update(H, s, y);
    verifyEqual(testCase, updated, updated', 'AbsTol', 1e-15);
    verifyLessThan(testCase, norm(updated * y - s), 1e-14);
    verifyGreaterThan(testCase, min(eig(updated)), 0);
end

function testPrescribedRecurrence(testCase)
    orbit = DFPExperiment.oracle(100, .03);
    summary = DFPExperiment.geometrySummary(orbit);
    verifySize(testCase, orbit.x, [201, 2]);
    verifyEqual(testCase, summary.armijo_failures, 0);
    verifyEqual(testCase, summary.strong_failures, 0);
    verifyLessThan(testCase, summary.max_secant_residual, 1e-12);
    verifyGreaterThan(testCase, summary.min_metric_eigenvalue, 0);
end

function testFminuncOriginalCoordinates(testCase)
    [obj, orbit] = DFPExperiment.finite(.03);
    orbit.x(1, :) = [10, -10];
    orbit.h0 = [2, .25; .25, 1];
    for method = {'dfp', 'bfgs'}
        r = run_fminunc(obj, orbit, method{1}, 100);
        verifyEqual(testCase, r.summary.status, 'gradient_tolerance');
        verifyLessThanOrEqual(testCase, r.summary.final_gradient_norm, 1e-10);
        verifyEqual(testCase, r.h0, orbit.h0, 'RelTol', 1e-14);
        verifyEqual(testCase, r.summary.counted_evaluations, r.output.funcCount);
    end
end

function testPrescribedInitialDirectionAndBudget(testCase)
    [obj, orbit] = DFPExperiment.finite(.03);
    r = run_fminunc(obj, orbit, 'dfp', 3);
    verifyEqual(testCase, r.summary.status, 'iteration_limit');
    verifyEqual(testCase, r.summary.iterations, 3);
    verifyEqual(testCase, r.h0, orbit.h0, 'RelTol', 1e-14);
    [~, g0] = DFPExperiment.valueGrad(obj, orbit.x(1, :)');
    step = [r.trace.x1(1); r.trace.x2(1)] - orbit.x(1, :)';
    expected = -r.trace.alpha(1)*orbit.h0*g0;
    verifyLessThan(testCase, norm(step - expected), 1e-15);
end

function testFminuncTraceAndStoppingTest(testCase)
    [obj, orbit] = DFPExperiment.finite(.03);
    r = run_fminunc(obj, orbit, 'bfgs', 200);
    for k = 1:height(r.trace)
        [f, g] = DFPExperiment.valueGrad(obj, [r.trace.x1(k); r.trace.x2(k)]);
        verifyEqual(testCase, f, r.trace.function_value(k), 'AbsTol', 1e-14);
        verifyEqual(testCase, norm(g), r.trace.gradient_norm(k), 'AbsTol', 1e-14);
    end
    [~, g] = DFPExperiment.valueGrad(obj, r.x);
    verifyEqual(testCase, strcmp(r.summary.status, 'gradient_tolerance'), norm(g) <= 1e-10);
    verifyEqual(testCase, r.summary.counted_evaluations, r.output.funcCount);
end

function testInterpolationAtNodes(testCase)
    [obj, orbit] = DFPExperiment.finite(.03);
    for k = [1, 10, size(obj.points, 1)]
        x = obj.points(k, :)';
        [f, g] = DFPExperiment.valueGrad(obj, x);
        verifyEqual(testCase, g, orbit.g(k, :)', 'AbsTol', 1e-13);
        verifyEqual(testCase, f, .5 * norm(x - obj.center_ref)^2, 'AbsTol', 1e-14);
    end
end

function testGradientInsideTransition(testCase)
    [obj, ~] = DFPExperiment.finite(.03);
    k = 10;
    x = obj.points(k, :)' + .6 * obj.radii(k) * [.6; .8];
    [~, gradient] = DFPExperiment.valueGrad(obj, x);
    h = 1e-4 * obj.radii(k);
    finiteDifference = zeros(2, 1);
    for j = 1:2
        delta = zeros(2, 1);
        delta(j) = h;
        finiteDifference(j) = (DFPExperiment.valueGrad(obj, x + delta) ...
                               - DFPExperiment.valueGrad(obj, x - delta)) / (2 * h);
    end
    verifyLessThan(testCase, norm(finiteDifference - gradient), 1e-5);
end

function testQuadraticOutsideSupports(testCase)
    [obj, ~] = DFPExperiment.finite(.03);
    x = [10; -10];
    [f, g] = DFPExperiment.valueGrad(obj, x);
    verifyEqual(testCase, g, x - obj.center_ref);
    verifyEqual(testCase, f, .5 * norm(x - obj.center_ref)^2, 'AbsTol', 1e-12);
end

function testEntryPointPreservesSession(testCase)
    folder = tempname;
    mkdir(folder);
    cleanup = onCleanup(@()rmdir(folder, 's'));
    oldDirectory = pwd;
    oldPath = path;
    run_experiments('probe', folder);
    verifyEqual(testCase, pwd, oldDirectory);
    verifyEqual(testCase, path, oldPath);
    env = jsondecode(fileread(fullfile(folder, 'raw', 'environment_probe.json')));
    verifyFalse(testCase, isfield(env, 'hostname'));
    verifyEqual(testCase, env.precision, 'IEEE 754 binary64');
end

function testMissingInputsFailClearly(testCase)
    folder = tempname;
    mkdir(folder);
    cleanup = onCleanup(@()rmdir(folder, 's'));
    mkdir(fullfile(folder, 'raw'));
    DFPExperiment.writeJSON(fullfile(folder, 'raw', 'verification.json'), struct('passed', true));
    verifyError(testCase, @()run_experiments('verify', folder), 'dfp:MissingInput');
    verifyFalse(testCase, isfile(fullfile(folder, 'raw', 'verification.json')));
end
