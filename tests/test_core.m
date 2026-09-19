function tests = test_core
    % TEST_CORE Check updates, line searches, interpolation, and the run entry point.
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
    checkUpdate(testCase, 'dfp');
end

function testBFGSUpdate(testCase)
    checkUpdate(testCase, 'bfgs');
end

function checkUpdate(testCase, method)
    H = [2, .2; .2, 1];
    s = [.3; -.2];
    y = [2, .1; .1, 3] * s;
    updated = DFPExperiment.update(H, s, y, method);
    verifyEqual(testCase, updated, updated', 'AbsTol', 1e-15);
    verifyLessThan(testCase, norm(updated * y - s), 1e-14);
    verifyGreaterThan(testCase, min(eig(updated)), 0);
end

function testUnknownUpdateIsRejected(testCase)
    verifyError(testCase, @()DFPExperiment.update(eye(2), [1; 0], ...
                                                  [1; 0], 'unknown'), 'dfp:UnknownMethod');
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

function testWeakWolfeRejectedFirstTrial(testCase)
    checkSearch(testCase, 'weak_unit', NaN);
end

function testCubicZoomRejectedFirstTrial(testCase)
    checkSearch(testCase, 'zoom_unit', NaN);
end

function testMoreThuenteRejectedFirstTrial(testCase)
    checkSearch(testCase, 'mt_unit', NaN);
end

function testHistoryInitialization(testCase)
    for policy = {'zoom_history', 'mt_history'}
        checkSearch(testCase, policy{1}, 42);
    end
end

function testUnknownSearchIsRejected(testCase)
    evaluate = @(x)quadratic(x, eye(2));
    x = [1; 2];
    [f, g] = evaluate(x);
    for policy = {'weak_typo', 'zoom_typo', 'mt_typo'}
        verifyError(testCase, @()wolfe_search(evaluate, x, f, g, -g, policy{1}, NaN), ...
                    'dfp:UnknownPolicy');
    end
end

function checkSearch(testCase, policy, previous)
    A = diag([2, 20]);
    evaluate = @(x)quadratic(x, A);
    x = [1; 2];
    [f, g] = evaluate(x);
    d = -g;
    [alpha, fp, gp, trials] = wolfe_search(evaluate, x, f, g, d, policy, previous);
    verifyLessThanOrEqual(testCase, fp, f + .25 * alpha * (g' * d));
    if startsWith(policy, 'weak')
        verifyGreaterThanOrEqual(testCase, gp' * d, .75 * (g' * d));
    else
        verifyLessThanOrEqual(testCase, abs(gp' * d), .75 * abs(g' * d));
    end
    if contains(policy, 'history')
        expected = min(1, 1.01 * 2 * (f - previous) / (g' * d));
        verifyEqual(testCase, trials(1, 1), expected, 'AbsTol', 1e-15);
    else
        verifyEqual(testCase, trials(1, 1), 1);
        verifyGreaterThan(testCase, size(trials, 1), 1);
    end
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

function [f, g] = quadratic(x, A)
    f = .5 * (x' * A * x);
    g = A * x;
end
