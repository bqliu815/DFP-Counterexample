function result = run_fminunc(obj, orbit, method, maxIterations, initialization)
% RUN_FMINUNC Minimize a finite interpolant with MATLAB's DFP or BFGS.
%   INITIALIZATION is 'prescribed' (default) or 'identity'. For the former,
%   x = x0 + L*z, L*L' = H0, gives the prescribed initial search direction.
%   MATLAB's internal line search, first-update scaling, and safeguards are
%   left unchanged. Gradients and the stopping test use the original x units.
% SPDX-License-Identifier: MIT

if nargin < 4, maxIterations = 5000; end
if nargin < 5, initialization = 'prescribed'; end
method = validatestring(method, {'dfp', 'bfgs'});
initialization = validatestring(initialization, {'prescribed', 'identity'});
validateattributes(maxIterations, {'numeric'}, {'scalar', 'integer', 'positive'});
x0 = orbit.x(1, :)';
L = eye(2);
if strcmp(initialization, 'prescribed')
    L = chol(orbit.h0, 'lower');
end
tolerance = 1e-10;
maxEvaluations = 100000;
rows = zeros(maxIterations + 1, 12);
nrows = 0;
evaluations = 0;
monitorEvaluations = 0;
previousX = x0;
[initialF, previousG] = DFPExperiment.valueGrad(obj, x0);
initialG = previousG;
previousF = initialF;
options = optimset('Algorithm', 'quasi-newton', 'HessUpdate', method, ...
    'GradObj', 'on', 'Display', 'off', 'TolFun', 0, 'TolX', 0, ...
    'MaxIter', maxIterations, 'MaxFunEvals', maxEvaluations, 'OutputFcn', @observe);
started = tic;
[z, fval, exitflag, output] = fminunc(@evaluate, zeros(2, 1), options);
elapsed = toc(started);
x = x0 + L*z;
[checkedF, g] = DFPExperiment.valueGrad(obj, x);
assert(abs(fval - checkedF) <= 1e-13 * max(1, abs(checkedF)));
status = 'solver_stopped';
if norm(g) <= tolerance
    status = 'gradient_tolerance';
elseif output.iterations >= maxIterations
    status = 'iteration_limit';
elseif output.funcCount >= maxEvaluations
    status = 'evaluation_limit';
end
names = {'iteration', 'function_value', 'gradient_norm', 'x1', 'x2', ...
         'g1', 'g2', 'alpha', 'function_evaluations', 'step_norm', ...
         'armijo_ratio', 'strong_ratio'};
result.trace = array2table(rows(1:nrows, :), 'VariableNames', names);
result.x0 = x0';
result.x = x;
result.h0 = L*L';
result.transform = L;
result.output = output;
result.exitflag = exitflag;
result.options = rmfield(options, 'OutputFcn');
result.summary = struct('solver', 'fminunc', 'epsilon0', obj.epsilon0, ...
    'method', method, 'initialization', initialization, 'status', status, ...
    'exitflag', exitflag, 'iterations', output.iterations, 'recorded_steps', nrows, ...
    'function_evaluations', output.funcCount, 'counted_evaluations', evaluations, ...
    'monitor_evaluations', monitorEvaluations, 'elapsed_seconds', elapsed, ...
    'prefix_steps', obj.steps, 'hessian_band', obj.band, ...
    'initial_gradient_norm', norm(initialG), 'final_gradient_norm', norm(g), ...
    'final_function_value', checkedF, 'gradient_tolerance', tolerance, ...
    'max_iterations', maxIterations, 'max_function_evaluations', maxEvaluations, ...
    'matlab_release', version('-release'));

    function [f, gz] = evaluate(z)
        [f, gx] = DFPExperiment.valueGrad(obj, x0 + L*z);
        gz = L'*gx;
        evaluations = evaluations + 1;
    end

    function stop = observe(z, values, state)
        currentX = x0 + L*z;
        [currentF, currentG] = DFPExperiment.valueGrad(obj, currentX);
        monitorEvaluations = monitorEvaluations + 1;
        stop = norm(currentG) <= tolerance || values.iteration >= maxIterations;
        if ~strcmp(state, 'iter') || values.iteration == 0
            return
        end
        step = currentX - previousX;
        q = -previousG'*step;
        armijo = NaN;
        strong = NaN;
        if q > 0
            armijo = (previousF - currentF)/q;
            strong = abs(currentG'*step)/q;
        end
        nrows = nrows + 1;
        rows(nrows, :) = [values.iteration, currentF, norm(currentG), ...
            currentX', currentG', values.lssteplength, values.funccount, ...
            norm(step), armijo, strong];
        previousX = currentX;
        previousF = currentF;
        previousG = currentG;
    end
end
