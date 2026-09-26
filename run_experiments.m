function outDir = run_experiments(mode, outDir)
    % RUN_EXPERIMENTS Reproduce the DFP counterexample experiments.
    %   Default: 'smoke'. Use 'all' for the complete paper protocol.
    %   Other modes: 'test', 'probe', 'full', 'certify', 'verify', 'figures'.
    %   RUN_EXPERIMENTS(MODE, OUTDIR) selects the output directory.
    %   See docs/REPRODUCIBILITY.md for stage dependencies and outputs.
    % SPDX-License-Identifier: MIT

    if nargin < 1
        mode = 'smoke';
    end
    if ~(ischar(mode) || (isstring(mode) && isscalar(mode)))
        error('dfp:UnknownMode', 'MODE must be a scalar string or character vector.');
    end
    mode = lower(char(mode));
    modes = {'smoke', 'all', 'test', 'probe', 'full', 'certify', 'verify', 'figures'};
    if ~any(strcmp(mode, modes))
        error('dfp:UnknownMode', 'Unknown mode: %s.', mode);
    end

    repoDir = fileparts(mfilename('fullpath'));
    oldPath = path;
    restorePath = onCleanup(@()path(oldPath));
    addpath(fullfile(repoDir, 'src'), fullfile(repoDir, 'experiments'));

    if nargin < 2 || isempty(outDir)
        if any(strcmp(mode, {'smoke', 'test', 'probe'}))
            outDir = fullfile(repoDir, 'results', mode);
        else
            outDir = fullfile(repoDir, 'results', 'paper');
        end
    end
    outDir = char(outDir);
    if ~isfolder(outDir)
        [ok, message] = mkdir(outDir);
        if ~ok
            error('dfp:OutputDirectory', '%s', message);
        end
    end
    [ok, attributes] = fileattrib(outDir);
    assert(ok, 'dfp:OutputDirectory', 'Cannot resolve the output directory.');
    outDir = attributes.Name;
    rawDir = fullfile(outDir, 'raw');
    if ~isfolder(rawDir)
        mkdir(rawDir);
    end

    if strcmp(mode, 'test')
        clearMarkers(rawDir, {'tests_complete.json'});
    elseif strcmp(mode, 'smoke')
        clearMarkers(rawDir, {'smoke.json'});
    end

    if any(strcmp(mode, {'smoke', 'full', 'all', 'test'}))
        requireToolbox('Optimization_Toolbox', 'fminunc', 'Optimization Toolbox');
        requireToolbox('Statistics_Toolbox', 'KDTreeSearcher', ...
                       'Statistics and Machine Learning Toolbox');
    end
    if any(strcmp(mode, {'certify', 'all'}))
        requireToolbox('Symbolic_Toolbox', 'sym', 'Symbolic Math Toolbox');
        % Start the symbolic engine before the long numerical stage.
        sym(0);
    end
    if any(strcmp(mode, {'figures', 'all'})) && ~usejava('jvm')
        error('dfp:GraphicsNeedsJVM', ...
              'Figure export requires MATLAB with the JVM; omit -nojvm.');
    end

    record_environment(outDir, mode);
    started = tic;
    if strcmp(mode, 'probe')
        fprintf('Environment recorded in %s\n', outDir);
        return
    elseif strcmp(mode, 'test')
        results = runtests(fullfile(repoDir, 'tests'));
        assertSuccess(results);
        DFPExperiment.writeJSON(fullfile(rawDir, 'tests_complete.json'), ...
                                struct('passed', true, 'tests', numel(results)));
        return
    elseif strcmp(mode, 'all')
        stages = {'full', 'certify', 'verify', 'figures'};
    else
        stages = {mode};
    end

    for k = 1:numel(stages)
        stage = stages{k};
        switch stage
            case 'full'
                required = {};
                stale = {'full_complete.json', 'certificates_complete.json', ...
                         'verification.json', 'figures_complete.json'};
            case 'certify'
                required = {'full_complete.json'};
                stale = {'certificates_complete.json', 'verification.json'};
            case 'verify'
                required = {'full_complete.json', 'certificates_complete.json'};
                stale = {'verification.json'};
            case 'figures'
                required = {'full_complete.json', 'geometry.mat', 'geometry.json', ...
                            'geometry_bfgs.mat', 'dfp_0p0025.mat', 'bfgs_0p0025.mat'};
                stale = {'figures_complete.json'};
            otherwise
                required = {};
                stale = {};
        end
        % A failed run must not leave an earlier success record.
        clearMarkers(rawDir, stale);
        requireFiles(rawDir, required);
        run_suite(stage, outDir);
    end
    fprintf('Completed %s in %.2f seconds. Results: %s\n', mode, toc(started), outDir);
end

function requireToolbox(feature, entryPoint, displayName)
    if ~license('test', feature) || isempty(which(entryPoint))
        error('dfp:MissingToolbox', '%s is required for this mode.', displayName);
    end
end

function requireFiles(folder, names)
    for k = 1:numel(names)
        if ~isfile(fullfile(folder, names{k}))
            error('dfp:MissingInput', ...
                  'Missing %s. Run the preceding reproduction stage first.', names{k});
        end
    end
end

function clearMarkers(folder, names)
    for k = 1:numel(names)
        file = fullfile(folder, names{k});
        if isfile(file)
            delete(file);
        end
    end
end
