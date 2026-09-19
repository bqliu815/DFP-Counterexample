function record_environment(outDir, mode)
    % RECORD_ENVIRONMENT Record MATLAB version, toolboxes, and execution settings.
    % SPDX-License-Identifier: MIT

    env.timestamp_utc = char(datetime('now', 'TimeZone', 'UTC', ...
                                      'Format', 'yyyy-MM-dd''T''HH:mm:ss''Z'''));
    env.mode = mode;
    env.matlab_version = version;
    env.matlab_release = version('-release');
    env.architecture = computer;
    env.toolboxes = ver;
    env.precision = 'IEEE 754 binary64';
    env.computational_threads = maxNumCompThreads;
    env.jvm = usejava('jvm');
    DFPExperiment.writeJSON(fullfile(outDir, 'raw', ...
                                     ['environment_', mode, '.json']), env);
end
