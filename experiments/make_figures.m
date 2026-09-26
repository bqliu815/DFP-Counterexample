function make_figures(outDir)
    % MAKE_FIGURES Export Figures 1 and 2 from the stored MATLAB trajectories.
    % SPDX-License-Identifier: MIT
    rawDir = fullfile(outDir, 'raw');
    figDir = fullfile(outDir, 'figures');
    if ~isfolder(figDir)
        mkdir(figDir);
    end
    S = load(fullfile(rawDir, 'geometry.mat'), 'o');
    o = S.o;
    S = load(fullfile(rawDir, 'geometry_bfgs.mat'), 'result');
    b = S.result;
    G = jsondecode(fileread(fullfile(rawDir, 'geometry.json')));
    theta = linspace(0, 2 * pi, 1200)';
    circle = o.centers(end, :) + G.limiting_radius_estimate * [cos(theta), sin(theta)];
    blue = [0, .447, .698];
    red = [.80, .239, .333];
    orange = [.725, .302, 0];
    defaultNames = {'defaultAxesFontName', 'defaultTextFontName', ...
                    'defaultAxesFontSize', 'defaultTextFontSize', 'defaultAxesTickLabelInterpreter', ...
                    'defaultTextInterpreter', 'defaultLegendInterpreter', 'defaultTextColor'};
    oldDefaults = get(groot, defaultNames);
    restoreDefaults = onCleanup(@()set(groot, defaultNames, oldDefaults));
    set(groot, 'defaultAxesFontName', 'Nimbus Sans', 'defaultTextFontName', 'Nimbus Sans', ...
        'defaultAxesFontSize', 8, 'defaultTextFontSize', 8, ...
        'defaultAxesTickLabelInterpreter', 'latex', 'defaultTextInterpreter', 'latex', ...
        'defaultLegendInterpreter', 'latex', 'defaultTextColor', 'k');
    fig = figure('Visible', 'off', 'Color', 'w', 'Units', 'inches', 'Position', [0, 0, 4.75, 2.55], ...
                 'PaperPositionMode', 'auto', 'Renderer', 'painters');
    for j = 1:2
        ax = axes(fig, 'Position', [.105 + (j - 1) * .47, .19, .395, .70]);
        set(ax, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'GridColor', [.5, .5, .5]);
        hold(ax, 'on');
        if j == 1
            h_method = plot(ax, o.x(:, 1), o.x(:, 2), 'Color', blue, 'LineWidth', .32, 'DisplayName', 'DFP');
            h_circle = plot(ax, circle(:, 1), circle(:, 2), '--', ...
                'Color', orange, 'LineWidth', 1, 'DisplayName', 'Circle estimate');
            plot(ax, o.x(1, 1), o.x(1, 2), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 3, 'HandleVisibility', 'off');
            plot(ax, o.x(end, 1), o.x(end, 2), 's', 'Color', blue, ...
                'MarkerFaceColor', 'w', 'MarkerSize', 3.5, 'HandleVisibility', 'off');
            title(ax, '(a) DFP: $100{,}000$ cycles', 'FontWeight', 'normal', 'FontSize', 8);
            ylabel(ax, '$x_2$');
        else
            xx = [b.x0; b.trace.x1, b.trace.x2];
            h_circle = plot(ax, circle(:, 1), circle(:, 2), '--', ...
                'Color', orange, 'LineWidth', 1, 'DisplayName', 'Circle estimate');
            h_method = plot(ax, xx(:, 1), xx(:, 2), '-o', 'Color', red, ...
                'LineWidth', .85, 'MarkerFaceColor', red, 'MarkerSize', 2.3, ...
                'DisplayName', 'BFGS (fminunc)');
            plot(ax, xx(1, 1), xx(1, 2), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 3, 'HandleVisibility', 'off');
            plot(ax, xx(end, 1), xx(end, 2), 'k*', 'MarkerSize', 6, 'HandleVisibility', 'off');
            title(ax, sprintf('(b) BFGS: %d iterations', height(b.trace)), 'FontWeight', 'normal', 'FontSize', 8);
            ax.YTickLabel = {};
        end
        lg = legend(ax, [h_method, h_circle], 'Location', 'northeast', 'FontSize', 7, 'Box', 'on');
        set(lg, 'Color', 'w', 'TextColor', 'k', 'EdgeColor', [.3, .3, .3]);
        lg.ItemTokenSize = [12, 8];
        axis(ax, 'equal');
        xlim(ax, [-1.07, 1.07]);
        ylim(ax, [-1.07, 1.07]);
        xticks(ax, -1:.5:1);
        yticks(ax, -1:.5:1);
        xlabel(ax, '$x_1$');
        ax.XTickLabelRotation = 0;
        grid(ax, 'on');
        box(ax, 'on');
        ax.GridAlpha = .18;
        ax.LineWidth = .5;
    end
    exportgraphics(fig, fullfile(figDir, 'Fig1.pdf'), 'ContentType', 'vector');
    exportgraphics(fig, fullfile(figDir, 'Fig1.png'), 'Resolution', 300);
    close(fig);

    S = load(fullfile(rawDir, 'dfp_0p0025.mat'), 'result');
    d = S.result;
    S = load(fullfile(rawDir, 'bfgs_0p0025.mat'), 'result');
    b = S.result;
    fig = figure('Visible', 'off', 'Color', 'w', 'Units', 'inches', 'Position', [0, 0, 3.6, 2.5], ...
                 'PaperPositionMode', 'auto', 'Renderer', 'painters');
    ax = axes(fig, 'Position', [.19, .25, .71, .69]);
    set(ax, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'GridColor', [.5, .5, .5]);
    hold(ax, 'on');
    plot(ax, d.trace.iteration, d.trace.gradient_norm, 'Color', blue, 'LineWidth', 1, 'DisplayName', 'DFP');
    plot(ax, b.trace.iteration, b.trace.gradient_norm, '--o', 'Color', red, 'LineWidth', 1, ...
         'MarkerSize', 3, 'MarkerIndices', 1:4:height(b.trace), 'MarkerFaceColor', red, 'DisplayName', 'BFGS');
    plot(ax, [1, max([d.summary.iterations, b.summary.iterations])], [1e-10, 1e-10], ':', 'Color', [.3, .3, .3], 'LineWidth', .8, 'DisplayName', 'Tolerance');
    set(ax, 'XScale', 'log', 'YScale', 'log', 'FontSize', 10, 'XLim', [1, 5000], 'YLim', [1e-13, 1.5], ...
        'XTick', [1, 10, 100, 1000, 5000], 'XTickLabel', {'$1$', '$10$', '$100$', '$1{,}000$', '$5{,}000$'}, ...
        'YTick', 10.^(-12:2:0), 'XMinorGrid', 'off', 'YMinorGrid', 'off');
    ax.XTickLabelRotation = 0;
    ax.YTickLabel = compose('$10^{%d}$', -12:2:0);
    xlabel(ax, 'Iteration $k$', 'FontSize', 10);
    ylabel(ax, '$\|g_k\|_2$', 'FontSize', 10, 'Interpreter', 'latex');
    lg = legend(ax, 'Location', 'northeast', 'FontSize', 9, 'Box', 'on');
    set(lg, 'Color', 'w', 'TextColor', 'k', 'EdgeColor', [.3, .3, .3]);
    grid(ax, 'on');
    box(ax, 'on');
    ax.GridAlpha = .18;
    exportgraphics(fig, fullfile(figDir, 'Fig2.pdf'), 'ContentType', 'vector');
    exportgraphics(fig, fullfile(figDir, 'Fig2.png'), 'Resolution', 300);
    close(fig);
    completed = struct('created', true, 'software', version, ...
                       'source', 'Stored MATLAB iterate data', 'smoothing', false, ...
                       'circle_radius', G.limiting_radius_estimate, 'bfgs_iterations', height(b.trace));
    DFPExperiment.writeJSON(fullfile(rawDir, 'figures_complete.json'), completed);
end
