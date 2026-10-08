function make_polar_figures
    % MAKE_POLAR_FIGURES Export the stored polar experiments in the public figure style.
    % SPDX-License-Identifier: MIT
    figDir = pwd;
    S = load('geometry.mat', 'geometry');
    o = S.geometry;
    S = load('bfgs.mat', 'bfgs');
    b = S.bfgs;
    theta = linspace(0, 2 * pi, 1200)';
    circle = [cos(theta), sin(theta)];
    axisLimit = max(1.2, ceil(10 * max(vecnorm(o.x, 2, 2))) / 10);
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
            h_method = plot(ax, o.x(:, 1), o.x(:, 2), 'Color', blue, 'LineWidth', .16, 'DisplayName', 'DFP');
            h_circle = plot(ax, circle(:, 1), circle(:, 2), '--', ...
                'Color', orange, 'LineWidth', 1, 'DisplayName', 'Unit circle');
            plot(ax, o.x(1, 1), o.x(1, 2), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 3, 'HandleVisibility', 'off');
            plot(ax, o.x(end, 1), o.x(end, 2), 's', 'Color', blue, ...
                'MarkerFaceColor', 'w', 'MarkerSize', 3.5, 'HandleVisibility', 'off');
            title(ax, '(a) DFP: $100{,}000$ cycles', 'FontWeight', 'normal', 'FontSize', 8);
        else
            xx = [b.trace.x1, b.trace.x2];
            h_circle = plot(ax, circle(:, 1), circle(:, 2), '--', ...
                'Color', orange, 'LineWidth', 1, 'DisplayName', 'Unit circle');
            h_method = plot(ax, xx(:, 1), xx(:, 2), '-o', 'Color', red, ...
                'LineWidth', .85, 'MarkerFaceColor', red, 'MarkerSize', 2.3, ...
                'DisplayName', 'BFGS');
            plot(ax, xx(1, 1), xx(1, 2), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 3, 'HandleVisibility', 'off');
            plot(ax, xx(end, 1), xx(end, 2), 'k*', 'MarkerSize', 6, 'HandleVisibility', 'off');
            title(ax, sprintf('(b) BFGS: %d iterations', b.summary.iterations), 'FontWeight', 'normal', 'FontSize', 8);
        end
        lg = legend(ax, [h_method, h_circle], 'Location', 'northeast', 'FontSize', 9, 'Box', 'on');
        set(lg, 'Color', 'w', 'TextColor', 'k', 'EdgeColor', [.3, .3, .3]);
        lg.ItemTokenSize = [12, 8];
        axis(ax, 'equal');
        xlim(ax, [-axisLimit, axisLimit]);
        ylim(ax, [-axisLimit, axisLimit]);
        xticks(ax, -1:.5:1);
        yticks(ax, -1:.5:1);
        xlabel(ax, '$x_1$');
        ylabel(ax, '$x_2$');
        ax.XTickLabelRotation = 0;
        grid(ax, 'on');
        box(ax, 'on');
        ax.GridAlpha = .18;
        ax.LineWidth = .5;
        if j == 1
            rectangle(ax, 'Position', [-.07, -1.12, .14, .14], ...
                'EdgeColor', [.25, .25, .25], 'LineWidth', .5);
        end
    end
    detail = axes(fig, 'Position', [.235, .40, .16, .24]);
    plot(detail, o.x(:, 1), o.x(:, 2), 'Color', blue, 'LineWidth', .4);
    axis(detail, 'equal');
    xlim(detail, [-.07, .07]); ylim(detail, [-1.12, -.98]);
    set(detail, 'XTick', [-.05, .05], 'YTick', [-1.10, -1.05, -1], ...
        'FontSize', 8, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', ...
        'LineWidth', .5, 'XTickLabelRotation', 0);
    box(detail, 'on');
    title(detail, 'Detail', 'FontSize', 8, 'FontWeight', 'normal');
    annotation(fig, 'arrow', [.3025, .3025], [.25417, .397], ...
        'Color', [.25, .25, .25], 'LineWidth', .55, ...
        'HeadLength', 4, 'HeadWidth', 4);
    exportgraphics(fig, fullfile(figDir, 'Fig1.pdf'), 'ContentType', 'vector');
    exportgraphics(fig, fullfile(figDir, 'Fig1.png'), 'Resolution', 300);
    close(fig);

    S = load('dfp.mat', 'dfp');
    d = S.dfp;
    S = load('bfgs.mat', 'bfgs');
    b = S.bfgs;
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
    lg = legend(ax, 'Location', 'northeast', 'FontSize', 10, 'Box', 'on');
    set(lg, 'Color', 'w', 'TextColor', 'k', 'EdgeColor', [.3, .3, .3]);
    grid(ax, 'on');
    box(ax, 'on');
    ax.GridAlpha = .18;
    exportgraphics(fig, fullfile(figDir, 'Fig2.pdf'), 'ContentType', 'vector');
    exportgraphics(fig, fullfile(figDir, 'Fig2.png'), 'Resolution', 300);
    close(fig);
    completed = struct('created', true, 'software', version, ...
                       'source', 'Stored MATLAB iterate data', 'smoothing', false, ...
                       'circle_radius', 1, 'bfgs_iterations', b.summary.iterations);
    PolarDFP.writeJSON(fullfile(figDir, 'figures_complete.json'), completed);
end
