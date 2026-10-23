function make_scaling_figure
% MAKE_SCALING_FIGURE Plot the stored experiments for three Holder exponents.
% SPDX-License-Identifier: MIT
colors = [0, .447, .698; .80, .239, .333; .725, .302, 0];
styles = {'-', '--', '-.'};
labels = {'$\beta=1/3$', '$\beta=1/2$', '$\beta=2/3$'};
fig = figure('Visible', 'off', 'Color', 'w', 'Units', 'inches', ...
    'Position', [0, 0, 4.75, 2.55], 'PaperPositionMode', 'auto', 'Renderer', 'painters');
cleanup = onCleanup(@()close(fig));
for panel = 1:2
    ax = axes(fig, 'Position', [.105 + (panel-1)*.49, .24, .35, .65]);
    hold(ax, 'on');
    set(ax, 'FontName', 'Nimbus Sans', 'FontSize', 8, ...
        'TickLabelInterpreter', 'latex', 'Color', 'w', 'XColor', 'k', 'YColor', 'k', ...
        'GridColor', [.5, .5, .5], 'GridAlpha', .18, 'LineWidth', .5);
    handles = gobjects(3,1);
    largest = 1;
    for k = 1:3
        saved = load(sprintf('exponent_case_%d.mat',k), 'diagnostic');
        o = saved.diagnostic;
        if panel == 1, y=o.ratio_radial; else, y=o.ratio_correction; end
        largest = max(largest, max(y));
        keep=o.cycle>0;
        handles(k) = plot(ax, o.cycle(keep), y(keep), styles{k}, ...
            'Color', colors(k,:), 'LineWidth', .95);
    end
    yline(ax, 1, ':', 'Color', [.3,.3,.3], 'LineWidth', .8, 'HandleVisibility', 'off');
    xlabel(ax, 'Cycle $j$', 'Interpreter', 'latex', 'FontSize', 8);
    ylabel(ax, 'Scaled ratio', 'Interpreter', 'latex', 'FontSize', 8);
    if panel == 1
        title(ax, '(a) Radial decrease', 'Interpreter', 'latex', 'FontSize', 8, 'FontWeight', 'normal', 'Color', 'k');
    else
        title(ax, '(b) Gradient correction', 'Interpreter', 'latex', 'FontSize', 8, 'FontWeight', 'normal', 'Color', 'k');
    end
    set(ax,'XScale','log','XMinorGrid','off');
    xlim(ax,[1,o.cycles]);
    xticks(ax,[1,1e2,1e4,1e6,1e7]);
    ylim(ax, [1-.06*(largest-1), 1+1.55*(largest-1)]);
    ax.XTickLabelRotation=0;
    lg=legend(ax, handles, labels, 'Location', 'northeast', 'Interpreter', 'latex', ...
        'FontSize', 8, 'Box', 'on');
    set(lg, 'Color', 'w', 'TextColor', 'k', 'EdgeColor', [.3,.3,.3]);
    lg.ItemTokenSize=[12,8];
    grid(ax, 'on');box(ax, 'on');
end
exportgraphics(fig, fullfile(pwd,'Fig3.pdf'), 'ContentType', 'vector');
exportgraphics(fig, fullfile(pwd,'Fig3.png'), 'Resolution', 300);
end
