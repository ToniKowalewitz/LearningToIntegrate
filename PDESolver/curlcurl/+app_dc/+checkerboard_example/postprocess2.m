plot_performance_characteristics_compare('2d', '2dfw', '2dnw');
plot_performance_characteristics_compare('3d', '3dfw', '3dnw');


function plot_performance_characteristics_compare(output_tag, tag1, tag2)

    matfile1 = load(sprintf('checkerboard-timings-%s.mat', tag1));
    matfile2 = load(sprintf('checkerboard-timings-%s.mat', tag2));
    timings1 = matfile1.timings;
    timings2 = matfile2.timings;

    t_normal1 = vertcat(timings1.t_normal);
    t_normal2 = vertcat(timings2.t_normal);
    num_dofs_inv1 = vertcat(timings1.num_dofs_inv);
    num_dofs_inv2 = vertcat(timings2.num_dofs_inv);
    num_obs1 = vertcat(timings1.num_obs);
    num_obs2 = vertcat(timings2.num_obs);

    figure();
    p = loglog(num_obs1.*num_dofs_inv1, t_normal1(:, :), 'x');
    p(1).Marker = 'x';
    p(2).Marker = '+';
    xlabel('$MN$', 'interpreter', 'latex');
    ylabel('$t_{\mathrm{norm}}$ [s]', 'interpreter', 'latex');

    hold on;

    p = loglog(num_obs2.*num_dofs_inv2, t_normal2(:, :), 'x');
    p(1).Marker = 'd';
    p(2).Marker = 's';

    [slope_x_, slope_y_] = get_slope(num_obs1.*num_dofs_inv1, t_normal1(:, end), 1);
    loglog(slope_x_, slope_y_, '--');

    x2 = num_dofs_inv1.*num_obs1;
    y2 = num_dofs_inv1.*num_obs1.^2;
    y2 = y2 * slope_y_(end) / y2(end);
    ymin = min(t_normal1(:));
    x2 = x2(y2>=ymin);
    y2 = y2(y2>=ymin);
    loglog(x2, y2, '-.');

    labels = {};
    for i = 1:size(t_normal1, 2)
        labels{end+1} = sprintf('$i=%d$ (%s)', i, tag_to_desc(tag1));  %#ok<AGROW>
    end
    for i = 1:size(t_normal2, 2)
        labels{end+1} = sprintf('$i=%d$ (%s)', i, tag_to_desc(tag2));  %#ok<AGROW>
    end
    labels{end+1} = 'slope $O(M N)$';
    labels{end+1} = 'slope $O(M^2 N)$';
    legend(labels, 'Location', 'northwest', 'interpreter', 'latex');

    saveas(gcf, sprintf('checkerboard-%s-normal.pdf', output_tag));
    export_tikz(sprintf('checkerboard-%s-normal.tex', output_tag));

end


function desc = tag_to_desc(tag)
    switch tag(3:4)
    case 'fw'
        desc = 'Algorithm~\ref{alg:gnstep-iterative}';
    case 'nw'
        desc = 'Algorithm~\ref{alg:gnstep-iterative-nw}';
    end
end


function [x_, y_] = get_slope(x, y, k)
    x_ = [min(x(:)), max(x(:))];
    y_ = x_.^k;
    y_ = y_ * max(y(:))/max(y_(:));
end


function export_tikz(filename)
    % Try exporting TikZ using matlab2tikz package or fail gracefully

    opts = {
        'showInfo'; false;
        'width'; '\figwidth';
        'height'; '\figheight';
        'parseStrings'; false;
        'extraTikzpictureOptions'; {'trim axis left', 'trim axis right'};
    };

    try
        matlab2tikz(filename, opts{:});
    catch ME
        switch ME.identifier
        case 'MATLAB:UndefinedFunction'
            warning('matlab2tikz package not installed; skipping TikZ export');
        otherwise
            rethrow(ME)
        end
    end
end
