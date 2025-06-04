compose_resistivity_figures_2d();
plot_performance_characteristics('2dfw');
plot_performance_characteristics('2dnw');
plot_performance_characteristics('3dfw');
plot_performance_characteristics('3dnw');


function compose_resistivity_figures_2d()

    files = dir('checkerboard-resistivity-2d*.fig');
    for i = 1:numel(files)
        f = openfig(files(i).name, 'invisible');

        fig = figure();
        axis('equal', 'tight');
        copyobj(allchild(f.Children(6)), gca);
        set(gca, 'ColorScale', 'log');
        colorbar();
        [~, basename, ~] = fileparts(files(i).name);
        saveas(fig, strcat(basename, '-true.pdf'));
        close(fig);

        fig = figure();
        axis('equal', 'tight');
        copyobj(allchild(f.Children(2)), gca);
        set(gca, 'ColorScale', 'log');
        colorbar();
        [~, basename, ~] = fileparts(files(i).name);
        saveas(fig, strcat(basename, '-inv.pdf'));
        close(fig);
    end

end


function plot_performance_characteristics(tag)

    matfile = load(sprintf('checkerboard-timings-%s.mat', tag));
    timings = matfile.timings;

    t_observe = vertcat(timings.t_observe);  %#ok<NASGU>
    t_normal = vertcat(timings.t_normal);
    t_woodbury1 = vertcat(timings.t_woodbury1);
    t_woodbury2 = vertcat(timings.t_woodbury2);
    t_woodbury3 = vertcat(timings.t_woodbury3);
    iter_normal = vertcat(timings.iter_normal);
    num_dofs_fw = vertcat(timings.num_dofs_fw);  %#ok<NASGU>
    num_dofs_inv = vertcat(timings.num_dofs_inv);
    num_obs = vertcat(timings.num_obs);
    num_electrodes = vertcat(timings.num_electrodes);
    misfit = vertcat(timings.misfit);  %#ok<NASGU>
    beta = vertcat(timings.beta);
    n = vertcat(timings.n);

    table = [n, num_electrodes, num_dofs_inv, num_obs, ...
             t_woodbury1, t_woodbury2, t_woodbury3, iter_normal, t_normal];
    filename = sprintf('checkerboard-timings-%s.tex', tag);
    write_latex(filename, table);

    [num_rows, num_cols] = size(t_normal);
    if size(n, 2) == 1
        n = [n, ones(size(n))];
    end
    dup = @(arr) repelem(arr, num_cols, 1);
    iter = repmat(transpose(1:num_cols), num_rows, 1);
    flat = @(arr) reshape(transpose(arr), [], 1);
    table = [dup(n), dup(num_electrodes), dup(num_dofs_inv), dup(num_obs), ...
             iter, ...
             flat(t_woodbury1), flat(t_woodbury2), flat(t_woodbury3), ...
             flat(iter_normal), flat(t_normal)];
    filename = sprintf('checkerboard-timings-flat-%s.tex', tag);
    write_latex(filename, table);

    assert(size(unique(beta, 'rows'), 1) == 1);
    beta = beta(1, :);
    labels = split(strip(sprintf('i=%d (\\beta=%f)  ', [1:numel(beta); beta])), '  ');

    [slope_x_, slope_y_] = get_slope(num_obs.*num_dofs_inv, t_normal(:, end), 1);
    figure();
    loglog(num_obs.*num_dofs_inv, t_normal(:, :), 'x', slope_x_, slope_y_, '--');
    xlabel('MN');
    ylabel('time for solving normal equations [s]');
    legend(vertcat(labels, 'C M N'), 'Location', 'northwest');
    saveas(gcf, sprintf('checkerboard-%s-normal.pdf', tag));

    [slope_x_, slope_y_] = get_slope(num_obs.*num_dofs_inv, t_woodbury1(:, end), 1);
    figure();
    loglog(num_obs.*num_dofs_inv, t_woodbury1(:, :), 'x', slope_x_, slope_y_, '--');
    xlabel('MN');
    ylabel('time for computing H for Woodbury correction [s]');
    legend(vertcat(labels, 'C M N'), 'Location', 'northwest');
    saveas(gcf, sprintf('checkerboard-%s-woodbury1.pdf', tag));

    [slope_x_, slope_y_] = get_slope(num_obs.^2.*num_dofs_inv, t_woodbury2(:, end), 1);
    figure();
    loglog(num_obs.^2.*num_dofs_inv, t_woodbury2(:, :), 'x', slope_x_, slope_y_, '--');
    xlabel('M^2 N');
    ylabel('time for assembling capacitance matrix for Woodbury correction [s]');
    legend(vertcat(labels, 'C M^2 N'), 'Location', 'northwest');
    saveas(gcf, sprintf('checkerboard-%s-woodbury2.pdf', tag));

    [slope_x_, slope_y_] = get_slope(num_obs, t_woodbury3(:, end), 3);
    figure();
    loglog(num_obs, t_woodbury3(:, :), 'x', slope_x_, slope_y_, '--');
    xlabel('M');
    ylabel('time for computing Cholesky factor for Woodbury correction [s]');
    legend(vertcat(labels, 'C M^3'), 'Location', 'northwest');
    saveas(gcf, sprintf('checkerboard-%s-woodbury3.pdf', tag));

    figure();
    semilogx(num_obs, iter_normal(:, :), 'x');
    xlabel('M');
    ylabel('number MINRES iterations');
    legend(labels, 'Location', 'northwest');
    saveas(gcf, sprintf('checkerboard-%s-minres.pdf', tag));

end


function [x_, y_] = get_slope(x, y, k)
    x_ = [min(x(:)), max(x(:))];
    y_ = x_.^k;
    y_ = y_ * max(y(:))/max(y_(:));
end


function write_latex(filename, arr)

    % Generate latex code
    prev = sympref('FloatingPointOutput', true);
    code = latex(sym(arr));
    sympref('FloatingPointOutput', prev);

    % Strip all prolog and epilog code (keep just numbers)
    code = replace(code, '\left(', '');
    code = replace(code, '\right)', '');
    code = regexprep(code, '\\begin{array}{c*}', '');
    code = replace(code, '\end{array}', '');
    code = strip(code);

    % Write to file
    [fid, msg] = fopen(filename, 'wt');
    if fid == -1
        error('Opening "%s" for writing failed with reason:\n%s', filename, msg);
    end
    fwrite(fid, code, 'char');
    fclose(fid);

end
