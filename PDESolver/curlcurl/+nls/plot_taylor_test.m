function plot_taylor_test(J, x0, f0, assemble_f)
    % Run (visual) Taylor test for derivative of a function
    %
    % For multiple, randomly chosen, direction vectors dx (of size
    % ||dx|| = ||x0||) evaluates the residual
    %
    %     ||assemble_f(x0+t*dx) - f0 - t*J*dx||
    %
    % for t = 1, 1/2, 1/4, 1/8, ..., and verifies (by plotting) that
    % the observed residual decays as C*t^2.
    %
    % INPUT PARAMETERS
    %   J          ... matrix [M, N]
    %   x0         ... vector [N, 1]
    %   f0         ... vector [M, 1]
    %   assemble_f ... function handle of signature:
    %
    %                      f = assemble_f(x);
    %
    %                  x ... vector [N, 1]
    %                  f ... vector [M, 1]

    % NB: For now hard-coded parameters
    num_directions = 3;
    num_steps = 8;

    scaling = vecnorm(x0);
    result = zeros(num_directions*num_steps, 2);
    for i = 1:num_directions
        dx = rand(size(x0));
        dx = scaling/vecnorm(dx)*dx;
        df = J*dx;
        t = 1;
        for j = 1:num_steps
            residual = assemble_f(x0+t*dx) - f0 - t*df;
            result((i-1)*num_steps + j, 1:2) = [t, vecnorm(residual)];
            t = 0.5*t;
        end
    end

    % Plot measured rate and reference quadratic rate
    plot_convergence(result(:, 1), result(:, 2));
end


function plot_convergence(x, y)
    [x_, y_] = get_slope(x, y, 2);
    loglog(x, y, '.', x_, y_, '--');
    set(gca, 'XDir', 'reverse');
    title('Taylor test');
    xlabel('t')
    legend({'||f(x0 + t*dx) - f0 - t*J*dx||', 'Ct^2'});
end


function [x_, y_] = get_slope(x, y, k)
    x_ = [min(x(:)), max(x(:))];
    y_ = x_.^k;
    y_ = y_ * max(y(:))/max(y_(:));
end
