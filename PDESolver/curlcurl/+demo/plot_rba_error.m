% Range on x-axis
[a, b] = deal(0, 10);
step = (b-a)/1024;

figure();

% Plot comparison of RBA and exp(-x) on (a, b)
subplot(3, 1, 1);
plot_rba(2:4, a:step:b);

% Plot error between RBA and exp(-x) on (a, b)
subplot(3, 1, 2);
plot_rba_err(2:18, a:step:b);

% Plot relative error between RBA and exp(-x) on (a, b)
subplot(3, 1, 3);
plot_rba_err_rel(2:18, a:step:b);


function plot_rba(ns, xs)
    % Loop over RBA orders
    for n = ns
        r = rba.create_rba_scalar(n);

        % Plot
        title = sprintf('r%d', n);
        semilogy(xs, r(xs), '.', 'DisplayName', title);
        hold('on');
    end

    % Plot exp
    semilogy(xs, exp(-xs), '-k', 'DisplayName', 'exp(-x)');
    legend('-DynamicLegend', 'Location', 'southwest');

    % Ignore warning
    w = warning('off', 'MATLAB:Axes:NegativeDataInLogAxis');
    shg();
    warning(w);
end


function plot_rba_err(ns, xs)
    % Loop over RBA orders
    for n = ns
        r = rba.create_rba_scalar(n);

        % Plot
        semilogy(xs, abs(exp(-xs)-r(xs)), '.');
        hold('on');
    end
end


function plot_rba_err_rel(ns, xs)
    % Loop over RBA orders
    for n = ns
        r = rba.create_rba_scalar(n);

        % Plot
        semilogy(xs, abs(exp(-xs)-r(xs))./exp(-xs), '.');
        hold('on');
    end
end
