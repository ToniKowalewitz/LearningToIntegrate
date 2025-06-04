function [Mtx, Mrx] = select_linearly_independent_measurements(Mtx, Mrx)
    % Remove linearly dependent measurements

    fprintf('Select linearly independent measurements: ');
    t = tic();
    fprintf('%d before, ', size(Mtx, 2));
    cols = app_dc.select_linearly_independent_measurements_impl(Mtx, Mrx);
    Mtx = Mtx(:, cols);
    Mrx = Mrx(:, cols);
    fprintf('%d after, ', size(Mtx, 2));
    fprintf('%f seconds\n', toc(t));

end
