fprintf('================Running 2D MT demo================\n');
run_inversion_demo([32, 32], [64, 64]);

fprintf('================Running 3D MT demo================\n');
% FIXME: Need better matching of mesh and anomaly/layers
%        to demonstrate actual inversion
run_inversion_demo([8, 8, 8], [12, 12, 12]);


function run_inversion_demo(mesh_resolution, mesh_resolution_synthetic)

    % Generate synthetic data
    [electrode_coords, measured_values, g, lambda1, lambda2, plot_mesh, plot_resistivity] = ...
        generate_synthetic_data(mesh_resolution_synthetic);

    % Prepare plot hook
    plot_data = create_data_plot_hooks(measured_values);

    % Build mesh
    [mesh, ~, sigma, ~, ~] = generate_mesh(mesh_resolution);

    % Report
    fprintf('Inversion mesh num cells = %d\n', mesh.num_entities(mesh.dim));
    fprintf('Inversion mesh num vertices = %d\n', mesh.num_entities(0));

    % Plot mesh
    plot_mesh(mesh);
    drawnow();

    % Build Nedelec space
    element = fe.create_nedelec_element(mesh.dim, 1);
    for d = element.get_dof_entity_dims()
        w = warning('off', 'Mesh:ExtraConnStored');  % Mute warning
        mesh.compute_connectivity(mesh.dim, d);
        warning(w);  % Restore warning
        mesh.clear_connectivity(d, 0);
    end
    dofmap = assembling.build_dofmap(mesh, element);

    % Prepare parameter space
    element_l2 = fe.create_p0_element(mesh.dim);
    dofmap_l2 = assembling.build_dofmap(mesh, element_l2);

    % Prepare transformed parameter
    % TODO: Add support for more general parameter transform
    m_ref = log(sigma);
    m = log(sigma);
    dm = zeros(size(m));

    % Prepare functions for assembling/solving (regularized) Gauss-Newton system
    switch mesh.dim
    case 2
        method = 'mixed';
    case 3
        method = 'krylov';
    end
    assemble_observation = app_mt.create_observation(dofmap, g, electrode_coords);
    solve_regularized_system = nls.create_regularization(measured_values, ...
                                                         dofmap_l2, m_ref, method);

    % Iteration parameters
    beta = 0.005;
    maxit = 4;

    fprintf('Entering Gauss-Newton loop...\n');
    fprintf('=============================\n');
    for i = 1:maxit

        % Assemble Gauss-Newton system
        [d, J] = assemble_observation(lambda1, lambda2, sigma);

        % Run Taylor test of Jacobian
        %figure(6);
        %nls.plot_taylor_test(J, m, d, @(m) assemble_observation(lambda1, lambda2, exp(m)));

        % Solve normal equations and update parameters
        dm(:) = solve_regularized_system(d, J, m, beta);
        m(:) = m + dm;
        sigma(:) = exp(m);

        misfit(i) = norm(d - measured_values);  %#ok<AGROW>

        % Report
        fprintf('i = %d, Beta = %f, ||dm|| = %f\n', i, beta, norm(dm, 2));
        fprintf('Misfit = %f\n', misfit(i));
        fprintf('Min sigma, max sigma: %f %f\n', min(sigma), max(sigma));

        % Update regularization
        if mod(i, 3) == 0
            beta = 0.5*beta;
        end

        % Plot current value of resistivity
        plot_resistivity(dofmap_l2, sigma, i);

        % Plot current value of data and misfit
        plot_data(d, misfit, i);

        % Update plots
        drawnow();
    end

end


function [electrode_coords, d, g, lambda1, lambda2, plot_mesh, plot_resistivity] = ...
        generate_synthetic_data(mesh_resolution)

    fprintf('Generating synthetic data...\n');
    t = tic();

    % Build mesh and true resistivity
    [mesh, sigma, ~, layers_rho, layers_d] = generate_mesh(mesh_resolution);

    % Prepare electrode coords
    num_electrodes_1d = 4;
    switch mesh.dim
    case 2
        electrode_coords = create_electrode_coords_2d(num_electrodes_1d);
    case 3
        electrode_coords = create_electrode_coords_3d(num_electrodes_1d);
    end

    % Prepare plot hooks
    switch mesh.dim
    case 2
        [plot_true, plot_mesh, plot_resistivity] = ...
            create_resistivity_plot_hooks_2d(electrode_coords);
    case 3
        [plot_true, plot_mesh, plot_resistivity] = ...
            create_resistivity_plot_hooks_3d();
    end

    % Plot resistivity
    plot_true(mesh, sigma);

    % Report
    fprintf('Synthetic mesh num cells = %d\n', mesh.num_entities(mesh.dim));
    fprintf('Synthetic mesh num vertices = %d\n', mesh.num_entities(0));

    % Build Nedelec space
    element = fe.create_nedelec_element(mesh.dim, 1);
    for d = element.get_dof_entity_dims()
        w = warning('off', 'Mesh:ExtraConnStored');  % Mute warning
        mesh.compute_connectivity(mesh.dim, d);
        warning(w);  % Restore warning
        mesh.clear_connectivity(d, 0);
    end
    dofmap = assembling.build_dofmap(mesh, element);

    % Problem parameters
    mu = 1.25663706212e-6;   % vacuum permeability [H/m]
    eps = 8.8541878128e-12;  % vacuum permittivity [F/m]
    switch mesh.dim
    case 2
        freq = 10 .^ linspace(0, 3, 10);  % frequencies [Hz]
    case 3
        freq = [1e3, 1e4];                % frequencies [Hz]
    end
    omega = 2*pi*freq;
    lambda1 = 1j*omega*mu;
    lambda2 = -omega.^2*mu*eps;

    % Prepare exact solutions for layered space
    g = get_reference_solutions(mesh.dim, freq, mu, layers_rho, layers_d);

    % Compute observation
    assemble_observation = app_mt.create_observation(dofmap, g, electrode_coords);
    d = assemble_observation(lambda1, lambda2, sigma);

    fprintf('Generated synthetic data: %f seconds\n', toc(t));
end


function g = get_reference_solutions(dim, frequencies, mu, layers_rho, layers_d)

    for j = 1:numel(frequencies)

        E_exact = @(z) app_mt.layered_space_plane_wave(frequencies(j), mu, layers_rho, layers_d, z);

        switch dim
        case 2
            g{1, j} = @(x, c) [E_exact(x(dim)), 0];  %#ok<AGROW>
        case 3
            g{1, j} = @(x, c) [E_exact(x(dim)), 0, 0];  %#ok<AGROW>
            g{2, j} = @(x, c) [0, E_exact(x(dim)), 0];  %#ok<AGROW>
        end

    end
end


function [mesh, sigma, sigma_ref, rho, d] = generate_mesh(mesh_resolution)

    % Layers of reference model
    rho = [1e14, 1e2, 1e5, 1e2];    % resistivities in layers [Ohm.m]
    d = [-inf, 0, 100, 250, +inf];  % depths of layers [m]

    % Create geometry/mesh for cube (-1000, 1000)^dim
    mesh = meshing.generate_unit_cube_mesh(mesh_resolution);
    shift = 0.5*ones(mesh.dim, 1);
    scale = 2e3;
    mesh.set_coordinates(scale*(mesh.vertex_coords - shift));

    % Construct cell markers cm from depths d
    cm = zeros(mesh.num_entities(mesh.dim), 1);
    coords = mesh.get_cell_centroids();
    z = coords(mesh.dim, :);
    for layer = 2:numel(d)
        cm(z >= d(layer-1) & z < d(layer)) = layer-1;
    end
    clear('z');
    assert(all(cm > 0));

    % Create anomaly
    z = coords(mesh.dim, :);
    r = vecnorm(coords);
    assert(numel(r) == mesh.num_entities(mesh.dim));
    rho_true = [rho, 1e1];
    anomaly_radius = 100;
    anomaly_marker_val = numel(rho_true);
    cm_with_anomaly = cm;
    cm_with_anomaly(z >= 0 & r < anomaly_radius) = anomaly_marker_val;
    clear('z', 'r');

    % Clean up
    clear('coords');
    mesh.clear_boundary_facets();

    % Reference value of sigma
    sigma = 1.0 ./ rho_true(cm_with_anomaly);
    sigma_ref = 1.0 ./ rho_true(cm);

    % We assume sigma is column
    sigma = sigma(:);
    sigma_ref = sigma_ref(:);
end


function coords = create_electrode_coords_1d(num_electrodes)
    coords = linspace(-50, 50, num_electrodes);
end


function coords = create_electrode_coords_2d(num_electrodes_1d)
    X = create_electrode_coords_1d(num_electrodes_1d);
    assert(isrow(X));
    coords = zeros(2, numel(X));
    coords(1, :) = X;
end


function coords = create_electrode_coords_3d(num_electrodes_1d)
    X = create_electrode_coords_1d(num_electrodes_1d);
    assert(isrow(X));
    n = numel(X);
    I = ones(1, n);
    coords = zeros(3, n^2);
    coords(1, :) = kron(I, X);
    coords(2, :) = kron(X, I);
end


function [p_true, p_mesh, p_current] = create_resistivity_plot_hooks_2d(electrode_coords)
    fig = figure(); subplot(1, 2, 1);

    function plot_true(mesh, sigma)
        set(0, 'CurrentFigure', fig);
        subplot(1, 2, 1);
        title('True conductivity');
        meshing.plot_mesh(mesh, 'cell_markers', sigma);
        set(gca(), 'colorscale', 'log');
        hold('on');
        plot(electrode_coords(1, :), electrode_coords(2, :), '.');
        hold('off');
    end

    function plot_mesh(mesh)
        figure();
        meshing.plot_mesh(mesh);
        savefig('mt2d-mesh.fig');
    end

    function plot_current(dofmap, sigma, i)
        set(0, 'CurrentFigure', fig);
        subplot(1, 2, 1); cax = caxis();
        subplot(1, 2, 2);
        cla();
        title(sprintf('Conductivity at iteration %d', i));
        meshing.plot_mesh(dofmap.mesh, 'cell_markers', sigma);
        set(gca(), 'colorscale', 'log');
        caxis(cax);
        hold('on');
        plot(electrode_coords(1, :), electrode_coords(2, :), '.');
        hold('off');

        savefig(sprintf('mt2d-conductivity-%04d.fig', i));
    end

    p_true = @plot_true;
    p_mesh = @plot_mesh;
    p_current = @plot_current;
end


function [p_true, p_mesh, p_current] = create_resistivity_plot_hooks_3d()

    xdmf = io.XDMF('mt3d-conductivity.xdmf');

    function plot_true(mesh, sigma)
        element = fe.create_p0_element(mesh.dim);
        dofmap = assembling.build_dofmap(mesh, element);
        xdmf_ = io.XDMF('mt3d-conductivity-true.xdmf');
        xdmf_.write(dofmap, sigma, 0);
    end

    function plot_mesh(~)
    end

    function plot_current(dofmap, sigma, i)
        xdmf.write(dofmap, sigma, i);
        xdmf.flush();
    end

    p_true = @plot_true;
    p_mesh = @plot_mesh;
    p_current = @plot_current;
end


function plot_data = create_data_plot_hooks(measured_values)
    fig = figure();
    subplot(2, 1, 1);

    function plot_data_(d, misfit, i)
        set(0, 'CurrentFigure', fig);
        clf();

        subplot(2, 1, 1);
        hold('on');
        plot(util.complex2real(measured_values), 'xr');
        plot(util.complex2real(d), 'ob');
        hold('off');
        title(sprintf('Data at iteration %d', i-1));
        legend({'measured', 'modeled'}, 'Location', 'bestoutside');

        subplot(2, 1, 2);
        semilogy(0:numel(misfit)-1, misfit, 'x-b');
        maxit = max(1, 10 * (1 + floor((numel(misfit)-2)/10)));
        xlim([0, maxit]);
        title(sprintf('Misfit at iteration %d', i-1));

        savefig('mt_fit.fig');
    end

    plot_data = @plot_data_;
end
