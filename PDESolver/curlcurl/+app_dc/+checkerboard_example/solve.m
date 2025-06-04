function timings = solve(dim, n, ref, ref_synth, tag, inv_solver_type)
    % Run electrical resistivity tomography demo with checkerboard pattern
    %
    % INPUT PARAMETERS
    %   dim             ... spatial dimension (2 or 3)
    %   n               ... checkerboard pattern, positive integer (dim=2)
    %                       or Vector[2, 1] of positive integers (dim=3)
    %   ref             ... inversion mesh refinement level
    %   ref_synth       ... synthetic mesh refinement level
    %   tag             ... basename for output files
    %   inv_solver_type ... type of regularization solver, e.g.,
    %                       'diag', 'mixed', 'krylov', or 'krylov-no-woodbury'
    %
    % OUTPUT PARAMETERS
    %   timings   ... struct array with various peformance data,
    %                 useful in postprocessing

    % Refinement parameters
    refinement_level = ref;
    refinement_level_synthetic = ref_synth;

    % Generate synthetic data
    [electrode_coords, Mtx, Mrx, sigma_ref, measured_values, plot_mesh, plot_resistivity] = ...
        generate_synthetic_data(dim, n, refinement_level_synthetic, tag);

    % Prepare plot hook
    plot_data = create_data_plot_hooks(measured_values, tag);

    % Build mesh
    switch dim
    case 2
        num_electrodes_1d = 1+4*n;
        h = 100/(num_electrodes_1d-1);
        [mesh, fm, pn] = meshing.generate_mesh2D(...
            'domain_c', [0, 0], ...
            'domain_r', 80, ...
            'point', electrode_coords.', ...
            'size_at_pt', h, ...
            'ref', refinement_level, ...
            'marker', [dim-1, -1]);
    case 3
        h = min(100./(4*n));
        [mesh, fm, pn] = meshing.generate_mesh3D(...
            'domain_c', [0, 0, 0], ...
            'domain_r', 80, ...
            'point', electrode_coords.', ...
            'size_at_pt', h, ...
            'ref', refinement_level, ...
            'marker', [dim-1, -1]);
    end
    bc_facet_tag = pn{1}{1}(ismember(pn{1}{2}, {'subsurface'}));

    % Report
    fprintf('Inversion mesh num cells = %d\n', mesh.num_entities(mesh.dim));
    fprintf('Inversion mesh num vertices = %d\n', mesh.num_entities(0));
    [bbmin, bbmax] = bounds(mesh.vertex_coords, 2);
    fprintf('Inversion mesh bounding box: %s\n', mat2str([bbmin, bbmax]));

    % Plot mesh
    plot_mesh(mesh, fm);
    drawnow();

    % Prepare forward solution space
    element = fe.create_lagrange_element(mesh.dim, 1);
    dofmap = assembling.build_dofmap(mesh, element);

    % Prepare parameter space
    element_l2 = fe.create_p0_element(mesh.dim);
    dofmap_l2 = assembling.build_dofmap(mesh, element_l2);

    % Prepare physical and transformed parameters
    % TODO: Add support for more general parameter transform
    sigma_ref = get_reference_parameters(dofmap_l2, sigma_ref);
    sigma = get_initial_parameters(dofmap_l2, sigma_ref);
    m_ref = log(sigma_ref);
    m = log(sigma);
    dm = zeros(size(m));

    % Prepare functions for assembling/solving (regularized) Gauss-Newton system
    switch dim
    case 2
        fwd_solver_type = 'backslash';
    case 3
        fwd_solver_type = 'backslash';
    end
    assemble_observation = app_dc.create_observation(electrode_coords, Mtx, Mrx, ...
                                                     dofmap, fm, bc_facet_tag, ...
                                                     fwd_solver_type);
    solve_regularized_system = nls.create_regularization(measured_values, ...
                                                         dofmap_l2, m_ref, ...
                                                         inv_solver_type);

    % Iteration parameters
    switch dim
    case 2
        beta = 1e-3;
    case 3
        beta = 1e3;
    end
    maxit = 2;

    t_observe = [];
    t_normal = [];
    iter_normal = [];
    t_woodbury1 = [];
    t_woodbury2 = [];
    t_woodbury3 = [];

    fprintf('Entering Gauss-Newton loop...\n');
    fprintf('=============================\n');
    for i = 1:maxit

        % Assemble Gauss-Newton system
        t = tic();
        [d, J] = assemble_observation(sigma);
        t_observe(i) = toc(t);  %#ok<AGROW>

        %% Dump sensitivity to file
        %plot_sensitivity(mesh, J, sigma, sprintf('%s-iter%02d', tag, i));

        %% Run Taylor test of Jacobian
        %figure(double(intmax-6));
        %nls.plot_taylor_test(J, m, d, @(m) assemble_observation(exp(m)));

        % Solve normal equations and update parameters
        t = tic();
        [dm(:), iter_normal(i), t_woodbury1(i), t_woodbury2(i), t_woodbury3(i)] = ...
            solve_regularized_system(d, J, m, beta);  %#ok<AGROW>
        m(:) = m + dm;
        sigma(:) = exp(m);
        t_normal(i) = toc(t);  %#ok<AGROW>

        misfit(i) = norm(d - measured_values);  %#ok<AGROW>

        % Report
        fprintf('i = %d, Beta = %f, ||dm|| = %f\n', i, beta, norm(dm, 2));
        fprintf('Misfit = %f\n', misfit(i));
        fprintf('Min sigma, max sigma: %f %f\n', min(sigma), max(sigma));

        % Plot current value of resistivity
        plot_resistivity(dofmap_l2, sigma, i);

        % Plot current value of data and misfit
        plot_data(d, misfit, i);

        % Update plots
        drawnow();

        beta_(i) = beta;  %#ok<AGROW>

        % Update regularization
        if mod(i, 3) == 0
            beta = 0.5*beta;
        end
    end

    timings = struct();
    timings.n = n;
    timings.num_electrodes = size(electrode_coords, 2);
    timings.num_dofs_fw = dofmap.dim;
    timings.num_dofs_inv = dofmap_l2.dim;
    timings.num_obs = numel(d);
    timings.t_observe = t_observe;
    timings.t_normal = t_normal;
    timings.t_woodbury1 = t_woodbury1;
    timings.t_woodbury2 = t_woodbury2;
    timings.t_woodbury3 = t_woodbury3;
    timings.iter_normal = iter_normal;
    timings.misfit = misfit;
    timings.beta = beta_;
end


function [electrode_coords, Mtx, Mrx, sigma_ref, d, plot_mesh, plot_resistivity] ...
        = generate_synthetic_data(dim, n, refinement_level, tag)
    % Generate synthetic observations for actual inversion
    %
    % INPUT PARAMETERS
    %   dim              ... spatial dimension (2 or 3)
    %   n                ... checkerboard pattern, positive integer (dim=2)
    %                        or vector [2, 1] of positive integers (dim=3)
    %   refinement_level ... integer, synthetic mesh refinement level
    %   tag              ... basename for output files
    %
    % OUTPUT PARAMETERS
    %   electrode_coords ... matrix [dim, num_electrodes], coordinates of electrodes
    %   Mtx, Mrx         ... matrix [num_electrodes, num_observations], the two matrices
    %                        define observations (apparent resistivities) in terms of
    %                        potential point values at tx electrodes for point sources
    %                        at rx electrodes
    %   sigma_ref        ... positive scalar, reference (background) conductivity
    %   d                ... vector [num_observations, 1], observed apparent resistivities
    %   plot_mesh        ... function handle, hook for plotting (see the code for details)
    %   plot_resistivity ... function handle, hook for plotting (see the code for details)

    fprintf('Generating synthetic data...\n');
    t = tic();

    % Create configuration of electrodes
    switch dim
    case 2
        num_electrodes_1d = 1+4*n;
        [Mtx, Mrx] = app_dc.create_electrode_configuration('pdp', [2, 4, 8], num_electrodes_1d);
        electrode_coords = create_electrode_coords_1d(num_electrodes_1d);
    case 3
        num_electrodes_1d = 1+4*max(n);
        [Mtx, Mrx] = app_dc.create_electrode_configuration('pdp', [2, 4, 8], num_electrodes_1d);
        [Mtx, Mrx] = create_electrode_configuration_2d_from_1d(Mtx, Mrx);
        electrode_coords = create_electrode_coords_2d(num_electrodes_1d);
    end
    [Mtx, Mrx] = app_dc.transform_to_apparent_resistivities(electrode_coords, Mtx, Mrx);

    % Build mesh
    switch dim
    case 2
        nx = n;
        ny = 1;
        h = 100/(num_electrodes_1d-1);
        [box_c, box_dx, box_dy] = checkerboard_coords_2d(nx, ny);
        [mesh, cm, fm, pn] = meshing.generate_checkerboard2D( ...
            'block', [box_c.', box_dx.', box_dy.'], ...
            'domain_r', 80, ...
            'point', electrode_coords.', ...
            'size_at_pt', h, ...
            'refinement', refinement_level, ...
            'marker', [dim, dim-1, -1]);
    case 3
        nx = n(1);
        ny = n(2);
        nz = 1;
        h = 50/(num_electrodes_1d-1);
        [box_c, box_dx, box_dy, box_dz] = checkerboard_coords_3d(nx, ny, nz);
        [mesh, cm, fm, pn] = meshing.generate_checkerboard3D( ...
            'block', [box_c.', box_dx.', box_dy.', box_dz.'], ...
            'domain_r', 80, ...
            'point', electrode_coords.', ...
            'size_at_pt', h, ...
            'refinement', refinement_level, ...
            'marker', [dim, dim-1, -1]);
    end
    bc_facet_tag = pn{2}{1}(ismember(pn{2}{2}, {'subsurface'}));

    % Divide cells into two groups: (1) background, (2) anomaly
    switch dim
    case 2
        % Odd cell markers are background, even cell markers are anomaly
        cm(mod(cm, 2) == 1) = 1;  % background
        cm(mod(cm, 2) == 0) = 2;  % anomaly
    case 3
        cm(cm ~= 1) = 2;  % anomaly
    end

    % Report
    fprintf('Synthetic mesh num cells = %d\n', mesh.num_entities(mesh.dim));
    fprintf('Synthetic mesh num vertices = %d\n', mesh.num_entities(0));
    [bbmin, bbmax] = bounds(mesh.vertex_coords, 2);
    fprintf('Synthetic mesh bounding box: %s\n', mat2str([bbmin, bbmax]));

    % Parameters
    rho_background = 3500;
    rho_anomaly = 7000;
    sigma = zeros(size(cm));
    sigma(cm==1) = 1/rho_background;
    sigma(cm==2) = 1/rho_anomaly;
    sigma_ref = 1/rho_background;

    % Prepare plot hooks
    switch dim
    case 2
        y_min = min(box_c(2, :)-2*box_dy);
        xlims = [-50, 50];
        ylims = [y_min, 0];
        [plot_true, plot_mesh, plot_resistivity] = ...
            create_resistivity_plot_hooks_2d(electrode_coords, xlims, ylims, tag);
    case 3
        [plot_true, plot_mesh, plot_resistivity] = ...
            create_resistivity_plot_hooks_3d(0, 0, 0, tag);
    end

    % Plot resistivity
    plot_true(mesh, sigma);

    % Build function space
    element = fe.create_lagrange_element(mesh.dim, 1);
    dofmap = assembling.build_dofmap(mesh, element);
    w = warning('off', 'Mesh:ExtraConnStored');  % Mute warning
    for d = element.get_dof_entity_dims()
        mesh.compute_connectivity(mesh.dim, d);
        mesh.clear_connectivity(d, 0);
    end
    warning(w);  % Restore warning

    % Compute observation
    switch dim
    case 2
        fwd_solver_type = 'backslash';
    case 3
        fwd_solver_type = 'cg/amg';
    end
    assemble_observation = app_dc.create_observation(electrode_coords, Mtx, Mrx, ...
                                                     dofmap, fm, bc_facet_tag, ...
                                                     fwd_solver_type);
    d = assemble_observation(sigma);

    fprintf('Generated synthetic data: %f seconds\n', toc(t));
end


function [box_c, box_dx, box_dy, pts] = checkerboard_coords_2d(nx, ny)
    xmin = -40;
    xmax = +40;
    dx = (xmax - xmin)/nx;
    dy = dx;
    ymax = -ny*dy;
    ymin = ymax - ny*dy;
    X = linspace(xmin, xmax, nx+1);
    Y = linspace(ymin, ymax, ny+1);
    pts = [
        kron(X, ones(size(Y)));
        kron(ones(size(X)), Y);
    ];
    dX = X(2:end)-X(1:end-1);
    dY = Y(2:end)-Y(1:end-1);
    X = X(1:end-1);
    Y = Y(1:end-1);
    box_c = zeros(2, nx*ny);
    box_dx = zeros(1, nx*ny);
    box_dy = zeros(1, nx*ny);
    box_c(1, :) = kron(ones(1, ny), X);
    box_c(2, :) = kron(Y, ones(1, nx));
    box_dx(1, :) = kron(ones(1, ny), dX);
    box_dy(1, :) = kron(dY, ones(1, nx));
end


function [box_c, box_dx, box_dy, box_dz, pts] = checkerboard_coords_3d(nx, ny, nz)
    xmin = -40;
    xmax = +40;
    ymin = -40;
    ymax = +40;
    dx = (xmax - xmin)/nx;
    dy = (ymax - ymin)/ny;
    dz = min(dx, dy);
    dz = 0.25*dz;
    zmax = 0;
    zmin = -nz*dz;
    X = linspace(xmin, xmax, nx+1);
    Y = linspace(ymin, ymax, ny+1);
    Z = linspace(zmin, zmax, nz+1);
    pts = [
        kron(X, ones(1, numel(Y)*numel(Z)));
        kron(kron(ones(1, numel(X)), Y), ones(1, numel(Z)));
        kron(ones(1, numel(X)*numel(Y)), Z);
    ];
    dX = X(2:end)-X(1:end-1);
    dY = Y(2:end)-Y(1:end-1);
    dZ = Z(2:end)-Z(1:end-1);
    dX = 0.5*dX;
    dY = 0.5*dY;
    dZ = 0.5*dZ;
    X = X(1:end-1);
    Y = Y(1:end-1);
    Z = Z(1:end-1);
    box_c = zeros(3, nx*ny*nz);
    box_dx = zeros(1, nx*ny*nz);
    box_dy = zeros(1, nx*ny*nz);
    box_dz = zeros(1, nx*ny*nz);
    box_c(1, :) = kron(ones(1, nz*ny), X);
    box_c(2, :) = kron(ones(1, nz), kron(Y, ones(1, nx)));
    box_c(3, :) = kron(Z, ones(1, ny*nx));
    box_dx(1, :) = kron(ones(1, nz*ny), dX);
    box_dy(1, :) = kron(ones(1, nz), kron(dY, ones(1, nx)));
    box_dz(1, :) = kron(dZ, ones(1, ny*nx));
end


function coords = create_electrode_coords_2d(num_electrodes_1d)
    X = create_electrode_coords_(num_electrodes_1d);
    assert(isrow(X));
    n = numel(X);
    I = ones(1, n);
    coords = zeros(3, n^2);
    coords(1, :) = kron(I, X);
    coords(2, :) = kron(X, I);
end


function coords = create_electrode_coords_1d(num_electrodes)
    coords = zeros(2, num_electrodes);
    coords(1, :) = create_electrode_coords_(num_electrodes);
end


function coords = create_electrode_coords_(num_electrodes)
    coords = linspace(-50, 50, num_electrodes);
end


function [Mtx, Mrx] = create_electrode_configuration_2d_from_1d(Mtx, Mrx)
    % Create 2D tensor-product configuration from 1D configuration
    assert(all(size(Mtx) == size(Mrx)));
    I = speye(size(Mtx, 1));
    Mtx = horzcat(kron(Mtx, I), kron(I, Mtx));
    Mrx = horzcat(kron(Mrx, I), kron(I, Mrx));
    assert(all(size(Mtx) == size(Mrx)));
end


function sigma_init = get_initial_parameters(dofmap_l2, sigma)
    sigma_init = zeros(dofmap_l2.dim, 1);
    sigma_init(:) = sigma;
end


function sigma_ref = get_reference_parameters(dofmap_l2, sigma)
    sigma_ref = zeros(dofmap_l2.dim, 1);
    sigma_ref(:) = sigma;
end


function [p_true, p_mesh, p_current] = create_resistivity_plot_hooks_2d(electrode_coords, xlims, ylims, tag)
    fig1 = figure(); subplot(3, 1, 1);
    fig2 = figure(); subplot(3, 1, 1);

    function plot_true(mesh, sigma)
        figure(fig1);
        subplot(3, 1, 1);
        title('True resistivity');
        meshing.plot_mesh(mesh, 'cell_markers', 1./sigma);
        xlim(xlims);
        ylim(ylims);
        hold('on');
        plot(electrode_coords(1, :), electrode_coords(2, :), '.');
        hold('off');

        figure(fig2);
        subplot(3, 1, 1);
        title('True resistivity');
        meshing.plot_mesh(mesh, 'cell_markers', 1./sigma);
    end

    function plot_mesh(mesh, fm)
        figure(fig1);
        subplot(3, 1, 2);
        meshing.plot_mesh(mesh);
        xlim(xlims);
        ylim(ylims);

        figure(fig2);
        subplot(3, 1, 2);
        meshing.plot_mesh(mesh, 'facet_markers', fm);
    end

    function plot_current(dofmap, sigma, i)
        figure(fig1);
        subplot(3, 1, 3);
        cla();
        meshing.plot_mesh(dofmap.mesh, 'cell_markers', 1./sigma);
        cax = caxis();
        title(sprintf('Resistivity at iteration %d', i));
        xlim(xlims);
        ylim(ylims);
        hold('on');
        plot(electrode_coords(1, :), electrode_coords(2, :), '.');
        hold('off');
        subplot(3, 1, 1);
        caxis(cax);

        savefig(sprintf('checkerboard-resistivity-detail-%s.fig', tag));

        figure(fig2);
        subplot(3, 1, 3);
        cla();
        meshing.plot_mesh(dofmap.mesh, 'cell_markers', 1./sigma);
        cax = caxis();
        title(sprintf('Resistivity at iteration %d', i));
        subplot(3, 1, 1);
        caxis(cax);

        savefig(sprintf('checkerboard-resistivity-%s.fig', tag));
    end

    p_true = @plot_true;
    p_mesh = @plot_mesh;
    p_current = @plot_current;
end


function [p_true, p_mesh, p_current] = create_resistivity_plot_hooks_3d(~, ~, ~, tag)

    xdmf = io.XDMF(sprintf('checkerboard-resistivity-%s.xdmf', tag));

    function plot_true(mesh, sigma)
        element = fe.create_p0_element(mesh.dim);
        dofmap = assembling.build_dofmap(mesh, element);
        xdmf_ = io.XDMF(sprintf('checkerboard-resistivity-true-%s.xdmf', tag));
        xdmf_.write(dofmap, 1./sigma, 0);
    end

    function plot_mesh(mesh, fm)  %#ok<INUSD>
    end

    function plot_current(dofmap, sigma, i)
        xdmf.write(dofmap, 1./sigma, i);
        xdmf.flush();
    end

    p_true = @plot_true;
    p_mesh = @plot_mesh;
    p_current = @plot_current;
end


function plot_data = create_data_plot_hooks(measured_values, tag)
    fig = figure();
    subplot(2, 1, 1);

    function plot_data_(d, misfit, i)
        figure(fig);
        clf();

        subplot(2, 1, 1);
        hold('on');
        plot(measured_values, 'xr');
        plot(d, 'ob');
        hold('off');
        title(sprintf('Data at iteration %d', i-1));
        legend({'measured', 'modeled'}, 'Location', 'bestoutside');

        subplot(2, 1, 2);
        semilogy(0:numel(misfit)-1, misfit, 'x-b');
        maxit = max(1, 10 * (1 + floor((numel(misfit)-2)/10)));
        xlim([0, maxit]);
        title(sprintf('Misfit at iteration %d', i-1));

        savefig(sprintf('checkerboard-data-%s.fig', tag));
    end

    plot_data = @plot_data_;
end


function plot_sensitivity(mesh, J, sigma, tag)  %#ok<DEFNU>
    xdmf = io.XDMF(sprintf('checkerboard-sensitivity-%s.xdmf', tag));
    element = fe.create_p0_element(mesh.dim);
    dofmap = assembling.build_dofmap(mesh, element);

    % Remove scaling due chain rule
    % FIXME: This is prone to cause bugs
    J = J ./ (sigma(:).');

    % Scale from cell integral into pointwise quantity
    vols = mesh.get_cell_volumes();
    J = J ./ (vols(:).');
    clear('vols');

    for i = 1:size(J, 1)
        xdmf.write(dofmap, J(i, :), i);
    end
    clear('xdmf');

    xdmf = io.XDMF(sprintf('checkerboard-sensitivity-%s-total.xdmf', tag));
    xdmf.write(dofmap, vecnorm(J, 1,   1), 1);
    xdmf.write(dofmap, vecnorm(J, 2,   1), 2);
    xdmf.write(dofmap, vecnorm(J, inf, 1), 3);
    clear('xdmf');
end
