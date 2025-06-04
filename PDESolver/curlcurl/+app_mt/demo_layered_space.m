% Compute 2D
for n = [1, 2, 4, 8, 16, 32]
    mesh_resolution = [n, n];
    fe_order = 1;
    tic(); err = run_layered_space_test(mesh_resolution, fe_order);
    fprintf('Solved MT layered space on cube mesh %s in  %f sec with L^2 error %f\n', ...
            mat2str(mesh_resolution), toc(), err);
end

% Compute 3D
for n = [1, 2, 4, 8, 16]
    mesh_resolution = [n, n, n];
    fe_order = 1;
    tic(); err = run_layered_space_test(mesh_resolution, fe_order);
    fprintf('Solved MT layered space on cube mesh %s in  %f sec with L^2 error %f\n', ...
            mat2str(mesh_resolution), toc(), err);
end


function err = run_layered_space_test(mesh_resolution, fe_order)

    % Create geometry/mesh for cube (-1000, 1000)^dim
    mesh = meshing.generate_unit_cube_mesh(mesh_resolution);
    shift = 0.5*ones(mesh.dim, 1);
    scale = 2e3;
    mesh.set_coordinates(scale*(mesh.vertex_coords - shift));

    % Build Nedelec space
    element = fe.create_nedelec_element(mesh.dim, fe_order);
    for d = element.get_dof_entity_dims
        w = warning('off', 'Mesh:ExtraConnStored');  % Mute warning
        mesh.compute_connectivity(mesh.dim, d);
        warning(w);  % Restore warning
    end
    dofmap = assembling.build_dofmap(mesh, element);

    % Problem parameters
    mu = 1.25663706212e-6;          % vacuum permeability [H/m]
    eps = 8.8541878128e-12;         % vacuum permittivity [F/m]
    freq = 1e3;                     % frequency [Hz]
    rho = [1e9, 1e2, 1e1, 1e2];     % resistivities in layers [Ohm.m]
    d = [-inf, 0, 100, 250, +inf];  % depths of layers [m]

    % Get exact solution of 1D problem
    E_exact = @(z) app_mt.layered_space_plane_wave(freq, mu, rho, d, z);
    zmin = min(mesh.vertex_coords(mesh.dim, :));
    zmax = max(mesh.vertex_coords(mesh.dim, :));
    zi = linspace(zmin, zmax, 401);
    E = E_exact(zi);

    % Plot exact solution
    figure();
    if ~verLessThan('matlab', '9.5')
        sgtitle('Exact solution of 1D layered space');
    end
    subplot(1, 4, 1);
    plot(real(E), zi, '-');
    xlabel('Re E');
    ylabel('z');
    subplot(1, 4, 2);
    plot(imag(E), zi, '-');
    xlabel('Im E');
    ylabel('z');
    subplot(1, 4, 3);
    plot(abs(E), zi, '-');
    xlabel('|E|');
    ylabel('z');
    subplot(1, 4, 4);
    plot(angle(E), zi, '-');
    xlabel('phase E');
    ylabel('z');
    xlim([-pi, pi]);
    xticks([-pi, 0, pi]);
    xticklabels({'-\pi', '0', '\pi'});

    % Construct cell markers cm from depths d
    cm = zeros(mesh.num_entities(mesh.dim), 1);
    coords = mesh.get_cell_centroids();
    z = coords(mesh.dim, :);
    for layer = 2:numel(d)
        cm(z >= d(layer-1) & z < d(layer)) = layer-1;
    end

    % Plot cell markers and resistivity distribution
    figure();
    subplot(1, 2, 1);
    scatter23(coords, [], cm, 'filled');
    title('cell markers');
    colorbar();
    subplot(1, 2, 2);
    scatter23(coords, [], rho(cm), 'filled');
    title('resistivity');
    if ~verLessThan('matlab', '9.4')
        set(gca(), 'colorscale', 'log');
    end
    colorbar();
    caxis([min(rho), max(rho)]);

    % Compute cell-wise PDE coefficient lambda from cell markers
    omega = 2*pi*freq;
    lambda = 1j*omega*(1./rho(cm)) - omega^2*eps;

    % Assemble linear system
    A = assembling.assemble_curl_curl(dofmap, -lambda, mu);
    b = zeros(dofmap.dim, 1);

    % Apply BCs
    w = warning('off', 'Mesh:ExtraConnStored');  % Mute warning
    mesh.compute_connectivity(mesh.dim, mesh.dim-1);
    warning(w);  % Restore warning
    mesh.compute_boundary_facets();
    switch mesh.dim
    case 2
        exact_sol = @(x) [E_exact(x(mesh.dim)), 0];
    case 3
        exact_sol = @(x) [E_exact(x(mesh.dim)), 0, 0];
    end
    [bc_dofs, bc_vals] = assembling.build_dirichlet_dofs(dofmap, {@(x) true}, {exact_sol});
    mesh.clear_boundary_facets();
    [A, b] = assembling.apply_dirichlet_bc(A, b, bc_dofs, bc_vals);

    % Solve
    x = A\b;
    E = assembling.create_fe_eval_function(dofmap, x);  %#ok<NASGU>

    % Output to file
    f = io.XDMF(sprintf('mt_layered_dim%d.xdmf', mesh.dim));
    f.write(dofmap, real(x), 1);
    f.write(dofmap, imag(x), 2);
    f.write(dofmap, abs(x), 3);    % FIXME: This might be wrong
    f.write(dofmap, angle(x), 4);  % FIXME: This might be wrong
    f.flush();

    % Compute L^2 error
    function val = exact_sol_wrapper(x)
        if isa(x, 'sym')
            % Bypass symbolic computation of curl, which does not work
            val = zeros(1, mesh.dim);
        else
            val = exact_sol(x);
        end
    end
    err = assembling.assemble_error_hcurl(dofmap, x, @exact_sol_wrapper, 0);
end


function varargout = scatter23(coords, varargin)
    switch size(coords, 1)
    case 2
        [varargout{1:nargout}] = scatter(coords(1, :), coords(2, :), varargin{:});
    case 3
        [varargout{1:nargout}] = scatter3(coords(1, :), coords(2, :), coords(3, :), varargin{:});
    end
end
