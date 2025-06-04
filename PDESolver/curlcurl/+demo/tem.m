run_demo(2, 1);
run_demo(2, 2);
run_demo(3, 1);
run_demo(3, 2);


function run_demo(dim, order)

    % Build mesh
    switch dim
    case 2
        mesh = meshing.generate_unit_cube_mesh([32, 32]);
    case 3
        mesh = meshing.generate_unit_cube_mesh([8, 8, 8]);
    end

    % Build Nedelec function space
    element = fe.create_nedelec_element(mesh.dim, order);
    dofmap = build_dofmap(mesh, element);

    % Mark edges with source
    markers = build_markers(mesh);

    % Assemble linear system
    [A, M, x0] = assemble_system(dofmap, markers);

    % Output files
    filename = sprintf('tem_dim%d_order%d.xdmf', dim, order);
    f = io.XDMF(filename, 'xml');

    % Timestepping parameters
    num_timesteps = 16;
    dt = 0.001;
    t = 0;
    x = x0;

    % Store initial value as XDMF
    f.write(dofmap, x, t);

    % Compute Cholesky factor
    fprintf('dim %d order %d time-stepping (Cholesky factor) ', dim, order);
    tic();
    [R, flag, P] = chol(dt*A + M);
    assert(flag == 0);
    fprintf(' %f sec\n', toc());

    % Do time-stepping
    fprintf('dim %d order %d time-stepping (%d steps) ', dim, order, num_timesteps);
    tic();
    for i = 1:num_timesteps
        t = t + dt;

        % Solve
        x = P*(R\(R'\(P'*(M*x))));

        % Store current step as XDMF
        f.write(dofmap, x, t);

        fprintf('.');
    end
    fprintf(' %f sec\n', toc());

    % Solve by rational best approximation (only 3 poles)
    num_poles = 3;
    fprintf('dim %d order %d RBA %d poles ', dim, order, num_poles);
    tic();
    r = rba.create_rba_matrix(num_poles, A, M, x0, @(A, b) A\b);
    x_rba = r(t);
    f.write(dofmap, x_rba, 2*t);  % Dummy time
    fprintf(' %f sec\n', toc());

    % Compare the two solutions
    fprintf('rel. diff(time-stepping, RBA): %f percent\n', 100*norm(x-x_rba)/norm(x));

end


function markers = build_markers(mesh)

    w = warning('off', 'Mesh:ExtraConnStored');  % Mute warning
    mesh.compute_connectivity(mesh.dim-1, 0);
    mesh.compute_connectivity(mesh.dim, 1);
    warning(w);  % Restore warning
    mesh.clear_connectivity(1, 0);
    mesh.compute_boundary_facets();

    % Helper functions
    tol = 1e-6;
    near = @(x, a) abs(x-a) < tol;
    between = @(x, a, b) x > a - tol && x < b + tol;
    function markers = mark_(mesh, markers, where, value)
        % If in dimension 3 require further that z == 0
        if mesh.dim == 3
            where_ = where;
            where = @(x, b) b && near(x(3), 0) && where_(x, b);
        end
        markers = meshing.mark_entities(mesh, 1, markers, where, value);
    end

    % Mark
    source_edges_char_fun = @(x, b) near(x(1), 0.25) && between(x(2), 0.25, 0.75);
    markers = mark_(mesh, [], source_edges_char_fun, 1);
    source_edges_char_fun = @(x, b) near(x(1), 0.75) && between(x(2), 0.25, 0.75);
    markers = mark_(mesh, markers, source_edges_char_fun, 3);
    source_edges_char_fun = @(x, b) near(x(2), 0.25) && between(x(1), 0.25, 0.75);
    markers = mark_(mesh, markers, source_edges_char_fun, 2);
    source_edges_char_fun = @(x, b) near(x(2), 0.75) && between(x(1), 0.25, 0.75);
    markers = mark_(mesh, markers, source_edges_char_fun, 4);
end


function dofmap = build_dofmap(mesh, element)

    % Build additional mesh connectivities
    w = warning('off', 'Mesh:ExtraConnStored');  % Mute warning
    for d = element.get_dof_entity_dims()
        mesh.compute_connectivity(mesh.dim, d);
        mesh.clear_connectivity(d, 0);
    end
    mesh.compute_connectivity(mesh.dim, mesh.dim-1);
    mesh.clear_connectivity(mesh.dim-1, 0);
    warning(w);  % Restore warning

    % Build dofmap
    dofmap = assembling.build_dofmap(mesh, element);
end


function [A, M, x0] = assemble_system(dofmap, markers)

    mesh = dofmap.mesh;

    % Dirichlet BC
    if mesh.dim == 2
        % Homogeneous Dirichlet on whole boundary
        [bc_dofs, ~] = assembling.build_dirichlet_dofs(dofmap, ...
            {@(x) true}, {@(x) [0, 0]});
    elseif mesh.dim == 3
        % Homogeneous Dirichlet on z>0 boundary
        tol = 1e-6;
        [bc_dofs, ~] = assembling.build_dirichlet_dofs(dofmap, ...
            {@(x) x(3) > tol}, {@(x) [0, 0, 0]});
    end

    % Clear unneeded connectivities
    for d = 2:mesh.dim-1
        mesh.clear_connectivity(mesh.dim, d);
        mesh.clear_connectivity(d, 0);
    end
    mesh.clear_boundary_facets();

    % Assemble operators
    [M, A] = assembling.assemble_hcurl_operators(dofmap);

    % Assemble right-hand side
    b = zeros(dofmap.dim, 1);
    b = b + assembling.assemble_edge_source(dofmap, markers, 1, @(x, tau) -tau(2), 0);
    b = b + assembling.assemble_edge_source(dofmap, markers, 2, @(x, tau) +tau(1), 0);
    b = b + assembling.assemble_edge_source(dofmap, markers, 3, @(x, tau) +tau(2), 0);
    b = b + assembling.assemble_edge_source(dofmap, markers, 4, @(x, tau) -tau(1), 0);
    mesh.clear_connectivity(mesh.dim, 1);

    % Apply BCs
    [A, b] = assembling.apply_dirichlet_bc(A, b, bc_dofs, 0);
    M = assembling.apply_dirichlet_bc(M, [], bc_dofs, 0);

    % Get initial value
    x0 = M\b;

end
