% Create mesh
N = 128;
mesh = meshing.generate_half_disc_mesh(N);
w = warning('off', 'Mesh:ExtraConnStored');  % Mute warning
mesh.compute_connectivity(mesh.dim, mesh.dim-1);
warning(w);  % Restore warning

% Create reference finite elements
element_l2 = fe.create_p0_element(mesh.dim);
element_hdiv = fe.create_raviart_thomas_element(mesh.dim, 1);

% Build dofmaps
dofmap_l2 = assembling.build_dofmap(mesh, element_l2);
dofmap_hdiv = assembling.build_dofmap(mesh, element_hdiv);

% Assemble matrices
[M, D] = assembling.assemble_hdiv_operators(dofmap_hdiv, dofmap_l2, @(x, c) 1, 0);

% Build facet markers with boundary conditions
tol = 1e-8;
where_neumann = @(x, b) b && x(2) < tol;
where_dirichlet = @(x, b) b && norm(x) > 1-tol;
mesh.compute_boundary_facets();
markers = meshing.mark_entities(mesh, mesh.dim-1, [], where_neumann, 1);
markers = meshing.mark_entities(mesh, mesh.dim-1, markers, where_dirichlet, 2);

% Prepare right-hand side vectors (with natural boundary condition)
rhs_hdiv = assembling.assemble_neumann_source(dofmap_hdiv, markers, 2, ...
                                              @(x, n, c) x(2), 1);
rhs_l2 = zeros(dofmap_l2.dim, 1);

% Compute essential boundary condition
[bc_dofs, bc_values] = assembling.build_dirichlet_dofs(dofmap_hdiv, ...
    markers, {1}, {@(x) [0, +1]});

% Apply essential boundary conditions while preserving symmetry
[M, rhs_hdiv] = assembling.apply_dirichlet_bc(M, rhs_hdiv, bc_dofs, bc_values);
rhs_l2 = rhs_l2 - D(:, bc_dofs)*bc_values;
D(:, bc_dofs) = 0;

% Solve directly
[~, x_l2_direct] = solve_direct(M, D, rhs_hdiv, rhs_l2);

figure;
meshing.plot_mesh(mesh, 'cell_markers', x_l2_direct);

% Form initial guess fulfilling BCs
x0_hdiv = zeros(dofmap_hdiv.dim, 1);
x0_l2 = zeros(dofmap_l2.dim, 1);
[~, x0_hdiv] = assembling.apply_dirichlet_bc([], x0_hdiv, bc_dofs, bc_values, ...
                                             'symmetric', false);

% Solve iteratively
[~, x_l2_iterative] = solve_iterative(M, D, rhs_hdiv, rhs_l2, x0_hdiv, x0_l2);

figure;
meshing.plot_mesh(mesh, 'cell_markers', x_l2_iterative);

fprintf('norm(direct-iterative) = %e\n', norm(x_l2_direct-x_l2_iterative));


function [x_hdiv, x_l2] = solve_direct(M, D, rhs_hdiv, rhs_l2)
    fprintf('Solving directly...');
    tic;

    % Form system matrix
    Z = sparse([], [], [], size(D, 1), size(D, 1));
    A = [M, D.'; D, Z];

    % Form right-hand side
    b = [rhs_hdiv; rhs_l2];

    % Solve directly
    x = A\b;

    % Split the solution into blocks
    x_hdiv = x(1:size(M, 1), 1);
    x_l2 = x(size(M, 1)+1:end, 1);

    fprintf('%f sec\n', toc);
end


function [x_hdiv, x_l2] = solve_iterative(M, D, rhs_hdiv, rhs_l2, x0_hdiv, x0_l2)
    fprintf('Solving iteratively...');
    tic;

    % Return gracefully if HSL_MI20 is not installed
    if ~solving.HSLMI20.is_installed()
        fprintf('HSL_MI20 not installed; returning zero\n');
        x_hdiv = x0_hdiv;
        x_l2 = x0_l2;
        return;
    end

    [N2, N1] = size(D);

    % Inverse diagonal of M
    Mdiag_inv = 1./diag(M);
    assert(~issparse(Mdiag_inv));

    % Preconditioner for Schur: AMG with D*inv(diag(M))*DT
    Mdiag_inv_mat = spdiags(Mdiag_inv, 0, N1, N1);
    assert(issparse(Mdiag_inv_mat));
    Sdiag = D*Mdiag_inv_mat*D.';
    assert(issparse(Sdiag));
    clear Mdiag_inv_mat;
    amg_S = solving.HSLMI20(Sdiag);

    % Preconditioner for Hdiv block: scaling by inverse diagonal of mass
    % TODO: Try Chebyshev+Jacobi with full M?
    function x = precondition_mass(b)
        x(:) = Mdiag_inv.*b;
    end

    % Block-diagonal preconditioner for the system
    tmp = zeros(N1+N2, 1);
    function x = precondition_system(b)
        tmp(N1+1:N1+N2) = amg_S.solve(b(N1+1:N1+N2));
        tmp(1:N1) = precondition_mass(b(1:N1));
        x = tmp;
    end

    % System operator
    tmp_y = zeros(N1+N2, 1);
    function y = Afun(z)
        tmp_y(1:N1, 1) = D.'*z(N1+1:N1+N2) + M*z(1:N1);
        tmp_y(N1+1:N1+N2, 1) = D*z(1:N1);
        y = tmp_y;
    end

    % Form system right-hand side and initial guess
    rhs = [rhs_hdiv; rhs_l2];
    x0 = [x0_hdiv; x0_l2];

    % Run preconditioned MINRES
    x = minres(@Afun, rhs, 1e-6, 100, @precondition_system, [], x0);

    % Split the solution into blocks
    x_hdiv = x(1:size(M, 1), 1);
    x_l2 = x(size(M, 1)+1:end, 1);

    fprintf('%f sec\n', toc);
end
