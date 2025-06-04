% This program solves the Poisson problem in the following form:
% u - \grad(p) = 0   in \Omega
%      \div(u) = f   in \Omega
% with boundary conditions:
%            p = u_0 on \Gamma_D
%     u\cdot n = g   on \Gamma_N
% and
% \Omega = [0,1]X[0,1]
% \Gamma_D = left and right edges
% \Gamma_N = bottom and top edges
% u_0 = 0
% g = sin(5*x(1))
% f = -10\exp(-((x - 0.5)^2 + (y - 0.5)^2) / 0.02)
% Reference: https://fenicsproject.org/docs/dolfin/2019.1.0/python/demos/mixed-poisson/demo_mixed-poisson.py.html
% Note the sign of f is flipped compared to the FEniCS demo.

% Create mesh
N = 128;
mesh = meshing.generate_unit_cube_mesh([N, N]);
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
mesh.compute_boundary_facets();
% left: 1, right: 2, bottom: 3, top: 4
markers = meshing.mark_entities(mesh, mesh.dim-1, [], @(x, c) x(1)<tol, 1);
markers = meshing.mark_entities(mesh, mesh.dim-1, markers, @(x, c) x(1)>1.0-tol, 2);
markers = meshing.mark_entities(mesh, mesh.dim-1, markers, @(x, c) x(2)<tol , 3);
markers = meshing.mark_entities(mesh, mesh.dim-1, markers, @(x, c) x(2)>1.0-tol , 4);
figure;
meshing.plot_mesh(mesh, 'facet_markers', markers);

% Assemble right-hand side
u_0 = @(x, n, c) 0;
rhs_hdiv = assembling.assemble_neumann_source(dofmap_hdiv, markers, 1, u_0, 0);
rhs_hdiv = rhs_hdiv ...
         + assembling.assemble_neumann_source(dofmap_hdiv, markers, 2, u_0, 0);
f = @(x, c) -10*exp(-((x(1) - 0.5)^2 + (x(2) - 0.5)^2) / 0.02);
rhs_l2 = assembling.assemble_volume_source(dofmap_l2, f, [], 0);

% Compute essential boundary condition
[bc_dofs, bc_values] = assembling.build_dirichlet_dofs(dofmap_hdiv, ...
    markers, {3, 4}, {@(x) [0, -sin(5*x(1))], @(x) [0, sin(5*x(1))]});

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
