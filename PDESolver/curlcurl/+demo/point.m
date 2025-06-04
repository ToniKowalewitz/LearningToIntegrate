%% Laplace point source problem

%%
% Finite elements
elements = [
    fe.create_lagrange_element(2, 1), ...
    fe.create_lagrange_element(2, 2), ...
    fe.create_lagrange_element(2, 3), ...
];

%%
% Manufactured solution
omega2 = 0;
x0 = [0.5; 0];
u = @(x) -1/pi*log(norm(x(:) - x0(:)));

%%
% Define Dirichlet boundary away from point source
tol = 1e-8;
boundary = @(x) ( x(2) > tol ) || ( ( x(1) < tol ) || ( x(1) > 1-tol ) );

%%
% Mesh resolutions for orders one to five
N{1} = [2 4 8 16 32];
N{2} = [2 4 8 16];
N{3} = [2 4 8];

%%
% Function for assembling (Laplace - omega2) and errors
assemble_mat = @(dofmap) assembling.assemble_laplace(dofmap, 1, omega2);
assemble_rhs = @(dofmap) assembling.assemble_point_sources(dofmap, x0);
degree_rise = 2;  % Rise the quadrature degree to account for non-polynomiality of u
assemble_err = @(dofmap, x) assembling.assemble_error_h1(dofmap, x, u, degree_rise);

%%
% Solve and evaluate errors on sequence of meshes using
% function <check_convergence.html |demo.common.check_convergence|>
demo.common.check_convergence(elements, assemble_mat, assemble_rhs, assemble_err, N, u, boundary);


%%
% Finite elements
elements = [
    fe.create_lagrange_element(3, 1), ...
    fe.create_lagrange_element(3, 2), ...
    fe.create_lagrange_element(3, 3), ...
];

%%
% Manufactured solution
omega2 = 0;
x0 = [0.5; 0.5; 0];
u = @(x) 1/(2*pi*norm(x(:) - x0(:)));

%%
% Define Dirichlet boundary away from point source
tol = 1e-8;
boundary = @(x) ( x(3) > tol ) || ( ( ( x(1) < tol ) || ( x(1) > 1-tol ) ) && ( ( x(2) < tol ) || ( x(2) > 1-tol ) ) );

%%
% Mesh resolutions for orders one to five
N{1} = [2 4 8 16];
N{2} = [2 4 8];
N{3} = [2 4];

%%
% Function for assembling (curl curl - omega2) and errors
assemble_mat = @(dofmap) assembling.assemble_laplace(dofmap, 1, omega2);
assemble_rhs = @(dofmap) assembling.assemble_point_sources(dofmap, x0);
degree_rise = 2;  % Rise the quadrature degree to account for non-polynomiality of u
assemble_err = @(dofmap, x) assembling.assemble_error_h1(dofmap, x, u, degree_rise);

%%
% Solve and evaluate errors on sequence of meshes using
% function <check_convergence.html |demo.common.check_convergence|>
demo.common.check_convergence(elements, assemble_mat, assemble_rhs, assemble_err, N, u, boundary);
