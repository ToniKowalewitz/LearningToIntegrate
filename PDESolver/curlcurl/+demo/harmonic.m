%% Harmonic function

%%
% Finite elements
elements = [
    fe.create_lagrange_element(2, 1), ...
    fe.create_lagrange_element(2, 2), ...
    fe.create_lagrange_element(2, 3), ...
    fe.create_lagrange_element(2, 4), ...
    fe.create_lagrange_element(2, 5), ...
];

%%
% Manufactured solution
omega2 = 0;
u = @(x) x(1)^2 - x(2)^2;

%%
% Define Dirichlet boundary on x=1 and y=1
tol = 1e-8;
boundary = @(x) x(1) > 1-tol || x(2) > 1-tol;

%%
% Mesh resolutions for orders one to five
N{1} = [2 4 8 16 32 64];
N{2} = [2 4 8 16 32];
N{3} = [2 4 8 16];
N{4} = [2 4 8];
N{5} = [2 4];

%%
% Function for assembling (curl curl - omega2) and errors
assemble_mat = @(dofmap) assembling.assemble_laplace(dofmap, 1, omega2);
assemble_rhs = @(dofmap) zeros(dofmap.dim, 1, 'double');
degree_rise = 0;  % Rise the quadrature degree to account for non-polynomiality of u
assemble_err = @(dofmap, x) assembling.assemble_error_h1(dofmap, x, u, degree_rise);

%%
% Solve and evaluate errors on sequence of meshes using
% function <check_convergence.html |demo.common.check_convergence|>
demo.common.check_convergence(elements, assemble_mat, assemble_rhs, assemble_err, N, u, boundary);
