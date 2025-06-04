%% Laplace cavity problem

%%
% Finite elements
elements = [
    fe.create_lagrange_element(3, 1), ...
    fe.create_lagrange_element(3, 2), ...
    fe.create_lagrange_element(3, 3), ...
    fe.create_lagrange_element(3, 4), ...
    fe.create_lagrange_element(3, 5), ...
];

%%
% Manufactured solution (plane wave)
omega2 = 1.5^2 + 1^2;
u = @(x) cos(-1.5*x(1) + x(2));
boundary = @(x) true;

%%
% Mesh resolutions for order 1 and 2
N{1} = [1 2 4 8 16];
N{2} = [1 2 4 8];
N{3} = [1 2 4];
N{4} = [1 2 4];
N{5} = [1 2];

%%
% Function for assembling (curl curl - omega2) and errors
assemble_mat = @(dofmap) assembling.assemble_laplace(dofmap, 1, omega2);
assemble_rhs = @(dofmap) zeros(dofmap.dim, 1, 'double');
degree_rise = 2;  % Rise the quadrature degree to account for non-polynomiality of u
assemble_err = @(dofmap, x) assembling.assemble_error_h1(dofmap, x, u, degree_rise);

%%
% Solve and evaluate errors on sequence of meshes using
% function <check_convergence.html |demo.common.check_convergence|>
demo.common.check_convergence(elements, assemble_mat, assemble_rhs, assemble_err, N, u, boundary);
