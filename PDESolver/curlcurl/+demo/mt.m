%% Curl-curl MT problem

%%
% Finite elements
elements = [
    fe.create_nedelec_element(2, 1), ...
    fe.create_nedelec_element(2, 2), ...
];

%%
% Manufactured solution (plane wave)
omega2 = +25i;
omega = -sqrt(omega2);
u = @(x) [-100i*omega*exp(1i*omega*x(2)), 0];
boundary = @(x) true;

%%
% Mesh resolutions for order 1 and 2
N{1} = [1 2 4 8 16 32 64];
N{2} = [1 2 4 8 16 32];

%%
% Function for assembling (curl curl - omega2)
assemble_mat = @(dofmap) assembling.assemble_curl_curl(dofmap, omega2, 1);
assemble_rhs = @(dofmap) zeros(dofmap.dim, 1, 'double');
degree_rise = 2;  % Rise the quadrature degree to account for non-polynomiality of u
assemble_err = @(dofmap, x) assembling.assemble_error_hcurl(dofmap, x, u, degree_rise);

%%
% Solve and evaluate errors on sequence of meshes using
% function <check_convergence.html |demo.common.check_convergence|>
demo.common.check_convergence(elements, assemble_mat, assemble_rhs, assemble_err, N, u, boundary);
