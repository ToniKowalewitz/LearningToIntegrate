function b = assemble_laplace_sensitivity(dofmap, U, M1, M2, lambda)
    % Assemble derivative of diffusion operator w.r.t. diffusivity
    %
    % Assemble dense matrix
    %
    %     b_kj = \int \psi_j \nabla u_k \cdot \nabla v_k ...
    %          + \lambda \int \psi_j u_k v_k
    %
    %                        (no summation over k or any index)
    %
    % where u and v are given by
    %
    %     u_k = U_im M1_mk \phi_i,    (sum over i, m),
    %     v_k = U_im M2_mk \phi_i,    (sum over i, m),
    %
    % lambda is a given (optional) scalar,
    % \phi_i are basis functions from H^1-conforming space
    % given by dofmap, and \psi_j are P_0 basis functions.

    assert(strcmp(dofmap.element.mapping, 'affine'));
    assert(size(U, 1) == dofmap.dim);
    assert(all(size(M1) == size(M2)), 'Expected same shapes of M1 and M2');
    assert(size(U, 2) == size(M1, 1), 'Expected as many columns of U as rows in M1');
    assert(nargin < 5 || isscalar(lambda));

    local_element_dim = dofmap.element.fe_space_dim;
    dim = dofmap.mesh.dim;
    num_cells = size(dofmap.mesh.cells, 2);

    % Precompute mass term
    if nargin >= 5 && lambda ~= 0
        mass_term = lambda*precompute_mass_integral(dofmap.element);
    else
        mass_term = 0;
    end

    % Pick quadrature rule
    quad_degree = 2*(dofmap.element.order - 1);
    [x, w] = fe.make_quadrature_rule(dim, quad_degree);
    num_quad_points = size(x, 1);

    % Tabulate basis gradients at quadrature points
    basis_grads = zeros(local_element_dim, dim, num_quad_points);
    for k = 1:num_quad_points
        basis_grads(:, :, k) = dofmap.element.tabulate_basis_grad(x(k, :));
    end

    % Fetch some data and preallocate temporaries
    cells = dofmap.mesh.cells;
    coords = dofmap.mesh.vertex_coords;
    cell_dofs = dofmap.cell_dofs;
    jac = zeros(dim, dim);
    detJ = zeros(1, 1, 'double');
    temp = zeros(local_element_dim, dim);
    A = zeros(local_element_dim, local_element_dim);
    U_cell = zeros(local_element_dim, size(U, 2), 'double');

    % Preallocate result
    b = zeros(size(M1, 2), num_cells);

    % Loop over cells
    for c = 1:num_cells

        % Compute geometric quantities
        jac(:, :) = coords(:, cells(1:dim, c)) - coords(:, cells(dim+1, c));
        detJ(1) = abs(det(jac));

        % Take precomputed mass term
        A(:) = mass_term*detJ;

        % Loop over quadrature points and assemble cell integral
        for k = 1:num_quad_points
            temp(:, :) = basis_grads(:, :, k)/jac;
            A = A + detJ*w(k)*(temp*temp.');
        end

        % Get cell coefficients
        U_cell(:) = U(cell_dofs(:, c), :);

        % Store to global vector
        b(:, c) = b(:, c) + dot((A*U_cell)*M1, U_cell*M2, 1).';

    end

end


function M = precompute_mass_integral(element)
    [x, w] = fe.make_quadrature_rule(element.simplex.dim, 2*element.order);
    M = zeros(element.fe_space_dim, element.fe_space_dim);
    for k = 1:size(x, 1)
        basis = element.tabulate_basis(x(k, :));
        M = M + w(k)*(basis*basis.');
    end
end
