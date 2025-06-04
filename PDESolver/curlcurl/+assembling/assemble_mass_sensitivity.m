function b = assemble_mass_sensitivity(dofmap, U, V)
    % Assemble derivative of mass matrix w.r.t. coefficient
    %
    % Assemble dense matrix
    %
    %     b_kj = \int \psi_j u_m \cdot v_n,
    %
    %     with rhs indices m and n lexicographically packed
    %     into lhs index k,
    %
    % where u and v are (vector-valued) functions given by
    %
    %     u_m = U_im \phi_i,    (sum over i),
    %     v_n = V_in \phi_i,    (sum over i),
    %
    % \phi_i are (dofmap-valued) basis functions from the space
    % given by dofmap, and \psi_j are (scalar) P_0 basis functions.

    switch dofmap.element.mapping
    case 'affine'
        b = assemble_mass_sensitivity_h1(dofmap, U, V);
    case 'covariant'
        b = assemble_mass_sensitivity_hcurl(dofmap, U, V);
    otherwise
        error('not implemented!');
    end

end


function b = assemble_mass_sensitivity_h1(dofmap, U, V)
    % Assemble derivative of mass matrix on H^1 w.r.t. coefficient
    %
    % Assemble dense matrix
    %
    %     b_kj = \int \psi_j u_m v_n,
    %
    %     with rhs indices m and n lexicographically packed
    %     into lhs index k,
    %
    % where u and v are functions given by
    %
    %     u_m = U_im \phi_i,    (sum over i),
    %     v_n = V_in \phi_i,    (sum over i),
    %
    % \phi_i are basis functions from H^1-conforming space
    % given by dofmap, and \psi_j are (scalar) P_0 basis functions.

    assert(strcmp(dofmap.element.mapping, 'affine'));
    assert(size(U, 1) == dofmap.dim);
    assert(size(V, 1) == dofmap.dim);

    local_element_dim = dofmap.element.fe_space_dim;
    dim = dofmap.mesh.dim;
    num_cells = size(dofmap.mesh.cells, 2);

    % Tabulate basis at quadrature points
    [x, w] = fe.make_quadrature_rule(dofmap.element.simplex.dim, 2*dofmap.element.order);
    M = zeros(dofmap.element.fe_space_dim, dofmap.element.fe_space_dim);
    for k = 1:size(x, 1)
        basis = dofmap.element.tabulate_basis(x(k, :));
        M = M + w(k)*(basis*basis.');
    end

    % Fetch some data and preallocate temporaries
    cells = dofmap.mesh.cells;
    coords = dofmap.mesh.vertex_coords;
    cell_dofs = dofmap.cell_dofs;
    A = zeros(local_element_dim, local_element_dim);
    jac = zeros(dim, dim);
    detJ = zeros(1, 1, 'double');
    U_cell = zeros(local_element_dim, size(U, 2), 'double');
    V_cell = zeros(local_element_dim, size(V, 2), 'double');

    % Preallocate result
    b = zeros(size(U, 2)*size(V, 2), num_cells);

    % Loop over cells
    for c = 1:num_cells

        % Compute geometric quantities
        jac(:, :) = coords(:, cells(1:dim, c)) - coords(:, cells(dim+1, c));
        detJ(1) = abs(det(jac));

        % Take precomputed mass term
        A(:) = M*detJ;

        % Get cell coefficients
        U_cell(:) = U(cell_dofs(:, c), :);
        V_cell(:) = V(cell_dofs(:, c), :);

        % Store to global vector
        b(:, c) = b(:, c) + reshape(V_cell.'*(A*U_cell), [], 1);

    end

end


function b = assemble_mass_sensitivity_hcurl(dofmap, U, V)
    % Assemble derivative of mass matrix on H(curl) w.r.t. coefficient
    %
    % Assemble dense matrix
    %
    %     b_kj = \int \psi_j u_m \cdot v_n,
    %
    %     with rhs indices m and n lexicographically packed
    %     into lhs index k,
    %
    % where u and v are vector-valued functions given by
    %
    %     u_m = U_im \phi_i,    (sum over i),
    %     v_n = V_in \phi_i,    (sum over i),
    %
    % \phi_i are (vector-valued) basis functions from H(curl)-conforming space
    % given by dofmap, and \psi_j are (scalar) P_0 basis functions.

    assert(strcmp(dofmap.element.mapping, 'covariant'));
    assert(size(U, 1) == dofmap.dim);
    assert(size(V, 1) == dofmap.dim);

    local_element_dim = dofmap.element.fe_space_dim;
    dim = dofmap.mesh.dim;
    num_cells = size(dofmap.mesh.cells, 2);

    % Pick quadrature rule
    quad_degree = 2*dofmap.element.order;
    [x, w] = fe.make_quadrature_rule(dim, quad_degree);
    num_quad_points = size(x, 1);

    % Tabulate basis at quadrature points
    basis = zeros(local_element_dim, dim, num_quad_points);
    for k = 1:num_quad_points
        basis(:, :, k) = dofmap.element.tabulate_basis(x(k, :));
    end

    % Fetch some data and preallocate temporaries
    cells = dofmap.mesh.cells;
    coords = dofmap.mesh.vertex_coords;
    cell_dofs = dofmap.cell_dofs;
    A = zeros(local_element_dim, local_element_dim);
    jac = zeros(dim, dim);
    jac_inv = zeros(dim, dim);
    temp = zeros(local_element_dim, dim);
    detJ = zeros(1, 1, 'double');
    U_cell = zeros(local_element_dim, size(U, 2), 'double');
    V_cell = zeros(local_element_dim, size(V, 2), 'double');

    % Preallocate result
    b = zeros(size(U, 2)*size(V, 2), num_cells);

    % Loop over cells
    for c = 1:num_cells

        % Compute geometric quantities
        jac(:, :) = coords(:, cells(1:dim, c)) - coords(:, cells(dim+1, c));
        detJ(1) = abs(det(jac));
        jac_inv(:, :) = inv(jac);

        % Zero cell matrix
        A(:, :) = 0;

        % Loop over quadrature points and assemble cell integrals
        for k = 1:num_quad_points
            temp(:, :) = basis(:, :, k)*jac_inv;
            A(:, :) = A + detJ*w(k)*(temp*temp.');
        end

        % Get cell coefficients
        U_cell(:) = U(cell_dofs(:, c), :);
        V_cell(:) = V(cell_dofs(:, c), :);

        % Store to global vector
        b(:, c) = b(:, c) + reshape(V_cell.'*(A*U_cell), [], 1);

    end

end
