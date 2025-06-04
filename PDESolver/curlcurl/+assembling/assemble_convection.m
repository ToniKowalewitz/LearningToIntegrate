function C = assemble_convection(dofmap, b, degree_rise)
    % Assemble convective term
    %
    %     \int b \cdot \nabla\phi_j \phi_i
    %
    % where \phi_k are (affine-mapped) basis functions,
    % b is a coefficient (mapped covariantly to reference) of signature
    %
    %     function val = b(x, c)
    %
    % where
    %
    %     x   ... physical coordinate, row vector, shape [1, dim]
    %     c   ... cell index
    %     val ... wind, row vector, shape [1, dim]

    assert(strcmp(dofmap.element.mapping, 'affine'));

    % Extract common mesh
    mesh = dofmap.mesh;
    dim = mesh.dim;
    num_cells = size(mesh.cells, 2);

    % Pick quadrature rule
    quad_degree = 2*dofmap.element.order - 1 + degree_rise;
    [x, w] = fe.make_quadrature_rule(dim, quad_degree);
    num_quad_points = size(x, 1);

    % Tabulate basis and basis curls at quadrature points
    local_element_dim   = dofmap.element.fe_space_dim;
    basis        = zeros(local_element_dim, 1,   num_quad_points);
    basis_grad   = zeros(local_element_dim, dim, num_quad_points);
    for k = 1:num_quad_points
        basis     (:, :, k) = dofmap.element.tabulate_basis     (x(k, :));
        basis_grad(:, :, k) = dofmap.element.tabulate_basis_grad(x(k, :));
    end

    % Fetch some data and preallocate temporaries
    cells = mesh.cells;
    coords = mesh.vertex_coords;
    cell_dofs = dofmap.cell_dofs;
    A_c = zeros(local_element_dim, local_element_dim);
    jac = zeros(dim, dim);
    detJ = zeros(1, 1, 'double');
    offset = zeros(dim, 1);

    % Preallocate assembly data
    nnz = num_cells*local_element_dim^2;
    I = zeros(nnz, 1, 'double');
    J = zeros(nnz, 1, 'double');
    V = zeros(nnz, 1, 'double');
    offsets = uint32(1:local_element_dim^2);

    % Loop over cells
    for c = 1:num_cells

        % Compute geometric quantities
        offset(:) = coords(:, cells(dim+1, c));
        jac(:, :) = coords(:, cells(1:dim, c)) - coords(:, cells(dim+1, c));
        detJ(1) = abs(det(jac));
        jac_inv(:, :) = inv(jac);

        % Zero cell matrices
        A_c(:, :) = 0;

        % Pull back coefficient to reference element
        assert(iscolumn(offset) && numel(offset) == dim);
        b_hat = @(xhat, c) b(jac*xhat.' + offset, c);

        % Loop over quadrature points and assembly cell integrals
        %num_quad_points = 1
        for k = 1:num_quad_points
            A_c(:, :) = A_c + detJ*w(k)*basis(:, 1, k)*b_hat(x(k, :), c)*jac_inv(:, :).'*basis_grad(:, :, k).';
        end

        % Compute global dof indices and store cell matrices
        [J(offsets), I(offsets)] = meshgrid(cell_dofs(:, c));

        V(offsets) = A_c;
        offsets(:) = offsets + numel(offsets);

    end

    % Assemble sparse matrices
    C = sparse(I, J, V, dofmap.dim,   dofmap.dim);

end
