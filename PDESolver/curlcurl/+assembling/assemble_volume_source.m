function b = assemble_volume_source(dofmap, f, g, degree_rise)
    % Assemble volume source on L^2
    %
    %   b_{i} = \int f \phi_i + g \cdot \nabla\phi_i
    %
    %   \phi_k - L^2-conforming basis functions (affine mapped)
    %
    % INPUT PARAMETERS
    %   dofmap   ... Dofmap of L^2 space
    %   f        ... Coefficient, function handle of signature:
    %
    %                     val = f(x, c);
    %
    %                   where:
    %
    %                     x   ... physical coordinate (row vector),
    %                     c   ... cell index,
    %                     val ... scalar value.
    %
    %   g        ... Coefficient, function handle of signature:
    %
    %                     val = g(x, c);
    %
    %                   where:
    %
    %                     x   ... physical coordinate (row vector of shape [1, dim]),
    %                     c   ... cell index,
    %                     val ... vector value (row vector of shape [1, dim]).
    %
    % NOTE
    %  One should pass in [] as f and/or g if the respective
    %  term is not required.
    %
    % OUTPUT PARAMETERS
    %  b ... RHS resulting from source term

    assert(strcmp(dofmap.element.mapping, 'affine'));

    assemble_f = ~isempty(f);
    assemble_g = ~isempty(g);

    % Fetch some data
    mesh = dofmap.mesh;
    dim = mesh.dim;

    % Pick quadrature rule
    quad_degree = dofmap.element.order + degree_rise;
    if ~assemble_f
        quad_degree = max(quad_degree-1, 0);
    end
    [x, w] = fe.make_quadrature_rule(dim, quad_degree);
    num_quad_points = size(x, 1);

    % Tabulate basis and basis grads at quadrature points
    local_element_dim = dofmap.element.fe_space_dim;
    basis      = zeros(local_element_dim,   1,   num_quad_points);
    basis_grad = zeros(local_element_dim, dim,   num_quad_points);
    for k = 1:num_quad_points
        basis(:, :, k) = dofmap.element.tabulate_basis(x(k, :));
        if assemble_g
            basis_grad(:, :, k) = dofmap.element.tabulate_basis_grad(x(k, :));
        end
    end

    % Fetch some data and preallocate temporaries
    num_cells = size(mesh.cells, 2);
    cells = mesh.cells;
    coords = mesh.vertex_coords;
    cell_dofs = dofmap.cell_dofs;
    jac = zeros(dim, dim);
    offset = zeros(dim, 1);
    b_e = zeros(local_element_dim, 1);

    % Allocate return value
    b = zeros(dofmap.dim, 1);

    % Loop over cells
    for c = 1:num_cells

        % Compute geometric quantities
        offset(:) = coords(:, cells(dim+1, c));
        jac(:, :) = coords(:, cells(1:dim, c)) - coords(:, cells(dim+1, c));
        detJ(1) = abs(det(jac));
        if assemble_g
            jac_inv(:, :) = inv(jac);
        end

        % Zero cell matrices
        b_e(:) = 0;

        % Pullback to reference element
        assert(iscolumn(offset) && numel(offset) == dim);
        f_hat = @(xhat, c) f(jac*xhat.' + offset, c);
        g_hat = @(xhat, c) g(jac*xhat.' + offset, c);

        % Loop over quadrature points and assembly cell integrals
        for k = 1:num_quad_points
            if assemble_f
                b_e(:) = b_e + basis(:, 1, k)*f_hat(x(k, :), c)*detJ*w(k);
            end
            if assemble_g
                b_e(:) = b_e + basis_grad(:, :, k)*jac_inv*g_hat(x(k, :), c).'*detJ*w(k);
            end
        end

        b(cell_dofs(:, c)) = b(cell_dofs(:, c)) + b_e;

    end

end
