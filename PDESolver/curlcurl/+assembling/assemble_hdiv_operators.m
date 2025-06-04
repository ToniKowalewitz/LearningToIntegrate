function [M, D] = assemble_hdiv_operators(dofmap_hdiv, dofmap_l2, sigma, degree_rise)
    % Assemble mass matrix and divergence on H(div)
    %
    %   M_{ij} = \int \phi_i \phi_j / sigma
    %   D_{ij} = \int \psi_i \div\phi_j
    %
    %   \phi_k - H(div)-conforming basis functions (Raviart-Thomas)
    %   \psi_k - L^2-conforming basis functions (Discontinuous Lagrange)
    %
    % INPUT PARAMETERS
    %   dofmap_hdiv ... Dofmap of H(div) space
    %   dofmap_l2   ... Dofmap of L^2 space
    %   sigma       ... Coefficient in the mass matrix,
    %                   function handle of signature:
    %
    %                     val = sigma(x, c);
    %
    %                   where:
    %
    %                     x   ... physical coordinate (row vector),
    %                     c   ... cell index,
    %                     val ... scalar value.
    %
    %  degree_rise  ... Quadrature degree contribution from 1/sigma
    %
    % OUTPUT PARAMETERS
    %
    %  M ... Mass matrix
    %  D ... Divergence matrix

    assert(strcmp(dofmap_hdiv.element.mapping, 'contravariant'));
    assert(strcmp(dofmap_l2.element.mapping,   'affine'));

    % Extract common mesh
    assert(dofmap_hdiv.mesh == dofmap_l2.mesh);
    mesh = dofmap_hdiv.mesh;
    dim = mesh.dim;
    num_cells = size(mesh.cells, 2);

    % Pick quadrature rule
    quad_degree = max([
        2*dofmap_hdiv.element.order + degree_rise;
        dofmap_hdiv.element.order - 1 + dofmap_l2.element.order;
    ]);
    [x, w] = fe.make_quadrature_rule(dim, quad_degree);
    num_quad_points = size(x, 1);

    % Tabulate basis and basis curls at quadrature points
    local_element_dim_hdiv = dofmap_hdiv.element.fe_space_dim;
    local_element_dim_l2   = dofmap_l2.element.fe_space_dim;
    basis_hdiv_divs = zeros(local_element_dim_hdiv, 1,   num_quad_points);
    basis_hdiv      = zeros(local_element_dim_hdiv, dim, num_quad_points);
    basis_l2        = zeros(local_element_dim_l2,   1,   num_quad_points);
    for k = 1:num_quad_points
        basis_hdiv_divs(:, :, k) = dofmap_hdiv.element.tabulate_basis_div(x(k, :));
        basis_hdiv     (:, :, k) = dofmap_hdiv.element.tabulate_basis    (x(k, :));
        basis_l2       (:, :, k) = dofmap_l2.element.tabulate_basis      (x(k, :));
    end

    % Fetch some data and preallocate temporaries
    cells = mesh.cells;
    coords = mesh.vertex_coords;
    cell_dofs_hdiv = dofmap_hdiv.cell_dofs;
    cell_dofs_l2   = dofmap_l2.cell_dofs;
    A_m = zeros(local_element_dim_hdiv, local_element_dim_hdiv);
    A_d = zeros(local_element_dim_l2,   local_element_dim_hdiv);
    jac = zeros(dim, dim);
    temp = zeros(local_element_dim_hdiv, dim);
    detJ = zeros(1, 1, 'double');
    detJinv = zeros(1, 1, 'double');
    offset = zeros(dim, 1);

    % Preallocate assembly data
    nnz_m = num_cells*local_element_dim_hdiv^2;
    nnz_d = num_cells*local_element_dim_hdiv*local_element_dim_l2;
    I_m = zeros(nnz_m, 1, 'double');
    I_d = zeros(nnz_d, 1, 'double');
    J_m = zeros(nnz_m, 1, 'double');
    J_d = zeros(nnz_d, 1, 'double');
    V_m = zeros(nnz_m, 1, 'double');
    V_d = zeros(nnz_d, 1, 'double');
    offsets_m = uint32(1:local_element_dim_hdiv^2);
    offsets_d = uint32(1:local_element_dim_hdiv*local_element_dim_l2);

    % Loop over cells
    for c = 1:num_cells

        % Compute geometric quantities
        offset(:) = coords(:, cells(dim+1, c));
        jac(:, :) = coords(:, cells(1:dim, c)) - coords(:, cells(dim+1, c));
        detJ(1) = abs(det(jac));
        detJinv(1) = 1/det(jac);

        % Zero cell matrices
        A_m(:, :) = 0;
        A_d(:, :) = 0;

        % Pull back coefficient to reference element
        assert(iscolumn(offset) && numel(offset) == dim);
        sigma_hat = @(xhat, c) sigma(jac*xhat.' + offset, c);

        % Loop over quadrature points and assembly cell integrals
        for k = 1:num_quad_points
            temp(:, :) = basis_hdiv(:, :, k)*jac.'*detJinv;
            A_m(:, :) = A_m + detJ*w(k)*(temp*temp.') / sigma_hat(x(k, :), c);
            A_d(:, :) = A_d + detJ*w(k)*basis_l2(:, 1, k)*basis_hdiv_divs(:, 1, k).'*detJinv;
        end

        % Compute global dof indices and store cell matrices
        [J_m(offsets_m), I_m(offsets_m)] = meshgrid(cell_dofs_hdiv(:, c));
        [J_d(offsets_d), I_d(offsets_d)] = meshgrid(cell_dofs_hdiv(:, c), ...
                                                    cell_dofs_l2(:, c));
        V_m(offsets_m) = A_m;
        V_d(offsets_d) = A_d;
        offsets_m(:) = offsets_m + numel(offsets_m);
        offsets_d(:) = offsets_d + numel(offsets_d);

    end

    % Assemble sparse matrices
    M = sparse(I_m, J_m, V_m, dofmap_hdiv.dim, dofmap_hdiv.dim);
    D = sparse(I_d, J_d, V_d, dofmap_l2.dim,   dofmap_hdiv.dim);

end
