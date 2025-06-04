function [M, C] = assemble_hcurl_operators(dofmap)
    % Assemble mass and curl-curl matrix on H(curl)

    assert(strcmp(dofmap.element.mapping, 'covariant'));

    local_element_dim = dofmap.element.fe_space_dim;
    dim = dofmap.mesh.dim;
    num_cells = size(dofmap.mesh.cells, 2);

    % Pick quadrature rule
    quad_degree = 2*dofmap.element.order;
    [x, w] = fe.make_quadrature_rule(dim, quad_degree);
    num_quad_points = size(x, 1);

    % FIXME: Move to the element interface
    if dim == 2
        curl_shp = 1;
    else
        curl_shp = 3;
    end

    % Tabulate basis and basis curls at quadrature points
    basis_curls = zeros(local_element_dim, curl_shp, num_quad_points);
    basis       = zeros(local_element_dim, dim, num_quad_points);
    for k = 1:num_quad_points
        basis_curls(:, :, k) = dofmap.element.tabulate_basis_curl(x(k, :));
        basis      (:, :, k) = dofmap.element.tabulate_basis     (x(k, :));
    end

    % Fetch some data and preallocate temporaries
    cells = dofmap.mesh.cells;
    coords = dofmap.mesh.vertex_coords;
    cell_dofs = dofmap.cell_dofs;
    A_m = zeros(local_element_dim, local_element_dim);
    A_c = zeros(local_element_dim, local_element_dim);
    jac = zeros(dim, dim);
    jac_inv = zeros(dim, dim);
    temp = zeros(local_element_dim, dim);
    detJ = zeros(1, 1, 'double');
    detJinv = zeros(1, 1, 'double');

    % Preallocate assembly data
    nnz = num_cells*local_element_dim^2;
    I = zeros(nnz, 1, 'double');
    J = zeros(nnz, 1, 'double');
    V_m = zeros(nnz, 1, 'double');
    V_c = zeros(nnz, 1, 'double');
    offsets = uint32(1:local_element_dim^2);

    % Loop over cells
    for c = 1:num_cells

        % Compute geometric quantities
        jac(:, :) = coords(:, cells(1:dim, c)) - coords(:, cells(dim+1, c));
        detJ(1) = abs(det(jac));
        detJinv(1) = 1/detJ;
        jac_inv(:, :) = inv(jac);

        % Zero cell matrices
        A_m(:, :) = 0;
        A_c(:, :) = 0;

        % Loop over quadrature points and assembly cell integrals
        for k = 1:num_quad_points
            if dim == 2
                A_c(:, :) = A_c + detJinv*w(k)*basis_curls(:, :, k)*basis_curls(:, :, k).';
            else
                temp(:, :) = basis_curls(:, :, k)*jac.';
                A_c(:, :) = A_c + detJinv*w(k)*(temp*temp.');
            end

            temp(:, :) = basis(:, :, k)*jac_inv;
            A_m(:, :) = A_m + detJ*w(k)*(temp*temp.');
        end

        % Compute global dof indices and store cell matrix
        [J(offsets), I(offsets)] = meshgrid(cell_dofs(:, c));
        V_m(offsets) = A_m;
        V_c(offsets) = A_c;
        offsets(:) = offsets + local_element_dim^2;

    end

    % Assemble sparse matrix
    M = sparse(I, J, V_m, dofmap.dim, dofmap.dim);
    C = sparse(I, J, V_c, dofmap.dim, dofmap.dim);

end
