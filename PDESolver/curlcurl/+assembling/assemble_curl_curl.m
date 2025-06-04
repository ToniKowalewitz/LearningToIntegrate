function M = assemble_curl_curl(dofmap, sigma, mue)
    % Assemble operator (\curl 1/mue \curl - \sigma)
    %
    %   curl-curl term:
    %   \int 1/mue \curl \phi_i \cdot \curl \phi_j
    %   mass matrix term:
    %       - \int \sigma \phi_i \phi_j
    %
    %   \phi_i       - H(curl) conforming basis functions
    %   \sigma, \mue - scalar or P_0 coefficients indexed by cell index
    %
    % INPUT PARAMETER
    %   dofmap ... Struct, containing mesh and FE element objects
    %              (the latter contain H^1 conforming basis functions)
    %              as well as the cell-2-DOF mapping.
    %   sigma  ... Vector of P_0 element values, indexed by cell index.
    %              Alternatively can be given as scalar (applying to all).
    %   mue    ... Vector of P_0 element values, indexed by cell index.
    %              Alternatively can be given as scalar (applying to all).
    %
    % OUTPUT PARAMETER
    %   M ... Sparse matrix [dofmap.dim x dofmap.dim], representing
    %         curl-curl term (and mass matrix, if sigma is given).

    % FIXME: Refactor this routine for repeated assembly (of mass)?!?

    assert(strcmp(dofmap.element.mapping, 'covariant'));

    local_element_dim = dofmap.element.fe_space_dim;
    dim = dofmap.mesh.dim;
    num_cells = size(dofmap.mesh.cells, 2);

    assert(isscalar(sigma) || (isvector(sigma) && length(sigma) == num_cells));
    assert(isscalar(mue) || (isvector(mue) && length(mue) == num_cells));

    % Determine which terms to assemble
    assemble_curlcurl = ~( isscalar(mue) && isinf(mue) );
    assemble_mass = ~( isscalar(sigma) && sigma == 0 );
    if ~assemble_curlcurl && ~assemble_mass
        M = sparse([], [], [], dofmap.dim, dofmap.dim);
        return
    end

    % Pick quadrature rule
    if assemble_mass
        quad_degree = 2*dofmap.element.order;
    else
        quad_degree = 2*(dofmap.element.order - 1);
    end
    [x, w] = fe.make_quadrature_rule(dim, quad_degree);
    num_quad_points = size(x, 1);

    % Wrap constant values as p0 function
    sigma_cell = wrap_as_p0(sigma, num_cells);
    mue_cell = wrap_as_p0(mue, num_cells);

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
    A = zeros(local_element_dim, local_element_dim);
    jac = zeros(dim, dim);
    jac_inv = zeros(dim, dim);
    temp = zeros(local_element_dim, dim);
    detJ = zeros(1, 1, 'double');
    detJinv = zeros(1, 1, 'double');

    % Preallocate assembly data
    nnz = num_cells*local_element_dim^2;
    I = zeros(nnz, 1, 'double');
    J = zeros(nnz, 1, 'double');
    V = zeros(nnz, 1, 'double');
    offsets = uint32(1:local_element_dim^2);

    % Loop over cells
    for c = 1:num_cells

        % Compute geometric quantities
        jac(:, :) = coords(:, cells(1:dim, c)) - coords(:, cells(dim+1, c));
        detJ(1) = abs(det(jac));
        detJinv(1) = 1/detJ;
        if assemble_mass
            jac_inv(:, :) = inv(jac);
        end

        % Zero cell matrix
        A(:, :) = 0;

        % Loop over quadrature points and assemble cell integrals
        for k = 1:num_quad_points
            if dim == 2 && assemble_curlcurl
                A(:, :) = A + detJinv*1/mue_cell(c)*w(k)*basis_curls(:, :, k)*basis_curls(:, :, k).';
            elseif assemble_curlcurl
                temp(:, :) = basis_curls(:, :, k)*jac.';
                A(:, :) = A + detJinv*1/mue_cell(c)*w(k)*(temp*temp.');
            end
            if assemble_mass
                temp(:, :) = basis(:, :, k)*jac_inv;
                A(:, :) = A - sigma_cell(c)*detJ*w(k)*(temp*temp.');
            end
        end

        % Compute global dof indices and store cell matrix
        [J(offsets), I(offsets)] = meshgrid(cell_dofs(:, c));
        V(offsets) = A;
        offsets(:) = offsets + local_element_dim^2;

    end

    % Assemble sparse matrix
    M = sparse(I, J, V, dofmap.dim, dofmap.dim);

end


function g = wrap_as_p0(f, num_cells)
    if isscalar(f)
        g = @(i) f;
    else
        assert(numel(f) == num_cells);
        g = f;
    end
end
