function A = assemble_robin_matrix(A, dofmap, markers, marker_value, g, degree_rise)
    % Assemble Robin exterior surface integral
    %
    %   \int_{markers==marker_value} g(x, n) \phi_i \phi_j ds,
    %
    % where \phi_j are H^1-conforming basis functions.
    % Integration domain are the exterior boundary facets where
    % the supplied markers take the provided value.
    %
    % INPUT PARAMETER
    %   A            ... Sparse matrix to add the values to; be sure
    %                    it has the right sparsity pattern
    %   dofmap       ... Struct, containing mesh and FE element objects
    %                    (the latter contain H^1 conforming basis functions)
    %                    as well as the cell-2-DOF mapping.
    %   markers      ... Facet markers
    %   marker_value ... Facet marker value over which integrate
    %   g            ... Coefficient, a function of signature
    %
    %                      val = g(x, n, c);
    %
    %                    where
    %
    %                      x ... physical coordinates
    %                      n ... outer normal
    %                      c ... cell index
    %
    %   degree_rise  ... Contribution from g to quadrature degree
    %
    % OUTPUT PARAMETER
    %   A            ... Sparse matrix

    assert(strcmp(dofmap.element.family, 'Lagrange'));
    assert(isscalar(marker_value));  % NB: Add test if supporting multiple markers

    % Fetch data
    local_element_dim = dofmap.element.fe_space_dim;
    dim = dofmap.mesh.dim;
    cells = dofmap.mesh.cells;
    coords = dofmap.mesh.vertex_coords;
    cell_dofs = dofmap.cell_dofs;
    num_facets = dofmap.element.simplex.num_entities(dim-1);
    normals = dofmap.element.simplex.normals;
    cell_facets = dofmap.mesh.get_connectivity(dim, dim-1);
    boundary_facets = dofmap.mesh.get_boundary_facets();

    % Pick quadrature rule on (dim-1)-simplex
    quad_degree = 2*dofmap.element.order + degree_rise;  % FIXME: -2?
    [x, w] = fe.make_quadrature_rule(dim-1, quad_degree);
    num_quad_points = size(x, 1);
    area = 1/factorial(dim-1);

    % Transform quadrature rule to reference facets
    xf = zeros(num_quad_points, dim, num_facets);
    wf = zeros(num_quad_points, num_facets);
    basis = zeros(local_element_dim, num_quad_points, num_facets);
    for f = 1:num_facets
        transform = dofmap.element.simplex.get_facet_transform(f);
        for k = 1:num_quad_points
            xf(k, :, f) = transform(x(k, :).');
            wf(k, f) = w(k)*dofmap.element.simplex.facet_area(f)/area;
            basis(:, k, f) = dofmap.element.tabulate_basis(xf(k, :, f));
        end
    end

    % Preallocate temporaries
    offset = zeros(dim, 1);
    jac = zeros(dim, dim);
    jac_inv = zeros(dim, dim);
    temp = zeros(local_element_dim, local_element_dim);
    detJ = zeros(1, 1, 'double');
    n = zeros(dim, 1, 'double');

    % Loop over cells with boundary facets
    for c = find(any(boundary_facets))

        % Find marked local facets on boundary
        local_facets = full(boundary_facets(:, c) & ...
                            markers(cell_facets(:, c)) == marker_value);

        % If no marked continue to another cell
        if all(~local_facets)
            continue
        end

        % Compute geometric quantities
        offset(:) = coords(:, cells(dim+1, c));
        jac(:, :) = coords(:, cells(1:dim, c)) - coords(:, cells(dim+1, c));
        detJ(1) = abs(det(jac));
        jac_inv(:, :) = inv(jac);

        % Pull back coefficient to reference element
        assert(iscolumn(offset) && numel(offset) == dim);
        g_hat = @(xhat, n, c) g(jac*xhat.' + offset, n, c);

        % Zero cell vector
        temp(:) = 0;

        % Loop over local marked facets and quadrature points
        for f = find(local_facets).'
            n(:) = jac_inv.'*normals(:, f);
            S = norm(n);
            n(:) = n/S;

            for k = 1:num_quad_points
                % TODO: Should loop only through facet-supported basis functions?
                temp = temp + ...
                    wf(k, f)*g_hat(xf(k, :, f), n, c)*basis(:, k, f)*basis(:, k, f).'*S*detJ;
            end
        end

        % Add cell contribution to global matrix
        c_dofs(:) = cell_dofs(:, c);
        A(c_dofs, c_dofs) = A(c_dofs, c_dofs) + temp;
    end
end
