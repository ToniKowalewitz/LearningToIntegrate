function b = assemble_edge_source(dofmap, markers, marker_value, g, degree_rise)
    % Assemble edge integral
    %
    %   \int_{markers==marker_value} g(x, \tau) \phi_j\cdot\tau ds,
    %
    % where
    %
    %   \phi_j          - H(curl) conforming basis functions
    %   \tau            - edge unit tangent (of unspecified direction)
    %
    % INPUT PARAMETER
    %   dofmap       ... Dofmap of H(curl) space
    %   markers      ... Edge markers
    %   marker_value ... Edge marker value over which integrate
    %   g            ... Scalar function of coordinates x and edge normal \tau;
    %                    assumed continuous over vertices
    %   degree_rise  ... Contribution from g to quadrature degree
    %
    % OUTPUT PARAMETER
    %   b            ... Vector

    assert(strcmp(dofmap.element.mapping, 'covariant'));
    assert(isscalar(marker_value));  % NB: Add test if supporting multiple markers

    % Fetch data
    local_element_dim = dofmap.element.fe_space_dim;
    dim = dofmap.mesh.dim;
    num_cells = size(dofmap.mesh.cells, 2);
    cells = dofmap.mesh.cells;
    coords = dofmap.mesh.vertex_coords;
    cell_dofs = dofmap.cell_dofs;
    num_edges = dofmap.element.simplex.num_entities(1);
    tangents = dofmap.element.simplex.edge_tangents;
    cell_edges = dofmap.mesh.get_connectivity(dim, 1);

    % Pick quadrature rule on edge
    quad_degree = dofmap.element.order + degree_rise;  % FIXME: -(dim-1)?
    [x, w] = fe.make_quadrature_rule(1, quad_degree);
    num_quad_points = size(x, 1);

    % Transform quadrature rule to reference facets
    xe = zeros(num_quad_points, dim, num_edges);
    we = zeros(num_quad_points, num_edges);
    basis = zeros(local_element_dim, num_quad_points, num_edges);
    for e = 1:num_edges
        transform = dofmap.element.simplex.get_edge_transform(e);
        for k = 1:num_quad_points
            xe(k, :, e) = transform(x(k, :).');
            we(k, e) = w(k)*dofmap.element.simplex.edge_length(e);
            basis(:, k, e) = dofmap.element.tabulate_basis(xe(k, :, e)) ...
                             *tangents(:, e);
        end
    end

    % Preallocate temporaries
    offset = zeros(dim, 1);
    jac = zeros(dim, dim);
    temp = zeros(local_element_dim, 1);
    tau = zeros(dim, 1, 'double');

    % Allocate return value
    b = zeros(dofmap.dim, 1);

    % Loop over cells
    for c = 1:num_cells

        % Find marked local edges
        local_edges = full(markers(cell_edges(:, c)) == marker_value);

        % If no marked continue to another cell
        if all(~local_edges)
            continue
        end

        % Compute geometric quantities
        offset(:) = coords(:, cells(dim+1, c));
        jac(:, :) = coords(:, cells(1:dim, c)) - coords(:, cells(dim+1, c));

        % Pull back coefficient to reference element
        assert(iscolumn(offset) && numel(offset) == dim);
        g_hat = @(xhat, tau) g(jac*xhat.' + offset, tau);

        % Zero cell vector
        temp(:) = 0;

        % Loop over local marked edges and quadrature points
        for e = find(local_edges).'

            % Compute physical edge tangents
            tau(:) = jac*tangents(:, e);
            tau(:) = tau/norm(tau);

            % Quadrature
            for k = 1:num_quad_points
                % TODO: Should loop only through edge-supported basis functions?
                temp(:) = temp(:) + ...
                    we(k, e)*g_hat(xe(k, :, e), tau)*basis(:, k, e);
            end
        end

        % Add cell contribution to global vector
        % NB: Assuming continuity on interior facet and
        %     overwriting contribution from neighboring cells
        b(cell_dofs(:, c)) = temp;
    end
end
