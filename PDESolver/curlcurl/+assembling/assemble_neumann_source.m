function b = assemble_neumann_source(dofmap, markers, marker_value, g, degree_rise)
    % Assemble Neumann exterior surface integral
    %
    % In affine case:
    %
    %   \int_{markers==marker_value} g(x, n) \phi_j ds,
    %
    % in contravariant case:
    %
    %   \int_{markers==marker_value} g(x, n) \phi_j\cdot n ds,
    %
    % where \phi_j are H^1-conforming (affine) or
    % H(div)-conforming (contravariant) basis functions.
    % Integration domain are the exterior boundary facets where
    % the supplied markers take the provided value.
    %
    % INPUT PARAMETER
    %   dofmap       ... Dofmap of H^1 or H(div) space
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
    %   b            ... Vector

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
    quad_degree = dofmap.element.order + degree_rise;  % FIXME: -1?
    [x, w] = fe.make_quadrature_rule(dim-1, quad_degree);
    num_quad_points = size(x, 1);
    area = 1/factorial(dim-1);

    % Scaling factors in integrals
    switch dofmap.element.mapping
    case 'affine'
        basis_multiplier = @(f) 1;
        integral_multiplier = @(S, detJ) S*abs(detJ);
    case 'contravariant'
        basis_multiplier = @(f) normals(:, f);
        integral_multiplier = @(S, detJ) sign(detJ);
    case 'covariant'
        % TODO: What is right here? Cross product with normal?
        error('not yet implemented!');
    otherwise
        error('unexpected mapping ''%s''', dofmap.element.mapping);
    end

    % Transform quadrature rule to reference facets
    xf = zeros(num_quad_points, dim, num_facets);
    wf = zeros(num_quad_points, num_facets);
    basis = zeros(local_element_dim, num_quad_points, num_facets);
    for f = 1:num_facets
        transform = dofmap.element.simplex.get_facet_transform(f);
        for k = 1:num_quad_points
            xf(k, :, f) = transform(x(k, :).');
            wf(k, f) = w(k)*dofmap.element.simplex.facet_area(f)/area;
            basis(:, k, f) = ...
                dofmap.element.tabulate_basis(xf(k, :, f))*basis_multiplier(f);
        end
    end

    % Preallocate temporaries
    offset = zeros(dim, 1);
    jac = zeros(dim, dim);
    temp = zeros(local_element_dim, 1);
    n = zeros(dim, 1, 'double');

    % Allocate return value
    b = zeros(dofmap.dim, 1);

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
        detJ = det(jac);

        % Pull back coefficient to reference element
        assert(iscolumn(offset) && numel(offset) == dim);
        g_hat = @(xhat, n) g(jac*xhat.' + offset, n, c);

        % Zero cell vector
        temp(:) = 0;

        % Loop over local marked facets and quadrature points
        for f = find(local_facets).'

            % Compute physical outer normal and the scaling factor S
            n(:) = jac.'\normals(:, f);
            S = norm(n);
            n(:) = n/S;

            % Quadrature
            for k = 1:num_quad_points
                % TODO: Should loop only through facet-supported basis functions?
                temp(:) = temp(:) + ...
                    wf(k, f)*g_hat(xf(k, :, f), n)*basis(:, k, f) ...
                    *integral_multiplier(S, detJ);
            end
        end

        % Add cell contribution to global vector
        b(cell_dofs(:, c)) = b(cell_dofs(:, c)) + temp;
    end
end
