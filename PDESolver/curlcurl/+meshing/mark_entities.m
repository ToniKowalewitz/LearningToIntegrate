function markers = mark_entities(mesh, dim, markers, where, value, varargin)
    % Create/modify entity markers according to provided region.
    %
    % SYNTAX
    %
    %   markers = mark_entities(mesh, dim, markers, where, value)
    %
    % INPUT PARAMETERS
    %   mesh    ... Instance of Mesh class.
    %   dim     ... Dimension of entities.
    %   markers ... Empty vector or vector to be modified
    %               by marking.
    %   where   ... Logical-valued function of '(x, b)', where
    %               'x' are coordinates, 'b' is true for boundary
    %               entity, otherwise false. The function is invoked
    %               on vertices of entities. Entities where the function
    %               evaluates 'true' are marked.
    %   value   ... Value to set the marked entities.
    %   varargin... Optional name-value pairs:
    %                   'check_barycenter', bool (default false):
    %                        If true, check also whether barycenter
    %                        of the entity fulfills 'where' condition;
    %                        otherwise just the entity vertices.
    %
    %  EXAMPLE
    %
    %    mesh = meshing.generate_half_disc_mesh(16);
    %    mesh.compute_connectivity(mesh.dim, mesh.dim-1);
    %    mesh.compute_boundary_facets();
    %    tol = 1e-10;
    %    boundary1 = @(x, b) b && norm(x) > 1-tol;
    %    markers = meshing.mark_entities(mesh, mesh.dim-1, [], boundary1, 1);
    %    boundary2 = @(x, b) b && x(2) < tol;
    %    markers = meshing.mark_entities(mesh, mesh.dim-1, markers, boundary2, 2);
    %    interior = @(x, b) ~b;
    %    markers = meshing.mark_entities(mesh, mesh.dim-1, markers, interior, -1);
    %    meshing.plot_mesh(mesh, 'facet_markers', markers);

    % Parse name-value arguments
    p = inputParser();
    p.addParameter('check_barycenter', false);
    p.parse(varargin{:});
    check_barycenter = p.Results.check_barycenter;

    % Fetch data of mesh
    coords = mesh.vertex_coords;
    num_cells = mesh.num_entities(mesh.dim);
    num_entities = mesh.num_entities(dim);
    c2v = mesh.cells;
    c2e = mesh.get_connectivity(mesh.dim, dim);
    boundary_facets = mesh.get_boundary_facets();

    % Fetch topological data of reference cell
    reference_simplex = fe.ReferenceSimplex(mesh.dim);
    if dim == mesh.dim
        num_facets_per_cell = reference_simplex.num_entities(mesh.dim-1);
        % NB: All facets are connected to the cell (entity index 1)
        % TODO: Implement this in ReferenceSimplex?
        f2e_ref = ones(1, num_facets_per_cell);
    else
        f2e_ref = reference_simplex.get_connectivity(mesh.dim-1, dim);
    end
    e2v_ref = reference_simplex.get_connectivity(dim, 0);
    num_vertices_per_entity = size(e2v_ref, 1);
    num_entities_per_cell = reference_simplex.num_entities(dim);

    % Preallocate dense result (if not passed in)
    if isempty(markers)
        markers = zeros(num_entities, 1);
    else
        assert(numel(markers) == num_entities, ...
               'Expected ''markers'' to be vector of length %d or empty vector', ...
               num_entities);
    end

    % Allocate temporaries
    entities_on_boundary = zeros(num_entities_per_cell, 1, 'logical');
    vertices_inside = zeros(num_vertices_per_entity, 1, 'logical');
    entity_vertices_coords = zeros(mesh.dim, num_vertices_per_entity);

    % Loop over mesh cells
    for c = 1:num_cells

        % Figure out which local entities are on boundary
        entities_on_boundary(:) = 0;
        entities_on_boundary(f2e_ref(:, boundary_facets(:, c))) = 1;

        % Loop over local entities
        for e = 1:num_entities_per_cell

            % Fetch coordinates of entity vertices
            entity_vertices_coords(:, :) = coords(:, c2v(e2v_ref(:, e), c));

            % Figure out which entity vertices are inside
            for v = 1:num_vertices_per_entity
                vertices_inside(v) = where(entity_vertices_coords(:, v), ...
                                           entities_on_boundary(e));
            end

            % Mark entity if all its vertices are inside
            if all(vertices_inside)
                if check_barycenter
                    % Additionally check barycenter if required
                    if where(mean(entity_vertices_coords(:, :), 2), ...
                             entities_on_boundary(e))
                        markers(c2e(e, c)) = value;
                    end
                else
                    markers(c2e(e, c)) = value;
                end
            end
        end
    end

    % Assemble out zeros
    [i, j, v] = find(markers);
    markers = sparse(i, j, v, num_entities, 1);
end
