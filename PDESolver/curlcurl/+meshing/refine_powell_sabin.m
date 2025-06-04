function mesh = refine_powell_sabin(mesh)
    % Compute Powell-Sabin refinement of given mesh

    if mesh.dim ~= 2
        error('Powell-Sabin refinement only supported for dimension 2');
    end

    % Fetch data
    num_cells = mesh.num_entities(mesh.dim);
    num_facets = mesh.num_entities(mesh.dim-1);
    num_vertices = mesh.num_entities(0);
    c2v = mesh.get_connectivity(mesh.dim, 0);
    c2f = mesh.get_connectivity(mesh.dim, mesh.dim-1);
    f2v_ref = fe.ReferenceSimplex(mesh.dim).get_connectivity(mesh.dim-1, 0);

    % Add new vertices at facets and cell incenters
    cell_incenters = mesh.get_cell_incenters();
    facet_points = compute_facet_points(mesh, cell_incenters);
    vertex_coords = zeros(2, num_vertices + num_facets + num_cells);
    vertex_coords(:, 1:num_vertices) = mesh.vertex_coords;
    vertex_coords(:, num_vertices+1:num_vertices+num_facets) = facet_points;
    vertex_coords(:, num_vertices+num_facets+1:end) = cell_incenters;

    % Compute new topology
    cells = zeros(3, 6*num_cells, 'uint32');
    for c = 1:num_cells

        % Indices to original vertices
        vv = c2v(:, c);

        % Indices to new facet points
        fv = c2f(:, c) + num_vertices;

        % Index to cell incenter
        cv = c + num_vertices + num_facets;

        cells(:, 6*c - 5) = [vv(f2v_ref(1, 1)), fv(1), cv];
        cells(:, 6*c - 4) = [vv(f2v_ref(2, 1)), fv(1), cv];
        cells(:, 6*c - 3) = [vv(f2v_ref(1, 2)), fv(2), cv];
        cells(:, 6*c - 2) = [vv(f2v_ref(2, 2)), fv(2), cv];
        cells(:, 6*c - 1) = [vv(f2v_ref(1, 3)), fv(3), cv];
        cells(:, 6*c - 0) = [vv(f2v_ref(2, 3)), fv(3), cv];

    end

    % Create new mesh object
    mesh = meshing.Mesh(2, vertex_coords, cells);

end


function coords = compute_facet_points(mesh, cell_incenters)
    % Compute new facet points for Powell-Sabin refinement
    %
    % For exterior facets these are facet midpoints, for interior
    % facets the intersection of adjacent cell incenters with
    % the facet are returned.

    assert(mesh.dim == 2);

    num_facets = mesh.num_entities(mesh.dim-1);
    f2c = mesh.get_connectivity(mesh.dim-1, mesh.dim);
    f2v = mesh.get_connectivity(mesh.dim-1, 0);
    vertex_coords = mesh.vertex_coords;

    % Allocate result
    coords = zeros(mesh.dim, num_facets);

    % Loop over facets
    for f = 1:num_facets
        exterior = any(f2c(:, f) == 0);
        if exterior
            % Facet midpoint
            coords(:, f) = sum(vertex_coords(:, f2v(:, f)), 2) / mesh.dim;
        else
            % Facet vertices
            fv = vertex_coords(:, f2v(:, f));

            % Incenters of adjacent cells
            cv = cell_incenters(:, f2c(:, f));

            % Intersection
            coords(:, f) = intersect_two_segments(fv(:, 1), fv(:, 2), cv(:, 1), cv(:, 2));
        end
    end

end


function coord = intersect_two_segments(x1, x2, x3, x4)
    % Find intersection of two segments (x1, x2), (x3, x4)
    %
    % See: https://en.wikipedia.org/w/index.php?title=Line%E2%80%93line_intersection&oldid=926557056

    t = ...
        ((x1(1)-x3(1))*(x3(2)-x4(2))-(x1(2)-x3(2))*(x3(1)-x4(1))) ...
        / ...
        ((x1(1)-x2(1))*(x3(2)-x4(2))-(x1(2)-x2(2))*(x3(1)-x4(1)));

    coord = x1 + t*(x2-x1);
end
