function plot_mesh(mesh, varargin)
    % Plot mesh, optionally with cell, facet, edge, and/or vertex markers
    %
    % SYNTAX
    %   plot_mesh(mesh[, 'Name', Value, ...])
    %
    % INPUT PARAMETER
    %   mesh           ... Object of Mesh class
    %
    % OPTIONAL PARAMETERS (NAME-VALUE PAIRS)
    %   cell_markers   ... Vector [num_cells x 1]
    %   facet_markers  ... (sparse) vector [num_facets x 1]
    %   edge_markers   ... (sparse) vector [num_edges x 1]
    %   vertex_markers ... (sparse) vector [num_vertices x 1]
    %   plane          ... selects only the cells intersecting a plane given by
    %                      the implicit equation:
    %
    %                          plane(x(1:3, :)) == 0;
    %
    %                      'plane' is a function handle of signature:
    %
    %                          double = function(double(1:3, :));
    %
    %                      for example, xz-plane intersecting point [0; 1/2; 0]:
    %
    %                          plane = @(x) 0*x(1,:) + 1*x(2,:) + 0*x(3,:) - 1/2;
    %
    %                      the function shall be vectorized along the second
    %                      axis of x.
    %   preferred_side ... real scalar; when positive (negative), cells on the
    %                      positive (negative) side of the plane are preferred,
    %                      i.e., cells which lie on the non-positive
    %                      (non-negative) side are ignored;
    %                      if zero (default), both sides are treated equally
    %
    % NOTES
    %   Facet, edge, and vertex markers will be typically sparse
    %   vectors with integer values but this is not enforced.
    %   In particular, zeros are ignored rather than plotted
    %   and legend features labels with all (integer) values.
    %
    %   On the other hand, cell markers are interpreted as
    %   continuous, double-valued field and it is accompanied
    %   by a color bar rather than a legend with values.
    %
    % EXAMPLE (2D)
    %   mesh = meshing.generate_unit_cube_mesh([3, 3]);
    %   mesh.compute_connectivity(1, 0);
    %   cm = mod(1:mesh.num_entities(mesh.dim), 3);
    %   em = mod(1:mesh.num_entities(1), 4);
    %   vm = mod(1:mesh.num_entities(0), 5);
    %   meshing.plot_mesh(mesh, 'cell_markers', cm, 'edge_markers', em, ...
    %                     'vertex_markers', vm);
    %
    % EXAMPLE (3D)
    %   mesh = meshing.generate_unit_cube_mesh([3, 3, 3]);
    %   mesh.compute_connectivity(1, 0);
    %   cm = mod(1:mesh.num_entities(mesh.dim), 3);
    %   em = mod(1:mesh.num_entities(1), 4);
    %   vm = mod(1:mesh.num_entities(0), 5);
    %   pl = @(x) 3*x(1,:) + 1*x(2,:) + 2*x(3,:) - 2;
    %   meshing.plot_mesh(mesh, 'cell_markers', cm, 'edge_markers', em, ...
    %                     'vertex_markers', vm, 'plane', pl);
    %
    % EXAMPLE (3D)
    %   mesh = meshing.generate_unit_cube_mesh([3, 3, 3]);
    %   pl = @(x) x(2,:) - 1/3;
    %   meshing.plot_mesh(mesh, 'plane', pl, 'preferred_side', +1);

    assert(isa(mesh, 'meshing.Mesh'));

    isvec = @(x) isempty(x) || isvector(x) && isnumeric(x) && isreal(x);
    isfun = @(x) isempty(x) || isa(x, 'function_handle');
    isnum = @(x) isscalar(x) && isnumeric(x) && isreal(x);

    parser_obj = inputParser();
    parser_obj.addParameter('cell_markers',   [], isvec);
    parser_obj.addParameter('facet_markers',  [], isvec);
    parser_obj.addParameter('edge_markers',   [], isvec);
    parser_obj.addParameter('vertex_markers', [], isvec);
    parser_obj.addParameter('plane',          {}, isfun);
    parser_obj.addParameter('preferred_side', 0,  isnum);
    parse(parser_obj, varargin{:});
    args = parser_obj.Results;

    if ~isempty(args.facet_markers) && mesh.dim == 3
        error('plotting of tetrahedra facet markers not yet implemented!');
    end

    % Facets are edges in 2D
    if ~isempty(args.facet_markers) && mesh.dim == 2
        if ~isempty(args.edge_markers)
            error('supply only one of "facet_markers" and "edge_markers" arguments!');
        else
            args.edge_markers = args.facet_markers;
        end
    end

    if ~isempty(args.plane) && mesh.dim ~= 3
        error('intersection with plane only implemented in 3D!');
    end

    hold on;
    axis equal;
    axis tight;

    annotated_plots = [];

    % Plot mesh and cell markers (if given)
    switch mesh.dim
    case 2
        assert(isempty(args.plane));
        plot_cell_markers_2d(mesh, args.cell_markers);
    case 3
        plot_cell_markers_3d(mesh, args.cell_markers, args.plane, args.preferred_side);
    end

    % Plot edge markers if given
    if ~isempty(args.edge_markers)
        plots = plot_edge_markers(mesh, args.edge_markers);
        annotated_plots = [annotated_plots(:); plots(:)];
    end

    % Plot vertex markers if given
    if ~isempty(args.vertex_markers)
        plots = plot_vertex_markers(mesh, args.vertex_markers);
        annotated_plots = [annotated_plots(:); plots(:)];
    end

    % Issue legend if appropriate
    if ~isempty(annotated_plots)
        legend(annotated_plots, 'Location', 'BestOutside');
    end

    hold off;
end


function plot_cell_markers_2d(mesh, cell_markers)
    if isempty(cell_markers)
        coords = mesh.vertex_coords;
        triplot(double(mesh.cells.'), coords(1, :), ...
                coords(2, :), 'Color', 'k');
    else
        patch('Faces', mesh.cells.', ...
              'Vertices', mesh.vertex_coords.', ...
              'FaceVertexCData', cell_markers(:), ...
              'FaceColor', 'flat', ...
              'EdgeColor', 'none');
        colorbar();
    end
end


function plot_cell_markers_3d(mesh, cell_markers, plane, preferred_side)
    num_cells = mesh.num_entities(mesh.dim);

    if isempty(plane)
        % Not using plane, using all cells
        intersected_cells = true(1, num_cells);
    else
        if preferred_side > 0
            % Cells with x<=0 ignored
            intersects_ = @(sdist) any(sdist>0) && any(sdist<=0);
        elseif preferred_side < 0
            % Cells with x>=0 ignored
            intersects_ = @(sdist) any(sdist>=0) && any(sdist<0);
        else
            % Cells on both sides treated equally
            intersects_ = @(sdist) any(sdist>=0) && any(sdist<=0);
        end

        % Find cells intersecting plane
        intersected_cells = false(1, num_cells);
        for c = 1:num_cells
            cell_coords = mesh.vertex_coords(:, mesh.cells(:, c));

            % Signed distance of each vertex in cell to the plane
            sd = plane(cell_coords);

            % Cell intersects plane if there are vertices on different side
            intersected_cells(c) = intersects_(sd);
        end
    end

    % Plot
    if isempty(cell_markers)
        cell_markers = ones(1, mesh.num_entities(mesh.dim));
        tetramesh(mesh.cells(:, intersected_cells).', ...
                  mesh.vertex_coords.', ...
                  cell_markers(intersected_cells));
    else
        tetramesh(mesh.cells(:, intersected_cells).', ...
                  mesh.vertex_coords.', ...
                  cell_markers(intersected_cells), ...
                  'EdgeColor', 'none');
        colorbar();
    end

    % Set view perpendicular to plane
    if ~isempty(plane)
        normal = plane(eye(3)) - plane(zeros(3));
        view(-normal);
    end
end


function plts = plot_edge_markers(mesh, edge_markers)
    coords = mesh.vertex_coords;
    edges_to_vertices = mesh.get_connectivity(1, 0);

    % Build unique marker values and mapping back to marked edges
    [marked_edges, ~, values] = find(edge_markers(:));
    [values, ~, inds] = unique(values);

    % Choose Matlab implementation
    switch mesh.dim
    case 2
        plot_ = @(edges, value, color) plot(...
            reshape(coords(1, edges_to_vertices(:, edges)), 2, []), ...
            reshape(coords(2, edges_to_vertices(:, edges)), 2, []), ...
            'DisplayName', sprintf('%d', value), ...
            'Color', color, 'LineWidth', 3);
    case 3
        plot_ = @(edges, value, color) plot3(...
            reshape(coords(1, edges_to_vertices(:, edges)), 2, []), ...
            reshape(coords(2, edges_to_vertices(:, edges)), 2, []), ...
            reshape(coords(3, edges_to_vertices(:, edges)), 2, []), ...
            'DisplayName', sprintf('%d', value), ...
            'Color', color, 'LineWidth', 3);
    end

    cmap = jet(numel(values));
    plts = gobjects(numel(values), 1);

    % Plot value-by-value
    for j = 1:numel(values)
        edges = marked_edges(inds == j);
        value = values(j);
        color = cmap(j, :);
        p = plot_(edges, value, color);

        % Remember only first object for each value to have concise legend
        plts(j) = p(1);
    end
end


function plts = plot_vertex_markers(mesh, vertex_markers)
    coords = mesh.vertex_coords;

    % Build unique marker values and mapping back to marked vertices
    [marked_vtx, ~, values] = find(vertex_markers(:));
    [values, ~, inds] = unique(values);

    % Choose Matlab implementation
    switch mesh.dim
    case 2
        plot_ = @(verts, value, color) plot(...
            coords(1, verts), coords(2, verts), ...
            '.', 'DisplayName', sprintf('%d', value), ...
            'Color', color, 'MarkerSize', 23);
    case 3
        plot_ = @(verts, value, color) plot3(...
            coords(1, verts), coords(2, verts), coords(3, verts), ...
            '.', 'DisplayName', sprintf('%d', value), ...
            'Color', color, 'MarkerSize', 23);
    end

    cmap = jet(numel(values));
    plts = gobjects(numel(values), 1);

    % Plot value-by-value
    for j = 1:numel(values)
        verts = marked_vtx(inds == j);
        value = values(j);
        color = cmap(j, :);
        p = plot_(verts, value, color);

        % Remember only first object for each value to have concise legend
        plts(j) = p(1);
    end
end
