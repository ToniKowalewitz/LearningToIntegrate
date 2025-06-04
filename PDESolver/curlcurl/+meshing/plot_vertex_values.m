function plot_vertex_values(mesh, vertex_values)
    % Plot 2D mesh with vertex elevation/coloring as trisurf
    %
    % SYNTAX
    %   plot_vertex_values(mesh, vertex_values)
    %
    % INPUT PARAMETER
    %   mesh          ... Mesh
    %   vertex_values ... vector with values at vertices

    if mesh.dim ~= 2
        error('not implemented');
    end

    assert(numel(vertex_values) == mesh.num_entities(0));

    axis equal;
    axis tight;

    trisurf(mesh.cells.', ...
            mesh.vertex_coords(1, :), mesh.vertex_coords(2, :), ...
            vertex_values(:));
end
