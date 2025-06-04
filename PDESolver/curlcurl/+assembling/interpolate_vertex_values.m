function vertex_values = interpolate_vertex_values(dofmap, x)
    % Interpolate FE function at mesh vertices.
    %
    % Second argument is a vector of expansion coefficients
    % of FE function on a space defined by the provided dofmap.

    % FIXME: This is a path working if the element has vertex dofs.
    % If not we have to take whole expansion and evaluate basis
    % members at vertices.

    assert(strcmp(dofmap.element.mapping, 'affine'));

    % Fetch data from dofmap
    num_vertices = dofmap.mesh.num_entities(0);
    cell_dofs = dofmap.cell_dofs;
    cell_vertex_conn = dofmap.mesh.get_connectivity(dofmap.mesh.dim, 0);
    local_vertex_dofs = dofmap.element.get_entity_dofs(0);

    % Sanity checks
    if isempty(local_vertex_dofs)
        error(['This element does not have vertex dofs; '...
               'Interpolation not implemented yet!'])
    end
    num_vertices_per_cell = dofmap.element.simplex.num_entities(0);
    assert(size(local_vertex_dofs, 1) == 1, 'Did not expect multiple dofs per vertex');
    assert(size(local_vertex_dofs, 2) == num_vertices_per_cell);
    assert(size(cell_vertex_conn, 1) == num_vertices_per_cell);

    % Build vertex values
    vertex_values = zeros(1, num_vertices);
    vertex_values(cell_vertex_conn(:, :)) = x(cell_dofs(local_vertex_dofs(:), :));
end
