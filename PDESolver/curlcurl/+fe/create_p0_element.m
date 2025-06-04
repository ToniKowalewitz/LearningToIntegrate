function element = create_p0_element(dim)
    % Create instance of P0 element.
    %
    % SYNTAX
    %   element = create_p0_element(dim)
    %
    % INPUT/OUTPUT PARAMETERS
    %
    %   dim     ... Dimension of the reference cell.
    %   element ... Instance of the FiniteElement class.

    reference_simplex = fe.ReferenceSimplex(dim);

    % Compute entity dofmaps
    entity_dofs = cell(dim+1, 1);
    for d = 0:dim-1
        num_entities_per_dim = reference_simplex.num_entities(d);
        entity_dofs{d+1} = reshape([], 0, num_entities_per_dim);
    end
    entity_dofs{dim+1} = 1;

    % Dual basis is evaluation at barycenter
    num_vertices = size(reference_simplex.vertex_coords, 2);
    barycenter = sum(reference_simplex.vertex_coords, 2).' / num_vertices;
    assert(isrow(barycenter));
    dual_basis_handle = @(v) v(barycenter);

    % Basis is a constant
    tabulate_basis_handle = @(p) 1;

    % Functions are scalar
    value_shape = 1;

    % Build nodal element from its defining data
    element = fe.FiniteElement('DiscontinuousLagrange', 'affine', 0, ...
                               reference_simplex, value_shape, entity_dofs,  ...
                               tabulate_basis_handle, dual_basis_handle);
end
