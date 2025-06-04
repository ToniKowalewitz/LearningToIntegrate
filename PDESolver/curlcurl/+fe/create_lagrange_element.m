function element = create_lagrange_element(dim, order)
    % Create instance of Lagrange element.
    %
    % SYNTAX
    %   element = create_lagrange_element(dim, order)
    %
    % INPUT/OUTPUT PARAMETERS
    %
    %   dim     ... Dimension of the reference cell.
    %   order   ... Order of the element.
    %   element ... Instance of the H1Element class.

    reference_simplex = fe.ReferenceSimplex(dim);

    % Compute compatible entity dofmaps and lattice points
    entity_dofs = cell(dim+1, 1);  % dofmap (row index is dimension of entities + 1)
    pts = zeros(0, dim);           % lattice points (row index is dof index)
    dof_cntr = 0;
    for d = 0:dim
        lat = reference_simplex.create_lattice_interior(d, order);
        num_dofs_per_entity = size(lat, 1);
        num_entities_per_dim = size(lat, 3);
        num_dofs_per_dim = num_dofs_per_entity*num_entities_per_dim;
        dofs = dof_cntr+1:dof_cntr+num_dofs_per_dim;
        entity_dofs{d+1} = reshape(dofs, num_dofs_per_entity, num_entities_per_dim);
        pts(dofs, :) = reshape(permute(lat, [1, 3, 2]), [num_dofs_per_dim, dim]);
        dof_cntr = dof_cntr + num_dofs_per_dim;
    end

    % Dual basis as evaluation at lattice points
    function dofs = dual_basis(v)
        for j = 1:size(pts, 1)
            dofs(:, j) = v(pts(j, :));  %#ok<AGROW>
        end
    end

    % Build monomial basis
    monomials_multiindices_barycentric = fe.barycentric_lattice(dim, order);
    monomials_multiindices = monomials_multiindices_barycentric(:, 1:dim);
    if verLessThan('matlab', '9.3')
        tabulate_basis_handle = @tabulate_basis_compat;
    else
        tabulate_basis_handle = @tabulate_basis;
    end

    function Phi = tabulate_basis(x)
        assert(isrow(x) && length(x) == size(monomials_multiindices, 2));
        Phi = prod(x.^monomials_multiindices, 2);
    end

    % NB: implementation compatible with R2017a
    function Phi = tabulate_basis_compat(x)
        assert(isrow(x));
        for j = 1:size(monomials_multiindices, 1)
            Phi(j, :) = prod(x.^monomials_multiindices(j, :));  %#ok<AGROW>
        end
    end

    % Functions are scalar
    value_shape = 1;

    % Build nodal element from its defining data
    element = fe.H1Element('Lagrange', 'affine', order, ...
                           reference_simplex, value_shape, entity_dofs,  ...
                           tabulate_basis_handle, @dual_basis);
end
