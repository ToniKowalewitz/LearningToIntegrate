function Q = assemble_point_sources(dofmap, points, derivative)
    % Assemble point sources
    %
    % SYNTAX
    %   Q = assemble_point_sources(dofmap, points, derivative)
    %
    % INPUT PARAMETERS
    %   dofmap     ... Struct describing the FE space
    %   points     ... Matrix [dim, num_points]
    %   derivative ... Optional character vector; possible values:
    %                  '', 'grad', 'curl', or 'div'.
    %
    % OUTPUT PARAMETER
    %  Q ... Sparse matrix, evaluation operator

    % Sanity check
    assert(size(points, 1) == dofmap.mesh.dim, ...
        'size(points, 1) has to match dofmap.mesh.dim.');

    if nargin < 3
        derivative = '';
    end

    % Fetch needed data
    num_points = size(points, 2);
    element = dofmap.element;
    cell_dofs = dofmap.cell_dofs;

    % Get pushforward of basis (derivatives)
    [tabulate, vs] = get_pushforward(element, dofmap.mesh, derivative);

    % Compute cells and barycentric coords for all points
    [found_cells, barycentric_coords] = dofmap.mesh.point_location(points.');
    assert(all(not(isnan(found_cells))), 'points not in mesh');

    % Transform coords to reference cell
    % NB: We use here definition of fe.ReferenceSimplex
    reference_coords = barycentric_coords(:, 1:end-1);

    % Allocate storage
    I = zeros(element.fe_space_dim, vs, num_points);
    V = zeros(element.fe_space_dim, vs, num_points);

    % Tabulate basis locally point-by-point, stores values and dof indices
    for j = 1:num_points
        c = found_cells(j);
        xhat = reference_coords(j, :);
        V(:, :, j) = tabulate(xhat, c);
        I(:, :, j) = double(cell_dofs(:, c)) * ones(1, vs);
    end

    % Assemble sparse matrix
    J = ones(element.fe_space_dim, 1)*(1:num_points*vs);
    Q = sparse(I(:), J(:), V(:), dofmap.dim, num_points*vs);
end


function [value, value_shape] = get_pushforward(element, mesh, derivative)

    % FIXME: This code could be part of FiniteElement?

    coords = mesh.vertex_coords;
    cells = mesh.cells;
    dim = mesh.dim;

    function val = transform_value_covariant(val, c)
        B = coords(:, cells(1:dim, c)) - coords(:, cells(dim+1, c));
        val = val/B;
    end

    function val = transform_value_contravariant(val, c)
        B = coords(:, cells(1:dim, c)) - coords(:, cells(dim+1, c));
        val = val*B.'/det(B);
    end

    function val = transform_value_div(val, c)
        B = coords(:, cells(1:dim, c)) - coords(:, cells(dim+1, c));
        val = val/det(B);
    end

    switch derivative
    case ''
        basis = @element.tabulate_basis;
        value_shape = element.value_shape;
        switch element.mapping
        case 'affine'
            value = @(xhat, c) basis(xhat);
        case 'covariant'
            value = @(xhat, c) transform_value_covariant(basis(xhat), c);
        case 'contravariant'
            value = @(xhat, c) transform_value_contravariant(basis(xhat), c);
        otherwise
            error('Unexpected element mapping');
        end
    case 'grad'
        basis = @element.tabulate_basis_grad;
        value_shape = element.value_shape*element.simplex.dim;
        switch element.mapping
        case 'affine'
            value = @(xhat, c) transform_value_covariant(basis(xhat), c);
        otherwise
            error('not implemented');
        end
    case 'curl'
        switch element.mapping
        case 'covariant'
            if element.simplex.dim == 2
                basis = @element.tabulate_basis_curl;
                value = @(xhat, c) transform_value_div(basis(xhat), c);
                % 2d curl turns vectors into scalars
                value_shape = element.value_shape/element.simplex.dim;
            elseif element.simplex.dim == 3
                basis = @element.tabulate_basis_curl;
                value = @(xhat, c) transform_value_contravariant(basis(xhat), c);
                % 3d curl turns vectors into vectors
                value_shape = element.value_shape;
            else
                error('unexpected dimension for curl');
            end
        otherwise
            error('not implemented');
        end
    case 'div'
        basis = @element.tabulate_basis_div;
        value_shape = element.value_shape/element.simplex.dim;
        switch element.mapping
        case 'contravariant'
            value = @(xhat, c) transform_value_div(basis(xhat), c);
        otherwise
            error('not implemented');
        end
    otherwise
        error('Unexpected derivative specification: %s', derivative);
    end

    assert(is_integer(value_shape));
end


function flag = is_integer(num)
    flag = floor(num) == ceil(num);
end
