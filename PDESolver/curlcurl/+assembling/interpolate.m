function x = interpolate(dofmap, value)
    % Interpolate function into FE space.
    %
    % SYNTAX
    %   x = interpolate(dofmap, value)
    %
    % INPUT PARAMETERS
    %   dofmap ... Struct representing FE space
    %   value  ... Function handle representing interpolated function
    %              with signature:
    %
    %                vals = value(x, c)
    %
    %              where
    %
    %                x    ... coordinates (in physical domain); row vector
    %                         of shape [dim, 1]
    %                c    ... cell index; function does not need to use it
    %                vals ... returned value
    %                            FIXME: of what shape?
    %
    % OUTPUT PARAMETER
    %   x ... Vector with expansion coefficients of FE function

    assert(isstruct(dofmap), 'Expected structure as argument 1');
    assert(isa(value, 'function_handle'), 'Expected function handle as arguent 2');

    % Fetch data
    element = dofmap.element;
    interpolation_operator = @element.evaluate_dual_basis;
    dim = dofmap.mesh.dim;
    cells = dofmap.mesh.cells;
    coords = dofmap.mesh.vertex_coords;
    cell_dofs = dofmap.cell_dofs;

    % Allocate result
    x = zeros(dofmap.dim, 1);

    % Preallocate temporaries
    B = zeros(dim, dim);
    b = zeros(dim, 1);
    dof_values = zeros(dofmap.element.fe_space_dim, 1, 'double');
    dof_indices = zeros(dofmap.element.fe_space_dim, 1, 'uint32');

    % Pullback to reference element
    if strcmp(element.mapping, 'affine')
        f_hat = @(B, b, c) @(xhat) value(B*xhat.' + b, c);
    elseif strcmp(element.mapping, 'covariant')
        f_hat = @(B, b, c) @(xhat) value(B*xhat.' + b, c)*B;
    elseif strcmp(element.mapping, 'contravariant')
        f_hat = @(B, b, c) @(xhat) det(B)*value(B*xhat.' + b, c)/B.';
    end

    % Loop over cells
    for c = 1:size(dofmap.cell_dofs, 2)

        % Populate geometry for pullback
        B(:, :) = coords(:, cells(1:dim, c)) - coords(:, cells(dim+1, c));
        b(:, :) = coords(:, cells(dim+1, c));

        % Compute dof values and indices; apply local 'interpolation_operator'
        % to pull-back 'f_hat' in order to get the expansion coefficients
        % 'dof_values'
        dof_values(:) = interpolation_operator(f_hat(B, b, c));
        dof_indices(:) = cell_dofs(:, c);

        % Store cell interpolant into global vector
        x(dof_indices) = dof_values;

    end

end
