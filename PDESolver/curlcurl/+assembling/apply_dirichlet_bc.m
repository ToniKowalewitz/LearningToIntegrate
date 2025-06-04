function [A, b] = apply_dirichlet_bc(A, b, dofs, values, varargin)
    % Apply Dirichlet BC to linear system Ax = b
    %
    % Function accepts name-value pairs:
    %
    %     'symmetric', bool (default true):
    %         If false, set BC rows to unit diagonal
    %         and set BC values into BC dofs of vector.
    %         If true, move BC columns times BC values
    %         to vector which results in symmetric system.
    %         It is correct only when both A, b are
    %         supplied (unless BC values are zero).
    %
    % Supplied A or b can be empty in which case only work on b or A,
    % respectively, is performed; this makes sense only if
    % 'symmetric'=false or BC is zero.

    % Parse name-value arguments
    p = inputParser();
    p.addParameter('symmetric', true);
    p.parse(varargin{:});
    symmetric = p.Results.symmetric;

    % Special case of zero BC
    iszero = isscalar(values) && values == 0;

    % Special case of empty BC
    if isempty(dofs)
        assert(isempty(values) || iszero);
        return
    end

    % Sanity check
    assert(iscolumn(dofs), 'provide column vector "dofs"');
    assert(iszero || size(values, 1) == numel(dofs), ...
           'number of rows in "values" and "dofs" has to match');
    assert(iszero || isempty(b) || size(values, 2) == size(b, 2), ...
           'number of columns in "values" and "b" has to match');
    assert(numel(dofs) == numel(unique(dofs)), 'provide unique "dofs"');

    % Move BC columns contribution to right-hand side to preserve symmetry
    if symmetric && ~iszero
        if ~isempty(A) && ~isempty(b)
            b(:, :) = b - A(:, dofs)*values;
        else
            error('Can only apply non-zero BC simultaneously to "A" and "b"');
        end
    end

    % Set BC values in vector
    if ~isempty(b)
        b(dofs, :) = values;
    end

    % Put ones on BC diagonal and zero BC rows and optionally BC columns
    if ~isempty(A)
        A = matrix_kernel(A, dofs, symmetric);
    end

end


function A = matrix_kernel(A, dofs, zero_columns)
    % Kernel for manipulating BC rows, columns, and diagonal
    %
    % Equivalent to the following naive implementation:
    %
    %   A(dofs, :) = 0;
    %   if zero_columns
    %       A(:, dofs) = 0;
    %   end
    %   A(sub2ind(size(A), dofs, dofs)) = 1;
    %
    % The naive implementation has quadratic complexity,
    % hence this specialized implementation is needed. See
    % https://gitlab.hrz.tu-chemnitz.de/geosax/curlcurl/issues/8

    % For dense A (if ever needed) the naive implementation should be used
    assert(issparse(A));

    % Disassemble matrix into coordinate representation
    [I, J, V] = find(A);
    [m, n] = size(A);

    % Free unneeded memory
    clear A;

    % Zero BC columns if required
    if zero_columns
        K = ismember(J, dofs);
        V(K) = 0;
    end

    % Zero BC rows
    K = ismember(I, dofs);
    V(K) = 0;

    % Set BC diagonal elements to one
    V(K & I==J) = 1;

    % Reassemble sparse matrix
    A = sparse(I, J, V, m, n);
end
