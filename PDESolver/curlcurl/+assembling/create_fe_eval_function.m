function fun = create_fe_eval_function(dofmap, vec, derivative)
    % Return function handle for evaluation of FE function
    %
    % SYNTAX
    %   u = create_fe_function(dofmap, vec)
    %
    % INPUT PARAMETERS
    %   dofmap     ... Struct describing the FE space
    %   vec        ... Column vector with expansion coefficients of
    %                  the FE function w.r.t. basis
    %   derivative ... Optional character vector; possible values:
    %                  '', 'grad', 'curl', or 'div'.
    %
    % OUTPUT PARAMETER
    %  fun ... Function handle with signature:
    %
    %            vals = fun(x, c)
    %            vals = fun(x)
    %
    %          Takes either point x and cell c or only point x
    %          and returns value of the FE function. Note that
    %          the second call (without cell c) is costly!

    assert(isstruct(dofmap));
    assert(iscolumn(vec));
    assert(numel(vec) == dofmap.dim);

    if nargin < 3
        derivative = '';
    end

    % Push forward the reference basis (or its derivative)
    basis = get_pushforward(dofmap.element, derivative);

    function c = find_colliding_cell(x)
        % Find first cell colliding with point x

        c = dofmap.mesh.point_location(x(:).');
        if isnan(c)
            error('point ''%s'' not in mesh', mat2str(x));
        end
    end

    % Fetch some data
    dim = dofmap.mesh.dim;
    cells = dofmap.mesh.cells;
    coords = dofmap.mesh.vertex_coords;
    cell_dofs = dofmap.cell_dofs;


    function vals = fun_(x, varargin)
        % Evaluate function at point x or point x and cell c

        % Get or find a colliding cell
        switch nargin
        case 2
            c = varargin{1};
            assert(isscalar(c), 'Expected scalar cell index as argument 2');
        case 1
            c = find_colliding_cell(x);
            assert(isscalar(c));
        otherwise
            error('Wrong number of arguments');
        end

        % Evaluate global basis
        jac = coords(:, cells(1:dim, c)) - coords(:, cells(dim+1, c));
        offset = coords(:, cells(dim+1, c));
        basis_vals = basis(jac, offset, x);

        % Get FE coefficients on cell c
        coeffs = vec(cell_dofs(:, c));

        % Form linear combination
        vals = coeffs.'*basis_vals;

    end

    fun = @fun_;

end


function basis = get_pushforward(element, derivative)

    % FIXME: This code could be part of FiniteElement?

    switch derivative
    case ''
        switch element.mapping
        case 'affine'
            basis = @(B, b, x) element.tabulate_basis((B\(x(:)-b(:))).');
        case 'covariant'
            basis = @(B, b, x) element.tabulate_basis((B\(x(:)-b(:))).')/B;
        case 'contravariant'
            basis = @(B, b, x) element.tabulate_basis((B\(x(:)-b(:))).')*B.'/det(B);
        otherwise
            error('Unexpected element mapping');
        end
    case 'grad'
        switch element.mapping
        case 'affine'
            basis = @(B, b, x) element.tabulate_basis_grad((B\(x(:)-b(:))).')/B;
        otherwise
            error('Unexpected element mapping for grad');
        end
    case 'curl'
        switch element.mapping
        case 'covariant'
            if element.simplex.dim == 2
                basis = @(B, b, x) element.tabulate_basis_curl((B\(x(:)-b(:))).')/det(B);
            else
                basis = @(B, b, x) element.tabulate_basis_curl((B\(x(:)-b(:))).')*B.'/det(B);
            end
        otherwise
            error('Unexpected element mapping for curl');
        end
    case 'div'
        switch element.mapping
        case 'contravariant'
            basis = @(B, b, x) element.tabulate_basis_div((B\(x(:)-b(:))).')/det(B);
        otherwise
            error('Unexpected element mapping for div');
        end
    otherwise
        error('Unexpected derivative specification');
    end

end
