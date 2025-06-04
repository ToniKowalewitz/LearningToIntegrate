function [dofs, values] = build_dirichlet_dofs(dofmap, varargin)
    % Build array of Dirichlet dof indices and values.
    %
    % SYNTAX
    %   [dofs, values] = build_dirichlet_dofs(dofmap, varargin)
    %
    % INPUT PARAMETER
    %   dofmap ... Struct, containing FE element and mesh object as well as
    %              cell-2-dof mapping.
    %
    % OPTIONAL PARAMETER
    %   nargin == 3:
    %   varargin{1} ... Cell array of functions taking coordinates and
    %                   evaluating as true on boundary facets
    %   varargin{2} ... Cell array of (complex-valued) functions in
    %                   physical coordinates which takes coordinates of
    %                   shape [dim, 1]
    %   nargin == 4:
    %   varargin{1} ... (Sparse) vector marking mesh facets
    %   varargin{2} ... Cell array of marker values where to apply BC
    %   varargin{3} ... Cell array of (complex-valued) functions in
    %                   physical coordinates which takes coordinates of
    %                   shape [dim, 1]
    %
    % OUTPUT PARAMETER
    %   dofs   ... Vector of Dirichlet boundary dof indices.
    %   values ... Vector of corresponding Dirichlet values (i.e.
    %              interpolated expansion coefficients) for dofs.
    %
    % REMARKS
    %   The concept of the cell-loop is:
    %   1. Check which facets of the cell are subject to BC
    %   2. Interpolate Dirichlet value on the cell using the element
    %      interpolation operator (including pull-back to reference element
    %      and transforming back for Piola-transformed elements)
    %   3. Store dof indices subject to BC and corresponding interpolated
    %      values from the previous step in the Map container
    %
    %   The idea of using MATLAB container is:
    %   In general a map `bc_dof -> bbc_value` is needed:
    %   var 1: `bc_dof, bc_value` could be stored as pairs into an array
    %          -> its size is not known a prior
    %          -> duplicate values should be eliminated (resulting from
    %             common dofs of neighbouring cells)
    %   var 2: `bc_dof, bc_value` could be stored in a map container
    %   (used)  -> immediatelly eliminates duplicates
    %           -> does not support complex, so real and imaginary parts
    %              (for complex-valued BC) are collected separately
    %
    % TODO: optimize this function when values not needed (BC is zero)

    % Check input and fetch information
    if nargin == 3
        % Boundary specified by bool-valued function of coords
        where_fun = varargin{1};
        value_fun = varargin{2};
        assert(iscell(where_fun) && iscell(value_fun));
        assert(isrow(where_fun) && isrow(value_fun));
        assert(numel(where_fun) == numel(value_fun));
        compute_local_boundary_facets = compute_from_wherefun(dofmap, ...
                                            where_fun);
    elseif nargin == 4
        % Boundary specified by facet markers and marker values
        facet_markers = varargin{1};
        marker_vals = varargin{2};
        value_fun = varargin{3};
        assert(isvector(facet_markers));
        assert(iscell(marker_vals) && iscell(value_fun));
        assert(isrow(marker_vals) && isrow(value_fun));
        assert(numel(marker_vals) == numel(value_fun));
        compute_local_boundary_facets = compute_from_markers(dofmap, ...
                                            facet_markers, marker_vals);
    else
        error(['Handling of Dirichlet BC supported for "where_fun" or ', ...
               'given facet markers, only.']);
    end
    element = dofmap.element;
    interpolation_operator = @element.evaluate_dual_basis;
    cells = dofmap.mesh.cells;
    coords = dofmap.mesh.vertex_coords;
    dim = dofmap.mesh.dim;
    facet_dofs = dofmap.element.facet_dofs;
    cell_dofs = dofmap.cell_dofs;

    % Preallocate temporaries
    boundary_cell_facets = zeros(1, dim+1, 'logical');  %#ok<PREALL>
    B = zeros(dim, dim);
    b = zeros(dim, 1);
    dof_values = zeros(dofmap.element.fe_space_dim, 1, 'double');
    dof_indices = zeros(dofmap.element.fe_space_dim, 1, 'uint32');
    dirichlet_map_re = containers.Map('KeyType', 'uint32', 'ValueType', 'double');
    dirichlet_map_im = containers.Map('KeyType', 'uint32', 'ValueType', 'double');

    % Pullback to reference element
    if strcmp(dofmap.element.mapping, 'affine')
        f_hat = @(B, b, i) @(xhat) value_fun{i}(B*xhat.' + b);
    elseif strcmp(dofmap.element.mapping, 'covariant')
        f_hat = @(B, b, i) @(xhat) value_fun{i}(B*xhat.' + b)*B;
    elseif strcmp(element.mapping, 'contravariant')
        f_hat = @(B, b, i) @(xhat) det(B)*value_fun{i}(B*xhat.' + b)/B.';
    end

    % Loop over cells with boundary facets
    for c = 1:size(cells, 2)

        % Bool array saying which cell facets are on boundary
        boundary_cell_facets = compute_local_boundary_facets(c);
        [~, marker_idx] = find(boundary_cell_facets);

        % Only if current cell obtains any facets on boundary proceed
        for i = unique(marker_idx).'
            % Populate geometry for pullback
            B(:, :) = coords(:, cells(1:dim, c)) - coords(:, cells(dim+1, c));
            b(:, :) = coords(:, cells(dim+1, c));

            % Compute dof values and indices for all dof of current cell
            dof_values(:) = interpolation_operator(f_hat(B, b, i));
            dof_indices(:) = cell_dofs(:, c);

            % Store into the map only those dofs which match boundary facets
            for f = find(boundary_cell_facets(:, i).')
                for d = facet_dofs(f, :)
                  dirichlet_map_re(dof_indices(d)) = real(dof_values(d));
                  dirichlet_map_im(dof_indices(d)) = imag(dof_values(d));
                end
            end
        end
    end

    % Build arrays
    dofs = cell2mat(dirichlet_map_re.keys);
    values = cell2mat(dirichlet_map_re.values) ...
           + cell2mat(dirichlet_map_im.values)*1i;

    % Reshape and check consistency
    dofs = dofs.';
    values = values.';
    assert(iscolumn(dofs) && iscolumn(values) || isempty(dofs));
end


function func = compute_from_wherefun(dofmap, where)

    % Fetch data
    cells = dofmap.mesh.cells;
    coords = dofmap.mesh.vertex_coords;
    dim = dofmap.mesh.dim;
    facet_vertices = dofmap.element.simplex.get_connectivity(dim-1, 0);
    boundary_facets = dofmap.mesh.get_boundary_facets();

    % Prepare function for computing local boundary facets on cell c
    function local_facets = compute_local_boundary_facets(c)
        % FIXME: Refactor this code!
        % Think only of mesh boundary facets
        % NB: Could be changed...
        local_facets = full(boundary_facets(:, c));
        bnd_facet_idx = find(local_facets).';
        local_facets = repmat(local_facets, 1, numel(where));

        % Keep only facets having all vertices on boundary
        % (according to supplied where function)
        for f = bnd_facet_idx
            facet_coords(:, :) = coords(:, cells(facet_vertices(:, f), c));
            local_facets(f, :) = cellfun(@all,  ...
                                         cellfun(@(w) ...
                                            {arrayfun(@(j) w(facet_coords(:, j)), 1:dim)},  ...
                                         where));
        end
    end

    func = @compute_local_boundary_facets;
end


function func = compute_from_markers(dofmap, facet_markers, facet_tag)

    % Get global cell-facet connectivity
    dim = dofmap.mesh.dim;
    cell_facets = dofmap.mesh.get_connectivity(dim, dim-1);

    % Prepare function for computing local boundary facets on cell c
    function local_facets = compute_local_boundary_facets(c)
        local_facets = full(facet_markers(cell_facets(:, c)) == [facet_tag{:}]);
    end

    func = @compute_local_boundary_facets;
end
