classdef Mesh < handle
    % Mesh class describing a simplicial mesh.

    properties (SetAccess = immutable)
        dim    % scalar, denoting spatial dimension
        cells  % matrix [dim+1 x k] of cell-vertex connectivity
    end

    properties (SetAccess = private)
        vertex_coords  % matrix [dim x n] of coordinates of vertices
    end

    properties (Access = private)
        topology_         % cell array with connectivity arrays
        num_entities_     % cell array with number of entities
        boundary_facets_  % cell array with list of boundary facets
        tri_              % cell array with triangulation instance
    end

    methods (Access = public)

        function obj = Mesh(dim, vertex_coords, cells)
            % Constructor.

            % Check input
            validateattributes(dim, {'numeric'}, ...
                               {'scalar', 'integer', 'positive'}, ...
                               mfilename, 'dim', 1);
            validateattributes(vertex_coords, {'double'}, ...
                               {'nrows', dim}, ...
                               mfilename, 'vertex_coords', 2);
            num_vertices = size(vertex_coords, 2);
            validateattributes(cells, {'uint32'}, ...
                               {'nrows', dim+1, '>=', 1, '<=', num_vertices}, ...
                               mfilename, 'cells', 3);

            % Check there are no unconnected vertices
            unconnected_vertices = setdiff(1:num_vertices, cells);
            if ~isempty(unconnected_vertices)
                warning('Mesh:UnconnectedVertices', ...
                        'Vertices %s not contained in any cell', ...
                        mat2str(unconnected_vertices));
            end

            % Ensure increasing order of vertices on every cell
            % NB: This is absolutely crucial to do (before anything else
            %     is computed from cells) for validity of dofmaps!
            cells = sort(cells);

            % Stored passed data
            obj.dim = dim;
            obj.vertex_coords = vertex_coords;
            obj.cells = cells;

            % Init private data structures
            % NB: Indices offset by one due to matlab's one-based indexing
            %     dim == 0 -> points  -> idx 1
            %     dim == 1 -> edges   -> idx 2
            %     dim == 2 -> faces   -> idx 3
            %     dim == 3 -> volumes -> idx 4
            obj.topology_ = cell(dim+1, dim+1);
            obj.topology_{dim+1, 1} = cells;
            obj.num_entities_{0  +1} = size(vertex_coords, 2);
            obj.num_entities_{dim+1} = size(cells, 2);
            obj.boundary_facets_ = cell(0, 0);
            obj.tri_ = cell(0, 0);
        end

        function conn = get_connectivity(obj, dim1, dim2)
            % Return connectivity dim1-dim2 assuming it has been computed.
            %
            % SYNTAX
            %   conn = get_connectivity(obj, dim1, dim2)

            if dim1 == dim2
                % Return identity mapping
                conn = 1:obj.num_entities(dim1);
            else
                conn = obj.topology_{dim1+1, dim2+1};
                if isempty(conn)
                    error('Mesh:ConnNotComputed', ...
                          'Connectivity %d-%d has not been computed', dim1, dim2);
                end
            end
        end

        function num = num_entities(obj, dim)
            % Return number of entities of given dimension.

            num = obj.num_entities_{dim+1};
            if isempty(num)
                error('Entities of dim %d have not been computed', dim);
            end
        end

        function compute_connectivity(obj, dim1, dim2)
            % Compute dim1-dim2 connectivity.
            %
            % The following cases are supported:
            %
            %  * dim-d        ... cell-(d-dimensional entities) connectivity
            %  * d-0          ... (d-dimensional entities)-vertex connectivity
            %  * d-d          ... identity
            %  * (dim-1)-dim  ... facet-cell connectity (contains zeros for
            %                     boundary facets)

            % Check if already computed
            if not(isempty(obj.topology_{dim1+1, dim2+1}))
                return
            end

            % Pick strategy for computing
            if dim1 == dim2
                % We don't store this special case - identity mapping
                return
            elseif (dim1 == obj.dim-1 && dim2 == obj.dim)
                % Compute facet-cell by inverting cell-facet
                obj.compute_connectivity(dim2, dim1);
                num_facets = obj.num_entities(dim1);
                obj.topology_{dim1+1, dim2+1} = ...
                    invert_cell_facet_connectivity(obj.topology_{dim2+1, dim1+1}, num_facets);
                return
            elseif (dim1 == obj.dim && dim2 == 0)
                % We always have cell-vertex connectivity in obj.cells
                obj.topology_{obj.dim+1, 1} = obj.cells;
                return
            elseif (dim1 == obj.dim && dim2 < obj.dim && dim2 > 0)
                % Compute as cell-dim2-0
                d = dim2;
            elseif (dim1 < obj.dim && dim1 > 0 && dim2 == 0)
                % Compute as cell-dim1-0
                d = dim1;
            else
                error('Computing %d %d connectivity not implemented', dim1, dim2);
            end

            % Compute and store
            warning('Mesh:ExtraConnStored', ...
                    'Computing and storing %d-%d-%d connectivity', obj.dim, d, 0);
            [obj.topology_{obj.dim+1, d+1}, obj.topology_{d+1, 1}] = ...
                obj.compute_connectivity_dim_d_0(d);

             % Store computed number of entities
            obj.num_entities_{d+1} = size(obj.topology_{d+1, 1}, 2);
        end

        function clear_connectivity(obj, dim1, dim2)
            % Clear dim1-dim2 connectivity.

            % Never delete cell-vertex connectivity
            if dim1 == obj.dim && dim2 == 0
                return
            end

            obj.topology_{dim1+1, dim2+1} = [];
        end

        function boundary_facets = get_boundary_facets(obj)
            % Get logical sparse array being true for boundary facets.
            %
            % Column index is cell index and row index is local facet index.

            if isempty(obj.boundary_facets_)
                error('Boundary facets have not been computed');
            end
            boundary_facets = obj.boundary_facets_{1, 1};
        end

        function boundary_facets = get_boundary_facets_indices(obj)
            % Get boundary (global) facet indices

            [local_facets, cells] = find(obj.get_boundary_facets());  %#ok<PROP>
            c2f = obj.get_connectivity(obj.dim, obj.dim-1);
            boundary_facets = c2f(sub2ind(size(c2f), local_facets, cells));  %#ok<PROP>
        end

        function markers = get_boundary_facet_markers(obj, markers, marker_value)
            % Build/modify facet markers with given value on boundary
            %
            % Result is a sparse column vector of length given by the
            % number of facets in the mesh. The resulting vector takes
            % the given value on boundary. One can provide an existing
            % marker vector to be modified or empty vector.

            % Compute global indices of boundary facets
            boundary_facets_indices = obj.get_boundary_facets_indices();

            num_facets = obj.num_entities(obj.dim-1);
            num_boundary_facets = numel(boundary_facets_indices);

            % Allocate result if needed
            if isempty(markers)
                markers = sparse([], [], [], num_facets, 1, num_boundary_facets);
            end
            assert(numel(markers) == num_facets);

            % Mark
            markers(boundary_facets_indices) = marker_value;
        end

        function compute_boundary_facets(obj)
            % Compute and store indices of boundary facets.
            %
            % The method exploits, that boundary facets are associated with
            % only one cell contrary to all interior facets.

            % Fetch data
            num_facets = obj.num_entities(obj.dim-1);
            cell_facet_connectivity = obj.get_connectivity(obj.dim, obj.dim-1);

            % Determine boundary facets by counting the number of cells
            % they appear in (1 = boundary, 2 = interior)
            boundary_facets = logical(accumarray(cell_facet_connectivity(:), ...
                                      1, [num_facets, 1], ...
                                      @(t) (sum(t)==1)*1, 0, true));

            % Switch to (local_facet, cell)-indexing
            boundary_facets = boundary_facets(cell_facet_connectivity);

            % Store the results
            obj.boundary_facets_{1, 1} = boundary_facets;
        end

        function clear_boundary_facets(obj)
            % Clear boundary facets.

            obj.boundary_facets_ = cell(0, 0);
        end

        function set_coordinates(obj, coords)
            % Set vertex coordinates

            assert(all(size(coords) == size(obj.vertex_coords)));

            % Clear internal state of mesh which depends on coordinates
            obj.clear_geometric_queries();

            obj.vertex_coords = coords;
        end

        function init_geometric_queries(obj)
            % Initialize internal data structure for geometric queries.
            %
            % Generates a second instance of the given mesh in form of the
            % MATLAB own triangluation data structure.
            % This allows for fast access to geometric queries (e.g.
            % location of points in mesh - see point_location method).

            if isempty(obj.tri_)
                % FIXME: cast to double and transpose cause copying?
                mesh_cells = double(obj.cells.');
                obj.tri_{1, 1} = triangulation(mesh_cells, obj.vertex_coords.');
            end
        end

        function clear_geometric_queries(obj)
            % Clear all internal data associated with geometric queries.

            obj.tri_ = cell(0, 0);
        end

        function varargout = point_location(obj, varargin)
            % Find cell (optionally barycentric coordinates) for given points.
            % This has same interface as matlab's triangulation.pointLocation().

            if isempty(obj.tri_)
                error('Geometry queries not initalized');
            end
            [varargout{1:nargout}] = obj.tri_{1, 1}.pointLocation(varargin{:});
        end

        function coords = get_cell_centroids(obj)
            % Provides coordinates of cell centroids.
            %
            % OUTPUT PARAMETER
            %   coords ... Matrix [obj.dim, size(obj.cells, 2)].

            num_cells = obj.num_entities(obj.dim);
            coords = zeros(obj.dim, num_cells);
            for c = 1:num_cells
                coords(:, c) = sum(obj.vertex_coords(:, obj.cells(:, c)), 2) / (obj.dim + 1);
            end
        end

        function coords = get_cell_incenters(obj)
            % Provides coordinates of cell incenters.
            %
            % OUTPUT PARAMETER
            %   coords ... Matrix [obj.dim, size(obj.cells, 2)].

            if obj.dim ~= 2
                error('Cell incenters not implemented for dim %d', obj.dim);
            end

            num_cells = obj.num_entities(obj.dim);
            coords = zeros(obj.dim, num_cells);
            for c = 1:num_cells
                verts = obj.vertex_coords(:, obj.cells(:, c));
                sides = [
                    verts(:, 2) - verts(:, 3), ...
                    verts(:, 3) - verts(:, 1), ...
                    verts(:, 1) - verts(:, 2), ...
                ];
                lengths = vecnorm(sides, 2, 1);
                % https://en.wikipedia.org/w/index.php?title=Incenter&oldid=917270547#Cartesian_coordinates
                coords(:, c) = verts * lengths.' / sum(lengths);
            end
        end

        function vol = get_cell_volumes(obj)
            % Provides the cell volumes.
            %
            % OUTPUT PARAMETER
            %   coords ... Matrix [obj.dim, size(obj.cells, 2)].

            dim = obj.dim;  %#ok<PROP>
            cells = obj.cells;  %#ok<PROP>
            coords = obj.vertex_coords;
            num_cells = size(obj.cells, 2);

            reference_simplex = fe.ReferenceSimplex(obj.dim);
            reference_vol = reference_simplex.cell_volume;
            assert(reference_vol > 0);

            vol = zeros(num_cells, 1);
            for c = 1:num_cells
                % https://en.wikipedia.org/w/index.php?title=Simplex&oldid=936472880#Volume
                vol(c) = reference_vol * ...
                    abs(det(coords(:, cells(1:dim, c)) - coords(:, cells(dim+1, c))));  %#ok<PROP>
            end
        end

    end

    methods (Access = private, Sealed = true)

        function [cell_ent_connect, ent_vertex_connect] = compute_connectivity_dim_d_0(obj, d)
            % Compute dim-d-0 global connectivity from provided local d-0 connectivity.
            %
            % The dim-to-d and d-to-0 connectivities are stored separately.
            %
            % cell   ... entity of dimension dim
            % ent    ... entity of dimension d
            % vertex ... entity of dimension 0

            % Get local d-0 connectivity of reference cell
            reference_simplex = fe.ReferenceSimplex(obj.dim);
            ent_vertex_connect_local = reference_simplex.get_connectivity(d, 0);

            % Fetch data
            cell_vertex_connect = obj.cells;
            vertices_per_ent = size(ent_vertex_connect_local, 1);
            ent_per_cell = size(ent_vertex_connect_local, 2);
            num_cells = obj.num_entities(obj.dim);

            % Compute ent-vertex connectivity cell-by-cell
            ent_vertex_connect = reshape(cell_vertex_connect(ent_vertex_connect_local(:), :), ...
                                    vertices_per_ent, ent_per_cell*num_cells).';

            % Gather cells togeter, pick unique,that gives us ent-numbering
            [ent_vertex_connect, ~, cell_ent_connect] = unique(ent_vertex_connect, 'rows');

            % Adapt data into desired shape
            cell_ent_connect = reshape(uint32(cell_ent_connect), ent_per_cell, num_cells);
            ent_vertex_connect = ent_vertex_connect.';
        end

    end

end


function f2c = invert_cell_facet_connectivity(c2f, num_facets)
    f2c = zeros(2, num_facets, 'uint32');
    for c = 1:size(c2f, 2)
        [~, i] = min(f2c(:, c2f(:, c)), [], 1);
        i = sub2ind(size(f2c), i.', c2f(:, c));
        f2c(i) = c;
    end
end
