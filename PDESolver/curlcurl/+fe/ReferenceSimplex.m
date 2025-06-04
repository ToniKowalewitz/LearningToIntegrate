classdef ReferenceSimplex < handle


    properties (GetAccess = public, SetAccess = private)
         dim            % topological dimension of the cell
         vertex_coords  % coordinates of vertices
         normals        % outer unit normals of facets
         edge_tangents  % unit tangents of edges
         cell_volume    % volume of cell
         facet_area     % volume of facets
         edge_length    % length of edges
    end


    methods (Access = public)

        function obj = ReferenceSimplex(dim)
            % Constructor.

            obj.init(dim);
        end

        function conn = get_connectivity(obj, dim1, dim2)
            if dim1 == dim2
                conn = 1:obj.num_entities(dim1);
            elseif dim1 < dim2
                conn = invert_connectivity(obj.topology_{dim2+1, dim1+1});
            else
                conn = obj.topology_{dim1+1, dim2+1};
            end

            if isempty(conn)
                error('Connectivity %d-%d has not been computed', dim1, dim2);
            end
        end

        function num = num_entities(obj, dim)
            % Get number of entites of given dimension.

            if dim == obj.dim
                num = 1;
            elseif dim == 0
                num = obj.dim+1;
            else
                num = size(obj.topology_{dim+1, 0+1}, 2);
            end
        end

        function fun = get_facet_transform(obj, facet)
            % Return affine transform (function handle) from
            % reference simplex of lower dimension to given facet

            jac = obj.get_facet_transform_jacobian(facet);
            vec = obj.get_facet_transform_offset(facet);
            fun = @(t) jac*t + vec;
        end

        function jac = get_facet_transform_jacobian(obj, facet)
            % Return Jacobian matrix of affine transform from
            % reference simplex of lower dimension to given facet

            jac = obj.facet_transform_(:, 1:obj.dim-1, facet);
        end

        function vec = get_facet_transform_offset(obj, facet)
            % Return offset vector of affine transform from
            % reference simplex of lower dimension to given facet

            vec = obj.facet_transform_(:, obj.dim, facet);
        end

        function fun = get_edge_transform(obj, edge)
            % Return affine transform (function handle) from
            % reference interval to given edge

            jac = obj.get_edge_transform_jacobian(edge);
            vec = obj.get_edge_transform_offset(edge);
            fun = @(t) jac*t + vec;
        end

        function jac = get_edge_transform_jacobian(obj, edge)
            % Return Jacobian matrix of affine transform from
            % reference interval to given edge

            jac = obj.edge_transform_(:, 1, edge);
        end

        function vec = get_edge_transform_offset(obj, edge)
            % Return offset vector of affine transform from
            % reference interval to given edge

            vec = obj.edge_transform_(:, end, edge);
        end

        function lat = create_lattice_interior(obj, d, q)
            % Create lattice in interior of entities of dimension d
            % with q divisions in every direction; the result is in
            % normal Cartesian coordinates.

            % Get barycentric lattice on simplex of dimession d
            bar = interior_barycentric_lattice(d, q);
            bar = bar/q;  % scale to unit simplex

            % Get topology
            entities_to_vertex = obj.get_connectivity(d, 0);
            num_entities = size(entities_to_vertex, 2);

            % Iterate over entities and transform local barycentric to cartesian
            lat = zeros(size(bar, 1), obj.dim, num_entities);
            for entity = 1:num_entities
                entity_vertices = entities_to_vertex(:, entity);
                entity_coords = obj.vertex_coords(:, entity_vertices);
                lat(:, :, entity) = bar*entity_coords.';
            end
        end

    end


    methods (Access = private, Sealed = true)
        function init(obj, dim)  %#ok<*PROPLC>
            topology = cell(dim+1, dim+1);

            switch dim
            case 1
                vertex_coords = [
                    1, 0;
                ];

                normals = [
                     1, -1;
                ];

                edge_tangents = [  %#ok<NBRAK>
                     -1;
                ];

                facet_transform(:, :, 1) = [  %#ok<NBRAK>
                     1;
                ];
                facet_transform(:, :, 2) = [  %#ok<NBRAK>
                     0;
                ];

                edge_transform(:, :, 1) = [
                     1, 0;
                ];

                cell_volume = 1;

                facet_area = [1; 1;];

                edge_length = [1;];  %#ok<NBRAK>

                cell_vertex = uint32(1:2).';

                topology{1+1, 0+1} = cell_vertex;

            case 2
                vertex_coords = [
                    1, 0, 0;
                    0, 1, 0;
                ];

                normals = [
                    sqrt(1/2),  0, -1;
                    sqrt(1/2), -1,  0;
                ];

                edge_tangents = [
                     -sqrt(1/2), -1,  0;
                      sqrt(1/2),  0, -1;
                ];

                facet_transform(:, :, 1) = [
                     1,  0;
                    -1,  1;
                ];
                facet_transform(:, :, 2) = [
                     1,  0;
                     0,  0;
                ];
                facet_transform(:, :, 3) = [
                     0,  0;
                     1,  0;
                ];

                edge_transform = facet_transform;

                cell_volume = 1/2;

                facet_area = [sqrt(2); 1; 1;];

                edge_length = facet_area;

                cell_vertex = uint32(1:3).';

                edge_vertex = uint32([
                    1, 1, 2;
                    2, 3, 3;
                ]);

                topology{2+1, 0+1} = cell_vertex;
                topology{1+1, 0+1} = edge_vertex;

            case 3
                vertex_coords = [
                    1, 0, 0, 0;
                    0, 1, 0, 0;
                    0, 0, 1, 0;
                ];

                normals = [
                    sqrt(1/3),  0,  0, -1;
                    sqrt(1/3),  0, -1,  0;
                    sqrt(1/3), -1,  0,  0;
                ];

                edge_tangents = [
                     -sqrt(1/2), -sqrt(1/2), -1,          0,  0,  0;
                      sqrt(1/2),          0,  0, -sqrt(1/2), -1,  0;
                              0,  sqrt(1/2),  0,  sqrt(1/2),  0, -1;
                ];

                facet_transform(:, :, 1) = [
                     1,  0,  0;
                     0,  1,  0;
                    -1, -1,  1;
                ];
                facet_transform(:, :, 2) = [
                     1,  0,  0;
                     0,  1,  0;
                     0,  0,  0;
                ];
                facet_transform(:, :, 3) = [
                     1,  0,  0;
                     0,  0,  0;
                     0,  1,  0;
                ];
                facet_transform(:, :, 4) = [
                     0,  0,  0;
                     1,  0,  0;
                     0,  1,  0;
                ];

                edge_transform(:, :, 1) = [
                     1,  0;
                    -1,  1;
                     0,  0;
                ];
                edge_transform(:, :, 2) = [
                     1,  0;
                     0,  0;
                    -1,  1;
                ];
                edge_transform(:, :, 3) = [
                     1,  0;
                     0,  0;
                     0,  0;
                ];
                edge_transform(:, :, 4) = [
                     0,  0;
                     1,  0;
                    -1,  1;
                ];
                edge_transform(:, :, 5) = [
                     0,  0;
                     1,  0;
                     0,  0;
                ];
                edge_transform(:, :, 6) = [
                     0,  0;
                     0,  0;
                     1,  0;
                ];

                cell_volume = 1/6;

                facet_area = [sqrt(3)/2; 1/2; 1/2; 1/2;];

                edge_length = [sqrt(2); sqrt(2); 1; sqrt(2); 1; 1;];

                cell_vertex = uint32(1:4).';

                face_vertex = uint32([
                    1, 1, 1, 2;
                    2, 2, 3, 3;
                    3, 4, 4, 4;
                ]);

                edge_vertex = uint32([
                    1, 1, 1, 2, 2, 3;
                    2, 3, 4, 3, 4, 4;
                ]);

                face_edge = uint32([
                    1, 1, 2, 4;
                    2, 3, 3, 5;
                    4, 5, 6, 6;
                ]);

                topology{3+1, 0+1} = cell_vertex;
                topology{2+1, 0+1} = face_vertex;
                topology{1+1, 0+1} = edge_vertex;
                topology{2+1, 1+1} = face_edge;

            otherwise
                error('Dimension %d not implemented', dim);
            end

            obj.dim = dim;
            obj.vertex_coords = vertex_coords;
            obj.normals = normals;
            obj.edge_tangents = edge_tangents;
            obj.cell_volume = cell_volume;
            obj.facet_area = facet_area;
            obj.edge_length = edge_length;
            obj.topology_ = topology;
            obj.facet_transform_ = facet_transform;
            obj.edge_transform_ = edge_transform;
        end
    end


    properties (Access = private)
        topology_         % cell array with connectivities
        facet_transform_  % simplex-to-facet affine-transform
        edge_transform_  % simplex-to-facet affine-transform
    end


end


function lat = interior_barycentric_lattice(dim, q)
    % Create interior simplex lattice in Z^dim

    simplex = fe.barycentric_lattice(dim, q);
    simplex_interior = simplex(all(simplex>0, 2), :);
    lat = simplex_interior;
end


function inverse = invert_connectivity(conn)
    n2 = max(conn(:));
    n1 = numel(conn) / n2;
    inverse = zeros(n1, n2, 'uint32');
    for j=1:n2
        inverse(:, j) = find(any(conn==j, 1));
    end
end
