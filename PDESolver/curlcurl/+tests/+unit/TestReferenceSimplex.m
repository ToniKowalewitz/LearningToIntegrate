classdef TestReferenceSimplex < matlab.unittest.TestCase

    properties (TestParameter)
        dim = {1, 2, 3};
    end

    methods (Test)

        function test_inverse_connectivity(t, dim)
            s = fe.ReferenceSimplex(dim);

            % Test consistency of all d-0 and 0-d connectivities
            for d = 1:dim
                conn = s.get_connectivity(d, 0);
                inverse = s.get_connectivity(0, d);
                verify_inverse(t, conn, inverse);
                verify_inverse(t, inverse, conn);
            end

            % For tetrahedra we have additionally implemented 2-1 connectivity
            if dim == 3
                conn = s.get_connectivity(2, 1);
                inverse = s.get_connectivity(1, 2);
                verify_inverse(t, conn, inverse);
                verify_inverse(t, inverse, conn);
            end

            % Connectivities we did not implement
            if dim == 2
                t.verifyError(@() s.get_connectivity(2, 1), ?MException);
                t.verifyError(@() s.get_connectivity(1, 2), ?MException);
            elseif dim == 3
                t.verifyError(@() s.get_connectivity(3, 1), ?MException);
                t.verifyError(@() s.get_connectivity(3, 2), ?MException);
                t.verifyError(@() s.get_connectivity(1, 3), ?MException);
                t.verifyError(@() s.get_connectivity(2, 3), ?MException);
            end
        end


        function test_composition_consistency(t, dim)
            s = fe.ReferenceSimplex(dim);

            % Test consistency of 2-1-0 with 2-0 on tetrahedron
            if dim == 3
                c21 = s.get_connectivity(2, 1);
                c10 = s.get_connectivity(1, 0);
                c20 = s.get_connectivity(2, 0);
                c210 = unique(reshape(c10(:, c21(:, :)), [], size(c20, 2)), 'rows');
                t.verifyEqual(c210, c20);
            end

            % Candidates for composed connectivities we did not implement
            if dim == 2
                t.verifyError(@() s.get_connectivity(2, 1), ?MException);
            elseif dim == 3
                t.verifyError(@() s.get_connectivity(3, 1), ?MException);
                t.verifyError(@() s.get_connectivity(3, 2), ?MException);
            end
        end


        function test_cell_volume(t, dim)
            s = fe.ReferenceSimplex(dim);
            coords = s.vertex_coords;

            % https://en.wikipedia.org/w/index.php?title=Simplex&oldid=936472880#Volume
            vol = abs(det(coords(:, 1:dim) - coords(:, dim+1)) / factorial(dim));

            t.verifyEqual(s.cell_volume, vol, 'RelTol', 1e-15);
        end


        function test_facet_area(t, dim)
            s = fe.ReferenceSimplex(dim);

            % Area of reference (dim-1)-simplex
            scaling_factor = factorial(dim-1);

            num_facets = s.num_entities(dim-1);
            for f = 1:num_facets
                area1 = prod(svd(s.get_facet_transform_jacobian(f))) / scaling_factor;
                area2 = s.facet_area(f);
                t.verifyEqual(area1, area2, 'RelTol', 1e-12);
            end
        end

        function test_edge_length(t, dim)
            s = fe.ReferenceSimplex(dim);

            num_edges = s.num_entities(1);
            for e = 1:num_edges
                et = s.get_edge_transform(e);
                t.verifyEqual(norm(et(1)-et(0)), s.edge_length(e));
            end
        end

        function test_facet_transform_offset(t, dim)
            s = fe.ReferenceSimplex(dim);
            num_facets = s.num_entities(dim-1);

            for f = 1
                x0 = s.get_facet_transform_offset(f);
                t.verifyEqual(x0, [zeros(dim-1, 1); 1]);
            end
            for f = 2:num_facets
                x0 = s.get_facet_transform_offset(f);
                t.verifyEqual(x0, [zeros(dim-1, 1); 0]);
            end
        end

        function test_facet_transform_jacobian(t, dim)  %#ok<INUSD>
            % FIXME: How to test it?
        end

        function test_edge_transform_offset(t, dim)
            s = fe.ReferenceSimplex(dim);

            num_edges = s.num_entities(1);
            edge_vertices = s.get_connectivity(1, 0);

            for e = 1:num_edges
                max_vertex_index = max(edge_vertices(:, e));
                t.verifyEqual(...
                    s.vertex_coords(:, max_vertex_index), ...
                    s.get_edge_transform_offset(e), ...
                    'AbsTol', 1e-12);
            end
        end

        function test_edge_transform_jacobian(t, dim)
            s = fe.ReferenceSimplex(dim);

            % Corner case
            if dim == 1
                t.verifyEqual(s.get_edge_transform_jacobian(1), 1);
                return
            end

            num_edges = s.num_entities(1);
            facet_edges = s.get_connectivity(dim-1, 1);

            for e = 1:num_edges
                facets_contain_this_edge = any(facet_edges == e, 1);

                % Test that image of edge jacobian is orthogonal to
                % adjacent facet normals
                t.verifyEqual(...
                    s.get_edge_transform_jacobian(e).'*s.normals ~= 0, ...
                    ~facets_contain_this_edge);

                % FIXME: How to test proper scaling of edge jacobian?
            end
        end

        function test_normals(t, dim)
            s = fe.ReferenceSimplex(dim);

            num_facets = s.num_entities(dim-1);
            facet_vertices = s.get_connectivity(dim-1, 0);

            for f = 1:num_facets
                n = s.normals(:, f);

                % Test norm
                t.verifyEqual(norm(n), 1, 'RelTol', 1e-12);

                % Test orthogonality to facet
                jac = s.get_facet_transform_jacobian(f);
                t.verifyEqual(jac.'*n, zeros(dim-1, 1), 'AbsTol', 1e-12);

                % Test outer direction
                facet_midpoint = mean(s.vertex_coords(:, facet_vertices(:, f)), 2);
                cell_midpoint = mean(s.vertex_coords, 2);
                cosine = (facet_midpoint-cell_midpoint).'*n;
                t.verifyTrue(cosine > 0.1);
            end
        end

        function test_edge_tangents(t, dim)
            s = fe.ReferenceSimplex(dim);

            num_edges = s.num_entities(1);
            edge_vertices = s.get_connectivity(1, 0);

            for e = 1:num_edges
                T = s.edge_tangents(:, e);

                % Test norm
                t.verifyEqual(norm(T), 1, 'RelTol', 1e-12);

                % Test direction
                vertex_index_hi = max(edge_vertices(:, e));
                vertex_index_lo = min(edge_vertices(:, e));
                direction = s.vertex_coords(:, vertex_index_hi) ...
                          - s.vertex_coords(:, vertex_index_lo);
                t.verifyEqual(T, direction/norm(direction), 'AbsTol', 1e-12);
            end
        end

    end

end


function verify_inverse(t, conn1, conn2)
    % Verify that composition of conn1 and conn2 contains identity mapping
    for j = 1:size(conn2, 2)
        col = intersect_cols(conn1(:, conn2(:, j)));
        t.verifyEqual(intersect(col, j), uint32(j));
    end
end


function result = intersect_cols(arr)
    result = arr(:, 1);
    for j = 2:size(arr, 2)
        result = intersect(arr(:, j), result);
    end
end
