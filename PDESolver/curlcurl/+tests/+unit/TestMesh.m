classdef TestMesh < matlab.unittest.TestCase

    properties (TestParameter)
        factory = {
            @() meshing.generate_half_disc_mesh(1                      );
            @() meshing.generate_half_disc_mesh(2                      );
            @() meshing.generate_half_disc_mesh(4                      );
            @() meshing.generate_half_disc_mesh(7                      );
            @() meshing.generate_unit_cube_mesh([1, 1   ]              );
            @() meshing.generate_unit_cube_mesh([2, 2   ]              );
            @() meshing.generate_unit_cube_mesh([5, 4   ]              );
            @() meshing.generate_unit_cube_mesh([8, 7   ]              );
            @() meshing.generate_unit_cube_mesh([1, 1, 1]              );
            @() meshing.generate_unit_cube_mesh([2, 4, 7]              );
            @() meshing.generate_rectangle_mesh(0, 1, 0, 2, 'Hmin', 1  );
            @() meshing.generate_rectangle_mesh(0, 1, 0, 2, 'Hmax', 0.2);
        };
    end

    methods (Test)

        function test_topology_ok(t, factory)
            mesh = factory();
            verify_mesh_topology(t, mesh);
        end

        function test_facet_cell_connectivity(t, factory)
            % Build mesh
            mesh = factory();

            % Supress warning about unneeded connectivities computed and stored
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            t.applyFixture(SuppressedWarningsFixture('Mesh:ExtraConnStored'));

            % Compute cell-facet and facet-cell connectivity
            mesh.compute_connectivity(mesh.dim, mesh.dim-1);
            mesh.compute_connectivity(mesh.dim-1, mesh.dim);
            c2f = mesh.get_connectivity(mesh.dim, mesh.dim-1);
            f2c = mesh.get_connectivity(mesh.dim-1, mesh.dim);

            % Compute connectivity needed to assemble facet integrals
            mesh.compute_boundary_facets();
            bf = mesh.get_boundary_facets_indices();

            % Test that boundary facets indices are unique
            t.verifyEqual(numel(bf), numel(unique(bf)));

            num_facets = mesh.num_entities(mesh.dim-1);
            num_boundary_facets = numel(bf);

            % Find column indices (of f2c) sorting f2c(:)
            [~, j] = sort(f2c(:));
            [~, sorted_facets] = ind2sub([2, num_facets], j);
            sorted_facets = uint32(sorted_facets);

            % Check that zero elements of f2c map to boundary
            t.verifyEqual(sorted_facets(1:num_boundary_facets), sort(bf));

            % Check that positive elements of f2c map to c2f
            t.verifyEqual(sorted_facets(num_boundary_facets+1:end), c2f(:));

        end

        function test_cell_incenters(t, factory)
            % Build mesh
            mesh = factory();

            % Only dimension 2 is implemented
            if mesh.dim ~= 2
                t.verifyError(@() mesh.get_cell_incenters(), ?MException);
                return
            end

            % Compute incenters
            incenters = mesh.get_cell_incenters();

            num_cells = mesh.num_entities(mesh.dim);
            e2v_ref = fe.ReferenceSimplex(mesh.dim).get_connectivity(1, 0);

            % Loop over cells
            for c = 1:num_cells

                % Loop over edges
                for e = 1:3
                    edge_verts = mesh.vertex_coords(:, mesh.cells(e2v_ref(:, e), c));

                    % Compute distance of incenter to edge e
                    dist(e) = distance_point_segment(incenters(:, c), edge_verts);  %#ok<AGROW>
                end

                % Check that distance of incenter to all edges is the same
                t.verifyEqual(dist(2), dist(1), 'RelTol', 1e-13);
                t.verifyEqual(dist(3), dist(1), 'RelTol', 1e-13);

            end

        end

        function test_cell_volumes(t, factory)
            % Build mesh
            mesh = factory();

            vols = mesh.get_cell_volumes();

            % Verify sign
            t.verifyTrue(all(vols > 0));

            % Test value cell-by-cell
            num_cells = mesh.num_entities(mesh.dim);
            coords = mesh.vertex_coords;
            cells = mesh.cells;
            dim = mesh.dim;
            for c = 1:num_cells
                % https://en.wikipedia.org/w/index.php?title=Simplex&oldid=936472880#Volume
                vol = abs(det(coords(:, cells(1:dim, c)) - coords(:, cells(dim+1, c))) / factorial(dim));
                t.verifyEqual(vols(c), vol, 'RelTol', 1e-15);
            end
        end

    end

end


function dist = distance_point_segment(pt, seg)
    % Return distance of point to segment
    assert(iscolumn(pt));
    a = vecnorm(pt - seg(:, 1));
    b = vecnorm(pt - seg(:, 2));
    c = vecnorm(seg(:, 2) - seg(:, 1));
    s = 0.5 * (a + b + c);
    dist = 2 * sqrt(s*(s-a)*(s-b)*(s-c)) / c;
end
