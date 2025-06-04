classdef TestDofmap < matlab.unittest.TestCase

    properties (TestParameter)
        factory = {
            % FIXME: Can we have shared fixture with finite elements?
            %        P5 on tetrahedron is expensive to create.
            {@meshing.generate_half_disc_mesh, 6        , @fe.create_lagrange_element,       1};
            {@meshing.generate_half_disc_mesh, 6        , @fe.create_lagrange_element,       2};
            {@meshing.generate_half_disc_mesh, 6        , @fe.create_lagrange_element,       3};
            {@meshing.generate_half_disc_mesh, 6        , @fe.create_lagrange_element,       4};
            {@meshing.generate_unit_cube_mesh, [2, 2   ], @fe.create_lagrange_element,       1};
            {@meshing.generate_unit_cube_mesh, [2, 2   ], @fe.create_lagrange_element,       2};
            {@meshing.generate_unit_cube_mesh, [2, 2   ], @fe.create_lagrange_element,       3};
            {@meshing.generate_unit_cube_mesh, [2, 2   ], @fe.create_lagrange_element,       4};
            {@meshing.generate_unit_cube_mesh, [2, 2, 2], @fe.create_lagrange_element,       1};
            {@meshing.generate_unit_cube_mesh, [2, 2, 2], @fe.create_lagrange_element,       2};
            {@meshing.generate_unit_cube_mesh, [2, 2, 2], @fe.create_lagrange_element,       3};
            {@meshing.generate_unit_cube_mesh, [2, 2, 2], @fe.create_lagrange_element,       4};
            {@meshing.generate_unit_cube_mesh, [2, 2, 2], @fe.create_lagrange_element,       5};
            {@meshing.generate_half_disc_mesh, 6        , @fe.create_nedelec_element,        1};
            {@meshing.generate_half_disc_mesh, 6        , @fe.create_nedelec_element,        2};
            {@meshing.generate_unit_cube_mesh, [2, 2   ], @fe.create_nedelec_element,        1};
            {@meshing.generate_unit_cube_mesh, [2, 2   ], @fe.create_nedelec_element,        2};
            {@meshing.generate_unit_cube_mesh, [2, 2, 2], @fe.create_nedelec_element,        1};
            {@meshing.generate_unit_cube_mesh, [2, 2, 2], @fe.create_nedelec_element,        2};
            {@meshing.generate_half_disc_mesh, 6        , @fe.create_raviart_thomas_element, 1};
            {@meshing.generate_unit_cube_mesh, [2, 2   ], @fe.create_raviart_thomas_element, 1};
            {@meshing.generate_unit_cube_mesh, [2, 2, 2], @fe.create_raviart_thomas_element, 1};
        };
    end

    methods (Test)

        function test_build_dofmap(t, factory)
            % Extract parameters
            [mesh_factory, mesh_resolution, element_factory, element_degree] = factory{:};

            % Build mesh
            mesh = mesh_factory(mesh_resolution);

            % Create element
            element = element_factory(mesh.dim, element_degree);

            % Supress warning about unneeded connectivities computed and stored
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            t.applyFixture(SuppressedWarningsFixture('Mesh:ExtraConnStored'));

            % Build connectivities needed to build dofmap and clear unneeded
            for d = element.get_dof_entity_dims()
                mesh.compute_connectivity(mesh.dim, d);
                mesh.clear_connectivity(d, 0);
            end

            % Test that dofmap building works
            dofmap = assembling.build_dofmap(mesh, element);

            % Test basic properties of dofmap
            t.verifyEqual(size(dofmap.cell_dofs, 1),  element.fe_space_dim);
            t.verifyEqual(size(dofmap.cell_dofs, 2),  mesh.num_entities(mesh.dim));
            t.verifyEqual(unique(dofmap.cell_dofs(:)), (uint32(1):uint32(dofmap.dim)).');

            % Verify Lagrange dofs coords
            verify_lagrange_dofs_coordinates(t, dofmap);
        end

    end

end


function verify_lagrange_dofs_coordinates(t, dofmap)
    % Check that the function fails on non-Lagrange
    if ~strcmp(dofmap.element.family, 'Lagrange')
        t.verifyError(@() assembling.get_lagrange_dofs_coordinates(dofmap), ...
                      ?MException);
        return;
    end

    % Call the tested function
    coords1 = assembling.get_lagrange_dofs_coordinates(dofmap);

    % Compute dofs coords manually
    coords2 = compute_lagrange_dofs_coords(dofmap);

    % Compare
    t.verifyEqual(coords1, coords2, 'AbsTol', 1e-14, 'RelTol', 1e-14);
end


function coords = compute_lagrange_dofs_coords(dofmap)
    dim = dofmap.mesh.dim;
    num_cells = dofmap.mesh.num_entities(dim);
    vertex_coords = dofmap.mesh.vertex_coords;
    cells = dofmap.mesh.cells;
    cell_dofs = dofmap.cell_dofs;
    local_dofs_coords = dofmap.element.evaluate_dual_basis(@(x) x);

    coords = zeros(dim, dofmap.dim);
    for c = 1:num_cells
        % Get affine mapping for current cell
        B = vertex_coords(:, cells(1:dim, c)) - vertex_coords(:, cells(dim+1, c));
        b = vertex_coords(:, cells(dim+1, c));

        % Transform local dofs coords to global and store
        coords(:, cell_dofs(:, c)) = B*local_dofs_coords + b;
    end
end
