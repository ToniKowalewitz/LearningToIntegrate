classdef TestAssembleLaplaceSensitivity < matlab.unittest.TestCase

    properties (TestParameter)
        factory = {
            % FIXME: Can we have shared fixture with finite elements?
            %        P5 on tetrahedron is expensive to create.
            {@meshing.generate_half_disc_mesh, 3        , @lag, 1, 0};
            {@meshing.generate_half_disc_mesh, 2        , @lag, 2, 0};
            {@meshing.generate_half_disc_mesh, 1        , @lag, 3, 0};
            {@meshing.generate_unit_cube_mesh, [2, 2   ], @lag, 1, 0};
            {@meshing.generate_unit_cube_mesh, [1, 1   ], @lag, 2, 0};
            {@meshing.generate_unit_cube_mesh, [2, 2, 2], @lag, 1, 0};
            {@meshing.generate_unit_cube_mesh, [1, 1, 1], @lag, 2, 0};
            {@meshing.generate_half_disc_mesh, 3        , @lag, 1, 432};
            {@meshing.generate_half_disc_mesh, 2        , @lag, 2, 432};
            {@meshing.generate_half_disc_mesh, 1        , @lag, 3, 432};
            {@meshing.generate_unit_cube_mesh, [2, 2   ], @lag, 1, 432};
            {@meshing.generate_unit_cube_mesh, [1, 1   ], @lag, 2, 432};
            {@meshing.generate_unit_cube_mesh, [2, 2, 2], @lag, 1, 432};
            {@meshing.generate_unit_cube_mesh, [1, 1, 1], @lag, 2, 432};
        };
    end

    methods (Test)

        function test_assemble_laplace_sensitivity(t, factory)
            % Extract parameters
            [mesh_factory, mesh_resolution, element_factory, element_degree, lambda] = factory{:};

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

            % Build dofmap
            dofmap = assembling.build_dofmap(mesh, element);

            % Compare with Laplacian cell-by-cell
            num_cells = mesh.num_entities(mesh.dim);
            sigma = sparse(num_cells, 1);
            I = eye(dofmap.dim);
            M1 = repmat(I, 1, dofmap.dim);
            M2 = reshape(repmat(I, dofmap.dim, 1), dofmap.dim, []);
            DA = assembling.assemble_laplace_sensitivity(dofmap, I, M1, M2, lambda);
            for c = 1:num_cells
                sigma(c) = 1;  %#ok<SPRIX>
                A = assembling.assemble_laplace(dofmap, sigma, -lambda*sigma);
                t.verifyEqual(norm(A(:)-DA*sigma, inf), 0, 'AbsTol', 1e-12);
                sigma(c) = 0;  %#ok<SPRIX>
            end
        end

    end

end


function varargout = lag(varargin)
    [varargout{1:nargout}] = fe.create_lagrange_element(varargin{:});
end
