classdef TestEvalFun < matlab.unittest.TestCase

    properties (TestParameter)
        factory_scalar = {
            {@meshing.generate_half_disc_mesh, 3        , @fe.create_lagrange_element,       1};
            {@meshing.generate_half_disc_mesh, 2        , @fe.create_lagrange_element,       2};
            {@meshing.generate_half_disc_mesh, 1        , @fe.create_lagrange_element,       3};
            {@meshing.generate_unit_cube_mesh, [2, 2   ], @fe.create_lagrange_element,       1};
            {@meshing.generate_unit_cube_mesh, [1, 1   ], @fe.create_lagrange_element,       2};
            {@meshing.generate_unit_cube_mesh, [2, 2, 2], @fe.create_lagrange_element,       1};
            {@meshing.generate_unit_cube_mesh, [1, 1, 1], @fe.create_lagrange_element,       2};
            {@meshing.generate_half_disc_mesh, 4        , @fe.create_discontinuous_lagrange_element, 0};
            {@meshing.generate_half_disc_mesh, 3        , @fe.create_discontinuous_lagrange_element, 1};
            {@meshing.generate_half_disc_mesh, 2        , @fe.create_discontinuous_lagrange_element, 2};
            {@meshing.generate_half_disc_mesh, 1        , @fe.create_discontinuous_lagrange_element, 3};
            {@meshing.generate_unit_cube_mesh, [4, 4   ], @fe.create_discontinuous_lagrange_element, 0};
            {@meshing.generate_unit_cube_mesh, [2, 2   ], @fe.create_discontinuous_lagrange_element, 1};
            {@meshing.generate_unit_cube_mesh, [1, 1   ], @fe.create_discontinuous_lagrange_element, 2};
            {@meshing.generate_unit_cube_mesh, [4, 4, 4], @fe.create_discontinuous_lagrange_element, 0};
            {@meshing.generate_unit_cube_mesh, [2, 2, 2], @fe.create_discontinuous_lagrange_element, 1};
            {@meshing.generate_unit_cube_mesh, [1, 1, 1], @fe.create_discontinuous_lagrange_element, 2};
        };
        factory_vector_valued_dim2 = {
            {@meshing.generate_half_disc_mesh, 3        , @fe.create_nedelec_element,        1};
            {@meshing.generate_half_disc_mesh, 2        , @fe.create_nedelec_element,        2};
            {@meshing.generate_unit_cube_mesh, [2, 2   ], @fe.create_nedelec_element,        1};
            {@meshing.generate_unit_cube_mesh, [1, 1   ], @fe.create_nedelec_element,        2};
            {@meshing.generate_half_disc_mesh, 3        , @fe.create_raviart_thomas_element, 1};
            {@meshing.generate_unit_cube_mesh, [2, 2   ], @fe.create_raviart_thomas_element, 1};
        };
        factory_vector_valued_dim3 = {
            {@meshing.generate_unit_cube_mesh, [2, 2, 2], @fe.create_nedelec_element,        1};
            {@meshing.generate_unit_cube_mesh, [1, 1, 1], @fe.create_nedelec_element,        2};
            {@meshing.generate_unit_cube_mesh, [2, 2, 2], @fe.create_raviart_thomas_element, 1};
        };
        func_scalar = {
            @(x, c) x(1);
            @(x, c) x(1)^2;
            @(x, c) sin(x(1));
            @(x, c) cos(x(1));
            @(x, c) exp(x(1));
            @(x, c) x(2);
            @(x, c) x(2)^2;
            @(x, c) sin(x(2));
            @(x, c) cos(x(2));
            @(x, c) exp(x(2));
            @(x, c) norm(x);
            @(x, c) norm(x)^2;
            @(x, c) sin(norm(x));
            @(x, c) cos(norm(x));
            @(x, c) exp(norm(x));
        }
        func_vector_valued_dim2 = {
            @(x, c) [x(1), x(2)];
            @(x, c) [x(2), -x(1)];
            @(x, c) [sin(x(1)), cos(x(2))];
            @(x, c) [cos(x(1)), sin(x(2))];
            @(x, c) [-cos(x(2)), sin(x(1))];
        }
        func_vector_valued_dim3 = {
            @(x, c) [x(1), x(2), x(3)];
            @(x, c) [x(2), -x(1), x(3)];
            @(x, c) [x(3), -x(1), x(2)];
            @(x, c) [sin(x(1)), cos(x(2)), exp(x(3))];
            @(x, c) [cos(x(1)), sin(x(2)), exp(x(3))];
            @(x, c) [-cos(x(2)), sin(x(1)), exp(x(3))];
        }
    end

    methods (Test)

        function test_scalar(t, factory_scalar, func_scalar)
            test_projection(t, factory_scalar, func_scalar);
        end

        function test_vector_valued_dim2(t, factory_vector_valued_dim2, func_vector_valued_dim2)
            test_projection(t, factory_vector_valued_dim2, func_vector_valued_dim2);
        end

        function test_vector_valued_dim3(t, factory_vector_valued_dim3, func_vector_valued_dim3)
            test_projection(t, factory_vector_valued_dim3, func_vector_valued_dim3);
        end

    end

    methods

        function test_projection(t, factory, func)
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

            % Build dofmap
            dofmap = assembling.build_dofmap(mesh, element);

            % Interpolate function to FE space
            vec = assembling.interpolate(dofmap, func);
            u = assembling.create_fe_eval_function(dofmap, vec);
            u2 = create_fe_eval_function_2(dofmap, vec);

            % We should have interpolated something non-trivial
            t.verifyTrue(norm(vec) > 0.1);

            % Check that interpolation is projection (evaluation on known cell)
            vec2 = assembling.interpolate(dofmap, u);
            t.verifyEqual(norm(vec-vec2), 0, 'AbsTol', 1e-13);

            % Evaluation on unknown cell has to fail before geometric queries are initialized
            t.verifyError(@() assembling.interpolate(dofmap, @(x, c) u(x)), ?MException);
            t.verifyError(@() assembling.interpolate(dofmap, @(x, c) u2(x)), ?MException);

            mesh.init_geometric_queries();

            % Check that interpolation is projection (evaluation on unknown cell)
            vec2 = assembling.interpolate(dofmap, @(x, c) u(x));
            t.verifyEqual(norm(vec-vec2), 0, 'AbsTol', 1e-13);

            % Check that interpolation is projection (using assemble_point_sources())
            vec2 = assembling.interpolate(dofmap, @(x, c) u2(x));
            t.verifyEqual(norm(vec-vec2), 0, 'AbsTol', 1e-13);
        end

    end

end


function func = create_fe_eval_function_2(dofmap, vec)
    % Wrap assembling.assemble_point_sources() to compare
    % it with assembling.create_fe_eval_function()
    function val = u(x)
        Q = assembling.assemble_point_sources(dofmap, x);
        val = vec.'*Q;
    end
    func = @u;
end
