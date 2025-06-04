classdef TestAssembleNeumann < matlab.unittest.TestCase

    properties (TestParameter)
        factory = {
            % FIXME: Can we have shared fixture with finite elements?
            %        P5 on tetrahedron is expensive to create.
            {@meshing.generate_half_disc_mesh, 3        , @fe.create_lagrange_element,       1};
            {@meshing.generate_half_disc_mesh, 2        , @fe.create_lagrange_element,       2};
            {@meshing.generate_half_disc_mesh, 1        , @fe.create_lagrange_element,       3};
            {@meshing.generate_unit_cube_mesh, [2, 2   ], @fe.create_lagrange_element,       1};
            {@meshing.generate_unit_cube_mesh, [1, 1   ], @fe.create_lagrange_element,       2};
            {@meshing.generate_unit_cube_mesh, [2, 2, 2], @fe.create_lagrange_element,       1};
            {@meshing.generate_unit_cube_mesh, [1, 1, 1], @fe.create_lagrange_element,       2};
        };

        factory_cubes = {
            % FIXME: Can we have shared fixture with finite elements?
            %        P5 on tetrahedron is expensive to create.
            {@meshing.generate_unit_cube_mesh, [2, 2   ], @fe.create_lagrange_element,       1};
            {@meshing.generate_unit_cube_mesh, [1, 1   ], @fe.create_lagrange_element,       2};
            {@meshing.generate_unit_cube_mesh, [2, 2, 2], @fe.create_lagrange_element,       1};
            {@meshing.generate_unit_cube_mesh, [1, 1, 1], @fe.create_lagrange_element,       2};
        };
    end

    methods (Test)

        function test_assemble_neumann(t, factory)
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

            % Compute connectivity needed to assemble facet integrals
            mesh.compute_connectivity(mesh.dim, mesh.dim-1);
            mesh.clear_connectivity(dofmap.mesh.dim-1, 0);
            num_facets = mesh.num_entities(mesh.dim-1);
            mesh.compute_boundary_facets();

            % Assemble Neumann source term
            g = @(x, n, c) 1;
            degree_rise = 0;
            marker_value = 42;
            markers = sparse([], [], [], num_facets, 1, 1);
            for f = 1:num_facets
                markers(f) = marker_value;  %#ok<SPRIX>
                b = assembling.assemble_neumann_source(dofmap, markers, marker_value, g, degree_rise);
                verify_exterior_facet_integral_(t, dofmap, markers, marker_value, b);
                markers(:) = 0;  %#ok<SPRIX>
            end
        end

        function test_assemble_neumann_normals(t, factory_cubes)
            % Extract parameters
            [mesh_factory, mesh_resolution, element_factory, element_degree] = factory_cubes{:};

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

            % Compute connectivity needed to assemble facet integrals and mark facets
            mesh.compute_connectivity(mesh.dim, mesh.dim-1);
            mesh.compute_boundary_facets();
            mesh.clear_connectivity(dofmap.mesh.dim-1, 0);

            tol = 1e-6;
            expected_normal = zeros(mesh.dim, 1);

            % Check that normal passed into g is correct
            for d = 1:mesh.dim
                expected_normal(:) = 0;

                % Left side
                markers = meshing.mark_entities(mesh, mesh.dim-1, [], @(x, b) b && x(d) < tol, 1);
                expected_normal(d) = -1;
                g = wrap_g_with_check_(@(x, n, c) 1, expected_normal, tol);
                assembling.assemble_neumann_source(dofmap, markers, 1, g, 0);

                % Right side
                markers = meshing.mark_entities(mesh, mesh.dim-1, [], @(x, b) b && x(d) > 1-tol, 1);
                expected_normal(d) = +1;
                g = wrap_g_with_check_(@(x, n, c) 1, expected_normal, tol);
                assembling.assemble_neumann_source(dofmap, markers, 1, g, 0);
            end
        end

    end

end


function verify_exterior_facet_integral_(t, dofmap, markers, marker_value, b)

    % We only check the value for P1 functions
    if ~strcmp(dofmap.element.family, 'Lagrange') || dofmap.element.order ~= 1
        return
    end

    % Fetch data
    facet_dofs = dofmap.element.facet_dofs;
    cell_dofs = dofmap.cell_dofs;
    c2f = dofmap.mesh.get_connectivity(dofmap.mesh.dim, dofmap.mesh.dim-1);
    boundary_facets = dofmap.mesh.get_boundary_facets();
    vertex_coords = dofmap.mesh.vertex_coords;

    % Allocate sparse result (guess <= 16 nonzeros)
    b_ = sparse([], [], [], dofmap.dim, 1, 16);

    % Compute unscaled basis function integral
    integral = basis_func_int_(dofmap.element.order, dofmap.mesh.dim-1);

    % Loop over exterior facets, scale by facet area, and store to appropriate DOF
    for f = find(markers == marker_value)

        % Find first adjacent cell and local facet index
        [lf, c] = find(c2f == f, 1);

        % Discard interior facets
        if ~boundary_facets(lf, c)
            continue
        end

        % Build dofs correspoding to the facet
        fdofs = cell_dofs(facet_dofs(lf, :), c);

        % Add contribution to integral
        area = facet_area_(vertex_coords(:, fdofs));
        b_(fdofs) = b_(fdofs) + integral*area;  %#ok<SPRIX>
    end

    t.verifyEqual(norm(b_-b), 0, 'AbsTol', 1e-14, 'RelTol', 1e-14);
end


function res = basis_func_int_(order, dim)
    if order == 0
        % Integral of constant func over unit volume
        res = 1;
    elseif order == 1
        % Integral of hat func over unit volume
        res = 1/(dim+1);
    else
        error('not impl');
    end
end


function area = facet_area_(verts)
    num_verts = size(verts, 2);

    if num_verts <= 1
        area = 1;
    elseif num_verts == 2
        area = norm(verts(:, 2) - verts(:, 1));
    elseif num_verts == 3
        % Heron's formula
        a = norm(verts(:, 2) - verts(:, 1));
        b = norm(verts(:, 3) - verts(:, 1));
        c = norm(verts(:, 3) - verts(:, 2));
        s = 0.5*(a+b+c);
        area = sqrt(s*(s-a)*(s-b)*(s-c));
    else
        error('not impl');
    end
end


function g_new = wrap_g_with_check_(g_orig, expected_normal, tol)
    g_new = @g_new_;

    function val = g_new_(x, n, c)
        if norm(expected_normal-n) > tol
            error('expected normal %s, got %s', ...
                  sprintf('%f ', expected_normal), ...
                  sprintf('%f ', n));
        end
        val = g_orig(x, n, c);
    end
end
