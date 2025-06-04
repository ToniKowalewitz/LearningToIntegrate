classdef TestAssembleMassSensitivity < matlab.unittest.TestCase

    properties (TestParameter)
        factory = {
            {@meshing.generate_half_disc_mesh, 3        , @nedelec, 1};
            {@meshing.generate_half_disc_mesh, 2        , @nedelec, 2};
            {@meshing.generate_unit_cube_mesh, [2, 2   ], @nedelec, 1};
            {@meshing.generate_unit_cube_mesh, [1, 1   ], @nedelec, 2};
            {@meshing.generate_unit_cube_mesh, [2, 2, 2], @nedelec, 1};
            {@meshing.generate_unit_cube_mesh, [1, 1, 1], @nedelec, 2};
            {@meshing.generate_half_disc_mesh, 3        , @lagrange, 1};
            {@meshing.generate_half_disc_mesh, 2        , @lagrange, 2};
            {@meshing.generate_half_disc_mesh, 2        , @lagrange, 3};
            {@meshing.generate_unit_cube_mesh, [2, 2   ], @lagrange, 1};
            {@meshing.generate_unit_cube_mesh, [1, 1   ], @lagrange, 2};
            {@meshing.generate_unit_cube_mesh, [1, 1   ], @lagrange, 3};
            {@meshing.generate_unit_cube_mesh, [2, 2, 2], @lagrange, 1};
            {@meshing.generate_unit_cube_mesh, [1, 1, 1], @lagrange, 2};
            {@meshing.generate_unit_cube_mesh, [1, 1, 1], @lagrange, 3};
        };
    end

    methods (Test)

        function test_assemble_mass_sensitivity(t, factory)
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

            % Run actual tests
            verify_values(t, dofmap)
            verify_conventions(t, dofmap)
        end

    end

end


function verify_values(t, dofmap)
    % Compare assemble_mass_sensitivity() and assemble_curl_curl() cell-by-cell
    num_cells = dofmap.mesh.num_entities(dofmap.mesh.dim);
    sigma = sparse(num_cells, 1);
    I = eye(dofmap.dim);
    DA = assembling.assemble_mass_sensitivity(dofmap, I, I);
    for c = 1:num_cells
        sigma(c) = 1.01;  %#ok<SPRIX>
        A = assemble_reference_operator(dofmap, sigma);
        t.verifyEqual(norm(A(:)-DA*sigma, inf), 0, 'AbsTol', 1e-15);
        sigma(c) = 0;  %#ok<SPRIX>
    end
end


function verify_conventions(t, dofmap)
    % Check the indexing convention and the complex conjugation convention
    % of assemble_mass_sensitivity()

    % Create several coefficients matrices (complex, constant in space)
    U = zeros(dofmap.dim, 2); U(:, 1) = 4+1i; U(:, 2) = 5+2i;
    V = zeros(dofmap.dim, 3); V(:, 1) = 5+1i; V(:, 2) = 6+2i; V(:, 3) = 7+3i;

    % Check cols of U and V are lexicographically packed in rows of A1
    A1 = sum(assembling.assemble_mass_sensitivity(dofmap, U, V), 2);
    A2 = V.'*(assemble_reference_operator(dofmap, 1)*U);
    t.verifyEqual(norm(A1 - A2(:), inf), 0, 'AbsTol', 1e-12);

    % Check that transposing the arguments leads to different result
    A1 = sum(assembling.assemble_mass_sensitivity(dofmap, V, U), 2);
    t.verifyGreaterThanOrEqual(norm(A1 - A2(:), inf), 5);
    p = reshape(reshape(1:size(U, 2)*size(V, 2), size(U, 2), size(V, 2)).', [], 1);
    t.verifyEqual(norm(A1(p) - A2(:), inf), 0, 'AbsTol', 1e-12);

    % Check how complex numbers are handled
    A1 = sum(assembling.assemble_mass_sensitivity(dofmap, U, conj(V)), 2);
    A2 = V'*(assemble_reference_operator(dofmap, 1)*U);
    t.verifyEqual(norm(A1 - A2(:), inf), 0, 'AbsTol', 1e-12);
end


function A = assemble_reference_operator(dofmap, sigma)
    switch dofmap.element.mapping
    case 'affine'
        A = assembling.assemble_laplace(dofmap, 0, -sigma);
    case 'covariant'
        A = assembling.assemble_curl_curl(dofmap, -sigma, 1/0);
    otherwise
        error('not implemented!');
    end
end


function varargout = nedelec(varargin)
    [varargout{1:nargout}] = fe.create_nedelec_element(varargin{:});
end


function varargout = lagrange(varargin)
    [varargout{1:nargout}] = fe.create_lagrange_element(varargin{:});
end
