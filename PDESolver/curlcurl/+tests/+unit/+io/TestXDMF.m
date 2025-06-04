classdef TestXDMF < matlab.unittest.TestCase

    properties (TestParameter)
        mesh_factory = {
            {@meshing.generate_half_disc_mesh, 1         };
            {@meshing.generate_half_disc_mesh, 8         };
            {@meshing.generate_half_disc_mesh, 32        };
            {@meshing.generate_unit_cube_mesh, [ 1, 1   ]};
            {@meshing.generate_unit_cube_mesh, [ 8, 8   ]};
            {@meshing.generate_unit_cube_mesh, [32,32   ]};
            {@meshing.generate_unit_cube_mesh, [ 1, 1, 1]};
            {@meshing.generate_unit_cube_mesh, [ 4, 4, 4]};
            {@meshing.generate_unit_cube_mesh, [ 8, 8, 8]};
        };

        element_factory = {
            {nedelec, 1};
            {nedelec, 2};
            {dg, 0};
            {dg, 1};
            {dg, 2};
            {cg, 1};
            {cg, 2};
        };
    end

    methods (Test)

        function test_write_function(t, mesh_factory, element_factory)

            % Extract parameters
            [mesh_factory, mesh_resolution] = mesh_factory{:};
            [element_factory, element_order] = element_factory{:};

            % Build mesh
            mesh = mesh_factory(mesh_resolution);

            % Build element
            element = element_factory(mesh.dim, element_order);

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

            % Check that these fail
            t.verifyError(@() t.t_write_function_internal_({'hdf5'}, dofmap), ?MException);  % not yet implemented
            t.verifyError(@() t.t_write_function_internal_({'quux'}, dofmap), ?MException);  % bogus encoding

            % Run actual test with various encodings
            t.t_write_function_internal_({        }, dofmap);
            t.t_write_function_internal_({'binary'}, dofmap);
            t.t_write_function_internal_({'xml'   }, dofmap);

        end

    end

    methods

        function t_write_function_internal_(t, xdmf_args, dofmap)

            % Initialize FE function dofs
            vec = zeros(dofmap.dim, 1);

            % Prepare filepaths
            filename = 'test_xdmf';
            path = fileparts(mfilename('fullpath'));
            xdmffile = fullfile(path, sprintf('%s.xdmf', filename));

            % Prepare XDMF handler
            xdmf = io.XDMF(xdmffile, xdmf_args{:});

            % Verify XDMF options
            t.verifyEqual(xdmf.filename, xdmffile);
            if isempty(xdmf_args)
                t.verifyMatches(xdmf.encoding, 'binary');  % binary is default
            else
                t.verifyMatches(xdmf.encoding, xdmf_args{1});
            end

            % Test first write of FE function
            vec(:) = 1;
            time = 30.0;
            xdmf.write(dofmap, vec, time);

            % Test second write of FE function
            vec(:) = 2;
            time = 60.0;
            xdmf.write(dofmap, vec, time);

            % Test flush
            xdmf.flush();

            % Test write after flush
            vec(:) = 3;
            time = 90.0;
            xdmf.write(dofmap, vec, time);

            % FIXME: How to test correctness of XDMF output robustly?
            % WORKAROUND: Run this (case-by-case) and check correctness in Paraview
            fun = get_interpolant(dofmap.element.simplex.dim, dofmap.element.value_shape);
            vec = assembling.interpolate(dofmap, fun);
            time = 120.0;
            xdmf.write(dofmap, vec, time);

            % Check that final flush by delete grows the file
            f0 = dir(xdmffile);
            clear xdmf;  % flush!
            f1 = dir(xdmffile);
            t.verifyGreaterThan(f1.bytes, f0.bytes);

        end

    end

end


function fun = nedelec()
    fun = @fe.create_nedelec_element;
end


function fun = dg()
    fun = @fe.create_discontinuous_lagrange_element;
end


function fun = cg()
    fun = @fe.create_lagrange_element;
end


function fun = get_interpolant(dim, value_shape)
    switch dim
    case 2
        switch value_shape
        case 1
            fun = @(x, c) sin(pi*x(1))*sin(2*pi*x(2));
        case 2
            fun = @(x, c) [sin(pi*x(1))*sin(2*pi*x(2)), ...
                           cos(pi*x(1))*sin(2*pi*x(2))];
        otherwise
            error('not implemented');
        end
    case 3
        switch value_shape
        case 1
            fun = @(x, c) sin(pi*x(1))*sin(2*pi*x(2))*sin(4*pi*x(3));
        case 3
            fun = @(x, c) [sin(pi*x(1))*sin(2*pi*x(2))*sin(4*pi*x(3)), ...
                           cos(pi*x(1))*sin(2*pi*x(2))*sin(4*pi*x(3)), ...
                           sin(pi*x(1))*cos(2*pi*x(2))*sin(4*pi*x(3))];
        otherwise
            error('not implemented');
        end
    otherwise
        error('not implemented');
    end
end
