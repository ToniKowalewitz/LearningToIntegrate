classdef TestMeshRead < matlab.unittest.TestCase

    properties (TestParameter)
        factory = {
            {2, '2D',          @test_read_msh   };
            {2, '2D_cm_fm',    @test_read_msh   };
            {3, '3D_cm_fm_vm', @test_read_msh   };
            {3, '3D_cm_fm_vm', @test_read_mat   };
            {3, '3D_cm_fm_vm', @test_read_tetgen};
        };
    end

    methods (Test)

        function test_read(t, factory)
            % Extract parameters
            [dim, filename, reader] = factory{:};

            % Prepare filepath
            path = fileparts(mfilename('fullpath'));

            % Read
            [mesh, pn, em] = verifyWarningFree(t, @() reader(dim, path, filename));

            % Apply some consistency checks
            assert(length(pn) == length(em));
            for ee = 1:length(pn)
                assert(length(em{ee}) == mesh.num_entities(ee-1));
            end
        end


        function test_read_mesh_unconnected_vertices_warning(t)
            % Check read_msh() with unconnected vertices emits warning

            path = fileparts(mfilename('fullpath'));
            mshfile = fullfile(path, '2D_unconnected_vertices.msh');
            mesh = verifyWarning(t, @() io.read_msh(mshfile), 'Mesh:UnconnectedVertices');  %#ok<NASGU>
        end
    end

end

function [mesh, pn, em] = test_read_msh(dim, path, filename)
    % Gmsh mesh creation from .geo and read from .msh.

    % Create (run Gmsh)
    geofile = prepare_filepath(path, filename, 'geo');
    cmd = sprintf('gmsh -%d -f msh2 %s', dim, geofile);
    util.run_sys_cmd(cmd);

    % Read
    mshfile = prepare_filepath(path, filename, 'msh');
    [mesh, pn, em{1:dim+1}] = io.read_msh(sprintf('%s', mshfile), [-1, 0:dim]);
end


function [mesh, pn, em] = test_read_mat(dim, path, filename)
    % Mesh read from .mat.

    % Read
    matfile = prepare_filepath(path, filename, 'mat');
    [mesh, pn, em{1:dim+1}] = io.read_mat(sprintf('%s', matfile), [-1, 0:dim]);
end


function [mesh, pn, em] = test_read_tetgen(dim, path, filename)
    % Mesh read from .mat.

    % Prepare TetGen command
    tetfile = prepare_filepath(path, filename, 'poly');
    cmd = sprintf('tetgen -pqiAef %s', dim, tetfile);

    % Prevent MATLAB<=2017b foisting incompatible libstdc++ upon us
    cmd = ['env --unset=LD_LIBRARY_PATH ', cmd];

    % Run command
    util.run_sys_cmd(cmd);

    % Read
    tetfile = prepare_filepath(path, filename, '1');
    [mesh, em{1:dim+1}] = io.read_tetgen(sprintf('%s', tetfile), 0:dim);
    % Dummy.
    pn = cell(size(em));
end


function matfile = prepare_filepath(path, filename, file_extension)
     matfile = fullfile(path, sprintf('%s.%s', filename, file_extension));
end
