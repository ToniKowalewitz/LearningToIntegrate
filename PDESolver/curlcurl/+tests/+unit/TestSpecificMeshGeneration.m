classdef TestSpecificMeshGeneration < matlab.unittest.TestCase

    properties (TestParameter)
        % NB: Second column is number of expected cell marker values
        factory = {
            {@checkerboard2D};
            {@checkerboard3D};
            {@halfspace_with_block3D};
            {@halfspace_with_block2D};
        };
    end

    methods (Test)

        function test_models_default(t, factory)

            [gen_mesh, args, ~] = factory{:}();

            % Generate mesh
            [mesh, pn, pm, fm, ~] = gen_mesh(args(end-1:end));

            % Test that topology of mesh is good
            verify_mesh_topology(t, mesh);

            % Check that markers span expected range
            t.verifyEqual(full(unique(fm)).', 0:2);
            t.verifyEqual(full(unique(pm)).', 0);
            t.assert_pn(pn, args{end});
        end

        function test_models(t, factory)

            [gen_mesh, args, num_cm] = factory{:}();

            % Generate mesh
            [mesh, pn, pm, fm, cm] = gen_mesh(args);

            % Test that topology of mesh is good
            verify_mesh_topology(t, mesh);

            % Check that coordinates span expected box
            switch mesh.dim
            case 2
                m = [-1e5; -1e5];
                M = [+1e5;    0];
            case 3
                m = [-1e5; -1e5; -1e5];
                M = [+1e5; +1e5;    0];
            end
            t.assert_in_box(mesh.vertex_coords, m, M);

            % Check that markers span expected range
            t.verifyEqual(full(unique(cm)).', 1:num_cm);
            t.verifyEqual(full(unique(fm)).', 0:2);
            t.verifyEqual(full(unique(pm)).', 0:1);
            t.assert_pn(pn, args{end});
        end

    end

    methods

        function assert_pn(t, pn, marker)
            % FIXME: Test correctness of pn?
            t.verifyTrue(iscell(pn));
            t.verifyTrue(length(pn) == length(marker)-1);
        end

        function assert_in_box(t, coordinates, min_, max_)
            % Check that coodinates are bounded in box
            [m, M] = bounds(coordinates, 2);
            t.verifyEqual(m, min_, 'RelTol', 3e-2, 'AbsTol', 5e-9);
            t.verifyEqual(M, max_, 'RelTol', 3e-2, 'AbsTol', 5e-9);
        end

    end
end


function [gen_mesh, args, num_cm] = checkerboard2D()
    width = 50;
    X1 = -100-width/2:  width: 75-width/2;
    X2 = -125-width/2:2*width:125-width/2;
    X3 =  -75-width/2:2*width: 75-width/2;
    blocks =          make_blocks(X1, 1);
    blocks = [blocks; make_blocks(X2, 2)];
    blocks = [blocks; make_blocks(X3, 2.5)];
    num_cm = 10;

    x_pt = linspace(-150, 150, size(X1, 2));
    y_pt = zeros(numel(x_pt), 1);
    points = [x_pt(:), y_pt(:)];

    args = {
        'point',        points,           ...
        'block',        blocks,           ...
        'domain_r',     1e5,              ...
        'keep_files',   true,             ...
        'size_at_pt',   5,                ...
        'marker',       [-1, 0, 1, 2],    ...
    };
    gen_mesh = @(inputs) meshing.generate_checkerboard2D(inputs{:});

    function bl = make_blocks(X, level)
        tmp = zeros(numel(X), 1);
        dxy = tmp + width;
        Y = tmp - level*width;
        bl = [X(:), Y, dxy, dxy];
    end
end


function [gen_mesh, args, num_cm] = checkerboard3D()
    width = 50;
    [x, y] = deal(-50-width/2:width:75-width/2);
    [X1, Y1] = ndgrid(x, y);
    X1 = X1(:);
    Y1 = Y1(:);
    X2 = X1((numel(X1)+1)/2);
    Y2 = Y1((numel(Y1)+1)/2);
    X1((numel(X1)+1)/2) = [];
    Y1((numel(Y1)+1)/2) = [];
    blocks =          make_blocks(X1, Y1, 1);
    blocks = [blocks; make_blocks(X2, Y2, 2)];
    num_cm = 10;

    [x_pt, y_pt] = deal(linspace(-50, 50, 5));
    [X_pt, Y_pt] = ndgrid(x_pt, y_pt);
    Z_pt = zeros(numel(X_pt), 1);
    points = [X_pt(:), Y_pt(:), Z_pt(:)];

    args = {
        'point',        points,           ...
        'block',        blocks,           ...
        'domain_r',     1e5,              ...
        'keep_files',   true,             ...
        'size_at_pt',   5,                ...
        'marker',       [-1, 0, 2, 3],    ...
    };
    gen_mesh = @(inputs) meshing.generate_checkerboard3D(inputs{:});

    function bl = make_blocks(X, Y, level)
        tmp = zeros(numel(X), 1);
        dxyz = tmp + width;
        Z = tmp - level*width;
        bl = [X, Y, Z, dxyz, dxyz, dxyz];
    end
end


function [gen_mesh, args, num_cm] = halfspace_with_block2D()
    num_pt = 21;
    pt = [linspace(-50, 50, num_pt); zeros(1, num_pt)];
    num_cm = 2;
    args = {
        'point',        pt.',             ...
        'block_c',      [15, -10],        ...
        'block_h',      20,               ...
        'block_w',      10,               ...
        'block_a',      45,               ...
        'domain_r',     1e5,              ...
        'keep_files',   true,             ...
        'size_at_pt',   5,                ...
        'marker',       [-1, 0, 1, 2],    ...
    };
    gen_mesh = @(inputs) meshing.generate_halfspace_with_block2D(inputs{:});
end


function [gen_mesh, args, num_cm] = halfspace_with_block3D()
    num_pt = 21;
    pt = [linspace(-50, 50, num_pt); zeros(1, num_pt); zeros(1, num_pt)];
    num_cm = 2;
    args = {
        'point',        pt.',             ...
        'block_o',      [10, -20, -10],   ...
        'block_dx',     10,               ...
        'block_dy',     40,               ...
        'block_dz',     20,               ...
        'block_a',      45,               ...
        'domain_r',     1e5,              ...
        'keep_files',   true,             ...
        'size_at_pt',   5,                ...
        'marker',       [-1, 0, 2, 3],    ...
    };
    gen_mesh = @(inputs) meshing.generate_halfspace_with_block3D(inputs{:});
end
