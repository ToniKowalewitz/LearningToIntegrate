classdef TestDefaultMeshGeneration < matlab.unittest.TestCase

    properties (TestParameter)
        factory = {
            {1  , 2};
            {1:2, 2};
            {1:3, 2};
            {1:4, 2};
            {1:5, 2};
            {1:6, 2};
            {1:7, 2};
            {1  , 3};
            {1:2, 3};
            {1:3, 3};
            {1:4, 3};
            {1:5, 3};
            {1:6, 3};
            {1:7, 3};
        };
    end

    methods (Test)

        function test_read(t, factory)

            import matlab.unittest.constraints.IssuesWarnings;

            % Get reference survey.
            [point, topo, point_info] = get_ref_input(factory{2});

            % Set up arguments.
            args = {{}, ...
                    {'point', point}, ...
                    {'topo', topo}, ...
                    {'size_at_pt', 5}, ...
                    {'domain_r', 500}, ...
                    {'ref', 1}, ...
                    {'domain_c', [500, repmat(10, 1, factory{2}-1)]}};
            test_args = [args{factory{1}}];
            marker = [-1, 0:factory{2}];

            % Warning ids.
            warn_args = {'MESHING:generate_mesh2D:UnsupportedObjects'; ...
                         'MESHING:generate_mesh2D:SurfPointOutsideDomain'; ...
                         'MESHING:generate_mesh2D:InsituPointOutsideDomain'; ...
                         'MESHING:generate_mesh3D:PointOutsideDomain'; ...
                         'MESHING:generate_mesh3D:TopoPointOutsideDomain'; ...
                         'Mesh:UnconnectedVertices'};

            % Create mesh.
            if factory{2} == 2
                if length(factory{1}) == 1
                    [mesh, phy, pm, fm, cm] = meshing.generate_mesh2D(...
                                                      test_args{:}, ...
                                                      'marker', marker);
                else
                    % Initialize multiple warnings from one function call.
                    if factory{1} < 7
                        warn_id = warn_args(1);
                    else
                        warn_id = warn_args(1:3);
                    end
                    verify_constraint = IssuesWarnings(warn_id, 'WhenNargoutIs', ...
                                               length(marker)+1, ...
                                               'RespectingOrder', true, ...
                                               'RespectingSet', true);

                    % Check verbosity.
                    t.verifyThat(@() meshing.generate_mesh2D(...
                                             test_args{:}, ...
                                             'marker', marker), ...
                                 verify_constraint);
                    [mesh, phy, pm, fm, cm] = verify_constraint.FunctionOutputs{:};
                end
            elseif factory{2} == 3
                if length(factory{1}) < 7
                    [mesh, phy, pm, lm, fm, cm] = ...
                    meshing.generate_mesh3D(test_args{:}, ...
                                            'marker', marker);
                else
                    warn_id = warn_args(4:end);
                    verify_constraint = IssuesWarnings(warn_id, ...
                                               'WhenNargoutIs', length(marker)+1, ...
                                               'RespectingOrder', true, ...
                                               'RespectingSet', true);
                    % Check verbosity.
                    t.verifyThat(@() meshing.generate_mesh3D(...
                                                test_args{:}, ...
                                                'marker', marker), ...
                                 verify_constraint);
                    [mesh, phy, pm, lm, fm, cm] = verify_constraint.FunctionOutputs{:};
                end
            else
                error(' ');
            end

            % Check mesh.
            t.verifyNotEmpty(mesh.vertex_coords);
            t.verifyNotEmpty(mesh.cells);
            t.verifyTrue(mesh.dim == factory{2});

            % Test that topology of mesh is good
            verify_mesh_topology(t, mesh);

            % Check marker.
            % cells.
            t.verifyTrue(unique(cm) == 1);
            % faces.
            t.verifyTrue(all(unique(fm) == [0; 1; 2]));
            if isempty(test_args)
                if factory{2} == 3
                    % edges (only valid in 3D case).
                    t.verifyTrue(all(unique(lm) == 0));
                end
                % points.
                t.verifyTrue(all(unique(pm) == 0));
            else
                if factory{2} == 3 && max(factory{1}) < 7
                    % edges (omitted if domain is shifted).
                    [ed_msh, ed_ref] = get_edge_id();
                    t.verifyTrue(all(unique(lm) == sort([0; ed_ref])));
                end
                if factory{2} == 2 && max(factory{1}) < 7
                    % points (omitted if domain is shifted).
                    [pt_msh, pt_ref] = get_point_id();
                    t.verifyTrue(all(unique(pm) == sort([0; pt_ref])));
                end
            end

            % Check marker info.
            t.verifyTrue(length(phy) == length(marker)-1);
            % cells.
            t.verifyTrue(unique(cm) == get_id_from_phy('domain', mesh.dim));
            % faces.
            t.verifyTrue(all([1; 2] == sort([...
                               get_id_from_phy('surface', mesh.dim-1); ...
                               get_id_from_phy('subsurface', mesh.dim-1)...
                                            ])));
            if isempty(test_args)
                % edges & points.
                t.verifyTrue(all(cellfun(@isempty, phy{1})));
            else
                if factory{2} == 3 && max(factory{1}) < 7
                    % edges.
                    t.verifyTrue(all(ed_ref == ed_msh) && ...
                                 (~isempty(ed_ref) || ~isempty(ed_msh)));
                end
                if factory{2} == 2 && max(factory{1}) < 7
                    % points
                    t.verifyTrue(all(pt_ref == pt_msh) && ...
                                 (~isempty(pt_ref) || ~isempty(pt_msh)));
                end
            end

            function id = get_id_from_phy(name, marker_num)
                % Extract stored ids for given marker names.
                vargout_idx = find(marker == marker_num)-1;
                id = phy{vargout_idx}{1}(ismember(phy{vargout_idx}{2}, ...
                                                  {name}));
            end

            function [id_msh, id_ref] = get_edge_id()
                % Extract reference ids and assign physical names.
                % (See meshing.generate_mesh3D.m / get_point)
                id_ref = [point_info.wire{:, 2}].';
                str = point_info.wire(:, 1);
                n = length(str);
                id_msh = zeros(n, 1);
                % Extract stored ids for given marker names.
                for i = 1:n
                    id_msh(i) = get_id_from_phy([str{i}, '_', ...
                                                 num2str(id_ref(i))], 1);
                end
            end

            function [id_msh, id_ref] = get_point_id()
                % Extract reference ids and assign physical names.
                % (See meshing.generate_mesh3D.m / get_point)
                id_ref = point_info.point;
                id_msh = get_id_from_phy('point', 0);
            end
        end
    end
end

function [point, topo, point_info] = get_ref_input(type)
    % Provide auxiliary information about a geophysical measurement.
    %
    % All objects (i.e. points, wires, loops) are solely described
    % by points.
    %
    % Data type
    %   point       ... Matrix [n x 4] in 3D or [n x 5] in 2D of doubles
    %   n           ... number of points
    %   with row indices
    %   1 - 2 or 3  ... x, y(, z) coordinate
    %   3 or 4      ... 0 or 1 for point located insitu or at surface
    %   4 or 5      ... tx, rx id
    %              point: each point gets an own id
    %              wire : multiple points points with same id form a wire
    %              loop : multiple points points with same id and the first
    %                     and last point coincide             form a loop
    %
    % Topography (topo) data type
    %   topo       ... Matrix [t x 2] or [t x 3] of doubles
    %   t          ... number of points
    %   with row indices
    %   1 - 2 or 3 ... x, y(, z) coordinate

    assert(any(type == [2, 3]));

    % Definitions.
    xyz = @(x, y) [x, y, sphere_cut(x, y)];

    % Topography.
    if type == 2
        X = linspace(-100, 100, 6);
        Z = sphere_cut(X, 0*X);
        topo = [X(:), Z(:)];
    else
        [X, Y] = ndgrid(linspace(-100, 100, 6));
        Z = sphere_cut(X, Y);
        topo = [X(:), Y(:), Z(:)];
    end

    % Objects.
    x = linspace(-30, 60, 5);
    if type == 2
        y = zeros(1, 5);
        z = sphere_cut(x, y);
        point = [xyz(-60, 0),      1,          1;       % point
                 -40, 0, -10,      0,          2;       % point
                 x(:), y(:), z(:), ones(5, 1), (3:7).'; % points
                 -10, 0, -15,      0,          12;      % wire (omitted)
                 -5, 0, -15,       0,          12];
        point = point(:, [1, 3:end]);
    else
        y = linspace(-70, 70, 5);
        z = sphere_cut(x, y);
        point = [xyz(-60, 25),     1,          1;       % point
                 -40, 25, -10,     0,          2;       % point
                 0,   0,  -5,      0,          3;
                 0,   0,  -10,     0,          4;
                 x(:), y(:), z(:), ones(5, 1), (5:9).'; % points
                 xyz(0, 45),       1,          10;      % wire
                 xyz(0, 50),       1,          10;
                 xyz(0, 55),       1,          10;
                 xyz(-50, 0),      1,          11;      % loop
                 xyz(-45, -5),     1,          11;
                 xyz(-50, -10),    1,          11;
                 xyz(-55, -5),     1,          11;
                 xyz(-50, 0),      1,          11;
                 -10, -10, -15,    0,          12;      % wire
                 -5,  -5,  -15,    0,          12;
                 0,   0,   -15,    0,          12];
    end
    point_info = struct();
    point_info.point = 1;
    point_info.wire = {'wire', 10;
                       'loop', 11;
                       'wire', 12};
end

function z = sphere_cut(x, y)
    % Z-coordinate of a flat surface (z=const.) comprising a sphere cut.

    % Fetch.
    [nx, ny] = size(x);
    assert(all([nx, ny] == size(y)));
    x = x(:);
    y = y(:);
    n = nx * ny;

    % Definitions.
    h = 10; % height above surface
    a = 50; % radius of topography
    r = (a^2 + h^2)/(2*h); % radius of sphere
    z0 = -r +h;            % center depht of sphere

    % Get elevation.
    z = zeros(n, 1);
    for ii = 1:n
        xy_r = norm([x(ii); y(ii)]);
        if xy_r < a
           z(ii) = sqrt(r^2 - x(ii)^2 - y(ii)^2) + z0;
        else
           z(ii) = 0;
        end
    end
    z = reshape(z, nx, ny);
end
