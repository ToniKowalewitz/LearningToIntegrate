function [mesh, varargout] = generate_mesh3D(varargin)
    % Stores mesh information of a mesh created by Gmsh.
    %
    % Creates a homogeneous half-sphere whose surface is formed by
    % incorporating points and topography information.
    %
    % Coordinate system: Cartesian right-handed
    %                    i.e. z-axis oriented upwards
    %
    % The center of the sphere is defined by
    %   x,y:             args.domain_c(1:2)
    %   if not available (max(t)-min(t))/2, t = x,y-coo of [point; topo]
    %   if not available [0, 0]
    %   z:               args.domain_c(3)
    %   if not available mean(t),           t = z-coo of [point; topo]
    %   if not available 0
    %
    % Supported (tested) Gmsh version: 4.5.6, MSH file format version 2
    %
    % SYNTAX
    %   [mesh, varargout] = generate_mesh3D(varargin)
    %
    % OPTIONAL PARAMETER
    %   point      ... Matrix [m x 3] or [m x 5], with:
    %                  [x, y, z(, is_on_surface, id)]
    %                  Given as [m x 3]: points are expected to describe
    %                                    single points on surface.
    %                  Given as [m x 5]: points can be placed at domain
    %                                    interior -> 'is_on_surface' column
    %                                    points can describe lines by
    %                                    assigning similar 'id' to multiple
    %                                    points   -> 'id'            column
    %   topo       ... Matrix [k x 3] of topography describing points.
    %   ref        ... Scalar, denoting the number of uniform mesh
    %                  refinement steps.
    %   domain_r   ... Scalar, denoting the radius of domain (sphere).
    %   domain_c   ... Vector [1 x 3] of center of domain (sphere).
    %   size_at_pt ... Scalar, denoting the cell sizes at points.
    %   keep_files ... Logical, denoting if .geo and .msh files shouldn't
    %                  be deleted.
    %   marker     ... Vector: -1 <= x_i <= 3, of geometric entities to be
    %                  additionally exported from mesh:
    %                   -1 -> vector of physical entity marker
    %                    0 -> vector of point    entity marker
    %                    1 -> vector of line     entity marker
    %                    2 -> vector of face     entity marker
    %                    3 -> vector of volume   entity marker
    %
    % OUTPUT PARAMETER
    %   mesh      ... Object from Mesh class.
    %   varargout ... Entity marker vectors, ordered by given input marker
    %                 sorting.
    %
    % REMARKS
    %   Representing spline surfaces from arbitrary point sets in Gmsh
    %   is not working in all cases!
    %   Please make sure that:
    %       - topography points are uniformly distributed across the survey
    %         area
    %       - or survey points with non-zero z-component are distributed
    %         uniformly

    %% Check input.

    % Define possible input keys and its properties checks.
    input_keys = {'ref', 'point', 'topo', 'domain_r', 'domain_c', ...
                  'size_at_pt', 'marker', 'keep_files'};
    assertRef = @(x) assert(isscalar(x) && ~islogical(x) && x >= 0, ...
        'ref - Scalar, denoting number of uniform ref steps, expected.');
    assertPos = @(x) assert(ismatrix(x) && (size(x, 2) == 3 || ...
        (size(x, 2) == 5 && all((x(:, 4) == 1 | x(:, 4) == 0))) ...
                         && all(x(:, 5) > 0)), ...
        ['point - Matrix [n x 3] (or [n x 5]) with columns ', ...
        '[x, y, z(, insitu, id)] expected.']);
    assertTopo = @(x) assert(ismatrix(x) && size(x, 2) == 3, ...
        'topo - Vector [n x 3] of surface points, expected.');
    assertScalar = @(x) isscalar(x) && x > 0;
    assertMarker = @(x) isvector(x) && (all(x >= -1) && all(x <= 3));

    % Create inputParser object and set possible inputs with defaults.
    parser_obj = inputParser();
    parser_obj.addParameter(input_keys{1}, 0, assertRef);
    parser_obj.addParameter(input_keys{2}, [], assertPos);
    parser_obj.addParameter(input_keys{3}, [], assertTopo);
    parser_obj.addParameter(input_keys{4}, 5e3, assertScalar);
    parser_obj.addParameter(input_keys{5}, [], assertPos);
    parser_obj.addParameter(input_keys{6}, 10, assertScalar);
    parser_obj.addParameter(input_keys{7}, [], assertMarker);
    parser_obj.addParameter(input_keys{8}, false, @islogical);

    % Exctract all properties from inputParser.
    parse(parser_obj, varargin{:});
    args = parser_obj.Results;

    % Sanity check.
    assert(args.domain_r/args.size_at_pt <= 5e4);

    % Prepare paths
    path = fileparts(mfilename('fullpath'));
    tmpfile = strcat(path, filesep, 'geocode', filesep, 'mesh3D.in');
    tmpfilepath_line = strcat(path, filesep, 'geocode', filesep, 'mesh3D_line2');
    geofile = strcat(path, filesep, 'mesh3D.geo');
    mshfile = strcat(path, filesep, 'mesh3D.msh');

    %% Fetch.

    % Standardize input.
    args = standardize_input(args);

    % Set domain center.
    args = set_domain_center(args);

    %% Create geo code parts.

    % Build geo code for points.
    point_geo = get_point(args);

    % Build geo code for surface.
    surface_geo = get_surface(args);

    % Build geo code for line objects on surface.
    line_geo = get_line(args, tmpfilepath_line);

    % Build geo code for refinement by splitting.
    cell_size_at_point_geo = num2str(args.size_at_pt, '%.17g');
    refine_geo = get_meshing(args);

    %% Build .geo from template.

    % Insert geo code parts into template.
    geo_template = fileread(tmpfile);
    geo_code = sprintf(geo_template, ...
                       point_geo, ...
                       surface_geo, ...
                       line_geo, ...
                       cell_size_at_point_geo, ...
                       refine_geo);

    % Write .geo file.
    f = fopen(geofile, 'w');
    fprintf(f, geo_code);
    fclose(f);

    %% Run Gmsh to generate msh.

    cmd = sprintf('gmsh %s -save -format msh2', geofile);
    util.run_sys_cmd(cmd);

    %% Read mesh.

    varargout = cell(size(args.marker));
    [mesh, varargout{:}] = io.read_msh(mshfile, args.marker);

    %% Clean up.

    if ~args.keep_files
        delete(geofile, mshfile);
    end
end

%% General helper.

function args = standardize_input(args)
    % Ensure format of given point info to be [m x 5].

    % Extend columns.
    if ~isempty(args.point)
        if size(args.point, 2) == 3
            n = size(args.point, 1);
            args.point = [args.point, ones(n, 1), (1:n).'];
        end
    end
end

function fs = make_fspec(num, spec)
    % Create string of 'num' comma separated format specifier 'spec'.

    fs = strjoin(repmat({spec}, num, 1), ', ');
end

function gc = make_new_pt_tag(gc)
    % Reset point tags (ensure them not to be reused).

    gc = [gc, 'pt_tag = pt_tag+1;', newline];
end

function r = offset2pt(inp, pt, dim)
    % Calculates the distance between all points in 'inp' and point 'pt'.

    assert(size(inp, 2) == size(pt, 2));
    n = size(inp, 1);
    r = zeros(n, 1);
    for ii = 1:n
        r(ii) = norm(inp(ii, 1:dim) - pt(1:dim));
    end
end

function args = set_domain_center(args)
    % Get domain center from point info.

    if isempty(args.domain_c)
        % Fetch.
        if ~isempty(args.point)
            point = args.point(:, 1:3);
            point_in_surf = logical(args.point(:, 4));
        else
            point = [];
            point_in_surf = logical([]);
        end

        % Is any point part of surface?
        point_at_surf = point(point_in_surf, :);
        if isempty(point_at_surf)
            % Set to default [0, 0, 0].
            args.domain_c = zeros(1, 3);
        else
            % Get centre of surface describing points.
            args.domain_c = [get_mid(point_at_surf(:, 1)), ...
                             get_mid(point_at_surf(:, 2)), ...
                             mean(point_at_surf(:, 3), 1)];
        end
    end

    % Remove all given points and topo points which lie outside the domain.
    if ~isempty(args.point)
        is_inside = vecnorm(args.point(:, 1:3)-args.domain_c, 2, 2) < ...
                    args.domain_r;
        is_inside = remove_entire_line_object(is_inside, args.point(:, 5));
        if ~all(is_inside)
            warning('MESHING:generate_mesh3D:PointOutsideDomain', ...
                    'Point outside domain extent (omited).');
            args.point = args.point(is_inside, :);
        end
    end
    if ~isempty(args.topo)
        is_inside = vecnorm(args.topo(:, 1:3)-args.domain_c, 2, 2) < ...
                    args.domain_r;
        if ~all(is_inside)
            warning('MESHING:generate_mesh3D:TopoPointOutsideDomain', ...
                    'Topography point outside domain extent (omited).');
            args.topo = args.topo(is_inside, :);
        end
    end
end

function is_inside = remove_entire_line_object(is_inside, id)
    % Remove entire wire if parts lie outside the domain extent.
    %
    % Note: Otherwise it may be, that remaining points form a wire of
    %       different shape.
    %
    % FIXME: Differ between points/wire on surface and insitu?

    id_unique = unique(id);
    for i = id_unique.'
        cur_id = id == i;
        is_inside(cur_id) = all(is_inside(cur_id));
    end
end

function out = get_mid(in)
    % Calculate centerpoint of 'in' from its max/min extent.

    [min_in, max_in] = bounds(in);
    out = mean([min_in, max_in]);
end

function pt = create_dummy_pt(c, r)
   % Create set of points on circle with centre c and radius r.

   alpha = linspace(0, 2, 9).' * pi;
   alpha(end) = [];
   x = c(1) + r * cos(alpha);
   y = c(2) + r * sin(alpha);
   z = c(3) + zeros(length(alpha), 1);
   pt = [x, y, z];
end

%% Geo code helper.

function point_geo = get_point(args)

    % Initialize.
    point = args.point;

    % Write TX points.
    point_geo = write_point_geo_with_id('', point, 'point');
end

function surface_geo = get_surface(args)

    % Fetch.
    topo = args.topo;
    domain_r = args.domain_r;
    domain_c = args.domain_c;

    % Extract points on surface.
    if ~isempty(args.point)
        is_surf = logical(args.point(:, 4));
        point = args.point(is_surf, 1:3);
    else
        point = [];
    end

    if ~isempty(topo)
        % Get center of all topography describing points.
        all_surf_pt = [point; topo];
        topo_c = [get_mid(all_surf_pt(:, 1)), ...
                  get_mid(all_surf_pt(:, 2)), ...
                  mean(all_surf_pt(:, 3))];

        % Sanity check.
        assert(abs(mean(topo_c(:, 3)) - args.domain_c(3)) < args.domain_r, ...
               ['Sphere cutting by surface fails. ', ...
                'Shortest distance between sphere center and surface ', ...
                'is larger than sphere radius.']);

        % Get sorted (counterclockwise w.r.t. horz. center) boundary of
        % topo points.
        % Note: starting point occures twice!
        bnd_idx = boundary(topo(:, 1), topo(:, 2));
        inner_idx = setdiff(1:size(topo, 1), bnd_idx);

        % Get offsets to center for boundary points.
        topo_horz_bnd_r = offset2pt(topo(bnd_idx, :), topo_c, 2);

        % Define dilation factor such that (small) surface part, formed by
        % topo points, is sufficiently enlarged to be entirely cutted by
        % the domain sphere.
        dilate_fac = floor((args.domain_r + norm(args.domain_c - topo_c)) / ...
                             min(topo_horz_bnd_r));

    else
        % Create dummy points on circle with radius > domain_r.
        dilate_fac = 0;
        topo = create_dummy_pt(domain_c, 1.5 * domain_r);
        bnd_idx = [1:size(topo, 1), 1];
        inner_idx = [];
    end

    % Write topo points and (append) lists.
    surface_geo = write_surface_geo(topo);
    if ~isempty(inner_idx)
        % Append global topography describing point list.
        surface_geo = [surface_geo, '// Append lists.', newline, ...
                sprintf(['topo_pt_list() += {', ...
                         make_fspec(length(inner_idx), 'pt_id+%d'), ...
                         '};'], inner_idx-1), newline];
    end

    % Write additional topo point info.
    surface_geo = [surface_geo, ...
                   '// Add surface parameter.', newline];
    surface_geo = [surface_geo, ...
                   sprintf(['surf_bnd_pt_list() = {', ...
                            make_fspec(length(bnd_idx), 'pt_id+%d')], ...
                            bnd_idx-1), '};', newline];
    surface_geo = [surface_geo, ...
                   sprintf('dilate_fac = %.17g', ...
                           dilate_fac), ...
                           ';', newline];
    surface_geo = [surface_geo, ...
                   sprintf('domain_r = %.17g', domain_r), ';', newline];
    surface_geo = [surface_geo, ...
                   sprintf(['domain_c() = {', ...
                            make_fspec(length(domain_c), '%.17g')], domain_c), ...
                            '};', newline];
    if exist('topo_c', 'var')
        surface_geo = [surface_geo, ...
                       sprintf(['topo_c() = {', ...
                                make_fspec(length(topo_c), '%.17g')], ...
                                           topo_c), ...
                                '};', newline];
    else
        surface_geo = [surface_geo, 'topo_c() = domain_c();', newline];
    end
end

function line_geo = get_line(args, tmp)

    % Initialize.
    point = args.point;
    line_geo = '';
    if isempty(point)
        return;
    end

    % Get domain radius and mean survey depth to define default
    % z-coordinate for all points, describing any line entity.
    % Note: pt_z need to be below the expected surface of the domain!
    pt_z = mean(point(:, 3)) - 1.5 * args.domain_r;

    % Create geo code.
    if ~isempty(point)
        line_geo = write_line_geo_with_id(line_geo, point, pt_z, tmp);
    end
end

function refine_geo = get_meshing(args)

    refine_geo = strjoin(repmat({'RefineMesh;'}, 1, args.ref), newline);
end

function gc = write_point_geo(gc, pt, type, is_topo)
    % Defines points and point lists in .geo syntax to .geo file template.

    % Fetch.
    list_name = ['pt_', type, '_list()'];

    % Point definitions.
    tmp_gc = cell(size(pt, 1)+1, 1);
    tmp_gc{1} = 'pt_id = newp;';
    for i = 1:size(pt, 1)
        tmp_gc{i+1} = sprintf(['Point(pt_id+%d) = {', ...
                               make_fspec(size(pt, 2), '%.17g'),'};'], ...
                              i-1, pt(i, :));
    end

    % Point list definition.
    tmp_gc = [tmp_gc; sprintf([list_name, ' = {', ...
                               make_fspec(size(pt, 1), ...
                               'pt_id+%d'),'};'], ...
                              0:size(pt, 1)-1)];

    % Physical entity definition.
    tmp_gc = [tmp_gc; ['Physical Point("', type,'", pt_tag) = ', ...
                       list_name, ';']];

    % Append global point list.
    tmp_gc = [tmp_gc; '// Append lists.'];
    tmp_gc = [tmp_gc; ['pt_list() += ', list_name, ';']];

    i = 0:size(pt, 1)-1;
    if any(is_topo)
        % Append global topography describing point list.
        tmp_gc = [tmp_gc; sprintf(['topo_pt_list() += {', ...
                                    make_fspec(length(i(is_topo)), ...
                                    'pt_id+%d'), '};'], i(is_topo))];
    end
    if any(~is_topo)
        % Append global insitu point list.
        tmp_gc = [tmp_gc; sprintf(['insitu_pt_list() += {', ...
                                    make_fspec(length(i(~is_topo)), ...
                                    'pt_id+%d'), '};'], i(~is_topo))];
    end

    % Write geo code.
    gc = [gc, strjoin(tmp_gc, newline), newline];
end

function gc = write_point_geo_with_id(gc, inp, type)
    % Wrapper for write_point_geo, providing unique point ids.

    if isempty(inp)
        return;
    end
    gc = [gc, newline, ['// Add ', type, ' point definitions.'], ...
          newline];

    % Identify points within inp.
    ids = inp(:, end);
    [id_occur, id_unique] = groupcounts(ids(:));
    is_pt = id_occur == 1;
    id_unique = id_unique(is_pt);
    pt = inp(any(ids == id_unique.'), :);

    % Create geo code.
    gc = make_new_pt_tag(gc);
    gc = write_point_geo(gc, pt(:, 1:3), type, logical(pt(:, 4)));
end

function gc = write_surface_geo(pt)
    % Defines points and point lists in .geo syntax to .geo file template.

    % Point definitions.
    tmp_gc = cell(size(pt, 1)+1, 1);
    tmp_gc{1} = 'pt_id = newp;';
    for i = 1:size(pt, 1)
        tmp_gc{i+1} = sprintf(['Point(pt_id+%d) = {', ...
                               make_fspec(size(pt, 2), '%.17g'),'};'], ...
                               i-1, pt(i, :));
    end

    % Write geo code.
    gc = ['// Add topo point definitions.', newline];
    gc = [gc, strjoin(tmp_gc, newline), newline];
end

function gc = write_line_geo(pt)
    % Defines points and point lists in .geo syntax to .geo file template.

    % Write point definitions.
    gc = cell(size(pt, 1), 1);
    for i = 1:size(pt, 1)
        gc{i} = sprintf(['Point(pt_id+%d) = {', ...
                          make_fspec(size(pt, 2), '%.17g'),'};'], ...
                          i-1, pt(i, :));
    end
    gc = [strjoin(gc, newline), newline];
end

function gc = write_line_geo_with_id(gc, inp, pt_z, tmp)
    % Wrapper for write_line_geo, providing unique line ids and names.
    %
    % Furthermore, the required .geo template file is chosen.

    % Find multiple occurences of point ids, denoting line entities.
    [id_occur, id_unique] = groupcounts(inp(:, 5));
    is_line = id_occur > 1;
    if ~any(is_line)
        return;
    end

    % Find respective points.
    id = id_unique(is_line);
    is_id = inp(:, 5) == id.';

    % Loop over single line entities.
    for ii = 1:size(is_id, 2)

        % Sanity check.
        is_at_surf = logical(inp(is_id(:, ii), 4));
        assert(all(is_at_surf) == any(is_at_surf), ...
               ['Line object parts are on surface and at domain ', ...
                'interior.']);

        % Get all points of current line.
        pt = inp(is_id(:, ii), 1:3);

        % Set line type.
        if (all(pt(1, :) == pt(end, :)))
            name = 'loop';
        else
            name = 'wire';
        end

        % Get line type.
        if pt(1, :) == pt(end, :)
            % Loop.
            assert(size(pt, 1) > 2);
            line_type = '1';
            pt(end, :) = [];
        else
            % Wire.
            line_type = '0';
        end

        % Set template file path.
        if all(is_at_surf)
            tmp_path = [tmp, 'surf.in'];

            % Set all depths to default value.
            pt(:, 3) = pt_z;
        else
            tmp_path = [tmp, 'vol.in'];
        end

        % Initialize list of line segments and type.
        gc = [gc, ['// Line ', num2str(id(ii)), '.'], newline]; %#ok<*AGROW>
        gc = [gc, 'wire_line_list() = {};', newline];
        gc = [gc, ['form_loop = ', line_type, ';'], newline];

        % Write points.
        n = size(pt, 1);
        pt_gc = write_line_geo(pt);
        pt_gc = [pt_gc, ['tmp_aux_pt() = {', ...
                         sprintf(make_fspec(n, 'pt_id+%d'), ...
                         (1:n)-1), '};']];

        % Insert geo code parts into template.
        geo_template = fileread(tmp_path);
        gc = [gc, sprintf(geo_template, pt_gc)];

        % Summarize all segments to a single physical curve.
        cur_id = num2str(id(ii));
        cur_name = [name, '_', cur_id];
        gc = [gc, ['Physical Curve("', cur_name, '", ', cur_id, ...
                   ') = wire_line_list();'], newline];
    end
end
