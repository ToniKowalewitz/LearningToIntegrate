function [mesh, varargout] = generate_mesh2D(varargin)
    % Stores mesh information of a mesh created by Gmsh.
    %
    % Creates a homogeneous disc whose surface curve is formed by
    % incorporating points and topography information.
    %
    % Coordinate system: Cartesian right-handed
    %                    i.e. y-axis oriented upwards
    %
    % The center of the disc is defined by
    %   x:               args.domain_c(1)
    %   if not available (max(t)-min(t))/2, t = x-coo of [point; topo]
    %   if not available 0
    %   y:               args.domain_c(2)
    %   if not available mean(t),           t = y-coo of [point; topo]
    %   if not available 0
    %
    % Supported (tested) Gmsh version: 4.5.6, MSH file format version 2
    %
    % SYNTAX
    %   [mesh, varargout] = generate_mesh2D(varargin)
    %
    % OPTIONAL PARAMETER
    %   point      ... Matrix [m x 2] or [m x 4], with:
    %                  [x, y(, is_on_surface, id)]
    %                  Given as [m x 2]: points are expected to describe
    %                                    single points on surface.
    %                  Given as [m x 4]: points can be placed at domain
    %                                    interior -> 'is_on_surface' column
    %                                    points can describe lines by
    %                                    assigning similar 'id' to multiple
    %                                    points   -> 'id'            column
    %   topo       ... Matrix [k x 2] of topography describing points.
    %   ref        ... Scalar, denoting the number of uniform mesh
    %                  refinement steps.
    %   domain_r   ... Scalar, denoting the radius of domain (disc).
    %   domain_c   ... Vector [1 x 2] of center of the domain (disc).
    %   size_at_pt ... Scalar, denoting the cell sizes at points.
    %   keep_files ... Logical, denoting if .geo and .msh files shouldn't
    %                  be deleted.
    %   marker     ... Vector: -1 <= x_i <= 2, of geometric entities to be
    %                  additionally exported from mesh:
    %                   -1 -> vector of physical entity marker
    %                    0 -> vector of point    entity marker
    %                    1 -> vector of line     entity marker
    %                    2 -> vector of face     entity marker
    %
    % OUTPUT PARAMETER
    %   mesh      ... object from Mesh class.
    %   varargout ... Entity markers with decreasing dimension.

    %% Check input.

    % Define possible input keys and its properties checks.
    input_keys = {'ref', 'point', 'topo', 'domain_r', 'domain_c', ...
                  'size_at_pt', 'marker', 'keep_files'};
    assertRef = @(x) assert(isscalar(x) && ~islogical(x) && x >= 0, ...
        'ref - Scalar, denoting number of uniform ref steps, expected.');
    assertPos = @(x) assert(ismatrix(x) && (size(x, 2) == 2 || ...
        (size(x, 2) == 4 && all((x(:, 3) == 1 | x(:, 3) == 0))) ...
                         && all(x(:, 4) > 0)), ...
        ['point - Matrix [n x 2] (or [n x 4]) with columns ', ...
        '[x, y(, insitu, id)] expected.']);
    assertTopo = @(x) assert(ismatrix(x) && size(x, 2) == 2, ...
        'topo - Vector [n x 2] of surface points, expected.');
    assertScalar = @(x) isscalar(x) && x > 0;
    assertMarker = @(x) isvector(x) && (all(x >= -1) && all(x <= 2));

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
    tmpfile = strcat(path, filesep, 'geocode', filesep, 'mesh2D.in');
    geofile = strcat(path, filesep, 'mesh2D.geo');
    mshfile = strcat(path, filesep, 'mesh2D.msh');

    %% Fetch.

    % Standardize input.
    args = standardize_input(args);

    %% Create geo code parts.

    % Build geo code for points.
    [point_geo, args] = get_point(args);

    % Build geo code for domain.
    face_geo = get_face(args);

    % Build geo code for refinement by splitting.
    cell_size_at_point_geo = num2str(args.size_at_pt, '%.17g');
    refine_geo = get_meshing(args);

    %% Build .geo from template.

    % Insert geo code parts into template.
    geo_template = fileread(tmpfile);
    geo_code = sprintf(geo_template, ...
                       point_geo, ...
                       face_geo, ...
                       cell_size_at_point_geo, ...
                       refine_geo);

    % Write .geo file.
    f = fopen(geofile, 'w');
    fprintf(f, geo_code);
    fclose(f);

    %% Run gmsh to generate msh.

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
    %
    % Note: A fifth colum is added to the [m x 4] representation,
    %       containing a point identifier flag:
    %           1 = point
    %           2 = topo

    % Extend columns.
    if ~isempty(args.point)
        n = size(args.point, 1);
        if size(args.point, 2) == 2
            args.point = [args.point, 1+zeros(n, 1), (1:n).'];
        end
        args.point = [args.point, 1+zeros(n, 1)];
    end
    if ~isempty(args.topo)
        n = size(args.topo, 1);
        args.topo = [args.topo, 1+zeros(n, 1), (1:n).', 2+zeros(n, 1)];
    end
end

function fs = make_fspec(num, spec)
    % Create string of 'num' comma separated format specifier 'spec'.

    fs = strjoin(repmat({spec}, num, 1), ', ');
end

function out = get_mid(in)
    % Calculate centerpoint of 'in' from its max/min extent.

    [min_in, max_in] = bounds(in);
    out = mean([min_in, max_in]);
end

function pt = remove_duplicates(pt)
    % Check for duplicate points, preferably keep point.

    % Get pairwise point distances and indices of coinciding points.
    match = check4duplicates(pt);
    if isempty(match)
        return;
    end
    % FIXME: also handle this case!
    assert(all(groupcounts(match(:)) == 1), ...
           'Some points occur more than twice.');

    % Remove non point, if possible.
    to_remove = false(size(pt, 1), 1);
    for m = 1:size(match, 1)
        is_point = pt(match(m, :), end) ~= 2;
        if any(is_point) && ~all(is_point)
            % Remove topo point.
            to_remove(match(m, ~is_point), :) = true;
        else
            % Remove first point.
            to_remove(match(m, 1), :) = true;
        end
    end
    pt(to_remove, :) = [];
end

function match = check4duplicates(pt)
    % Check for duplicate lines in pt(:, 1:2).

    match = [];
    for i = 2:size(pt, 1)
        for j = 1:i-1
            if vecnorm(pt(i, 1:2) - pt(j, 1:2)) <= eps
                match = [match; i, j]; %#ok<AGROW>
            end
        end
    end
end

function pt = fetch_point_info_from_args(args)
    % Extracts coodrinates from args and sort it in an array.

    % Merge point lists from all input args.
    tmp = {args.point, args.topo};
    pt = [];
    for ii = 1:length(tmp)
        if isempty(tmp{ii})
            continue;
        end

        % Remove line objects.
        id_occur = groupcounts(tmp{ii}(:, end-1));
        is_ln = id_occur > 1;
        if any(is_ln)
            warning('MESHING:generate_mesh2D:UnsupportedObjects', ...
                    ['"point" contains unsupported line ', ...
                     'objects - these are omitted.']);
            tmp{ii} = tmp{ii}(~is_ln, :);
        end
        pt = [pt; tmp{ii}];  %#ok<AGROW>
    end
end

function [pt, args] = ensure_gmsh_surf(pt, args)
    % Add/remove info from pt array and args to ensure valid Gmsh surface.

    % Get all surface points.
    if isempty(pt)
        is_surf = logical([]);
    else
        is_surf = logical(pt(:, 3));
    end

    % If no surface points given, define dummy at (given) center.
    if ~any(is_surf)
        if isempty(args.domain_c)
            args.domain_c = [0, 0];
        end
        tmp = [args.domain_c, 0, 0, 0];
        dummy = true;
    else
        % Get surface and insitu points.
        tmp = pt(is_surf, :);
        dummy = false;

        % Set domain center w.r.t. surface.
        % FIXME: What would be a better choice?
        if isempty(args.domain_c)
            args.domain_c = [get_mid(tmp(:, 1)), mean(tmp(:, 2))];
        end
    end
    pt_c = args.domain_c;

    % Check if points at surface exist that lie outside domain extent.
    is_valid_surf = (1:size(tmp, 1)).';
    surf_out = vecnorm((tmp(:, 1:2) - pt_c).') > ...
                       args.domain_r;
    if any(surf_out)
        is_valid_surf(surf_out) = [];
        warning('MESHING:generate_mesh2D:SurfPointOutsideDomain', ...
                'Surface point outside domain extent (omited).');
    end

    % Add auxiliary points (outside domain) for meshing surface.
    get_x = @(idx) sqrt((1.5 * args.domain_r)^2-(tmp(idx, 2)-pt_c(2))^2);
    idx_min = find(tmp(:, 1) == min(tmp(:, 1)), 1, 'first');
    idx_max = find(tmp(:, 1) == max(tmp(:, 1)), 1, 'first');
    new_x_min = pt_c(1) - get_x(idx_min);
    new_x_max = pt_c(1) + get_x(idx_min);
    tmp_min = [new_x_min, tmp(idx_min, 2), 1, max(tmp(:, 4))+1, 3];
    tmp_max = [new_x_max, tmp(idx_max, 2), 1, max(tmp(:, 4))+1, 3];

    % Check if at least auxiliary points vertical position lies within the
    % domain vertical extent to be able to cut disc by curve.
    assert(all(abs(tmp([idx_min, idx_max], 2) - pt_c(2)) < args.domain_r));

    % Check if points at domain interior lie outside domain extent.
    is_valid_insitu = find(~is_surf);
    if ~isempty(is_valid_insitu)
        insitu_out = vecnorm((pt_c - pt(is_valid_insitu, 1:2)).') > ...
                          args.domain_r;
        if any(insitu_out)
            is_valid_insitu(insitu_out) = [];
            warning('MESHING:generate_mesh2D:InsituPointOutsideDomain', ...
                    'Insitu point outside domain extent (omited).');
        end
    end

    % Sort surface points w.r.t. x-coordinate (required by Gmsh).
    tmp = [tmp(is_valid_surf, :); tmp_min; tmp_max];
    if dummy
       tmp(1, :) = [];
    end
    [~, map] = sort(tmp(:, 1));
    pt = [tmp(map, :); pt(is_valid_insitu, :)];

    % Remove duplicates (required by Gmsh).
    pt = remove_duplicates(pt);
end

function gc = write_pt(pt)
    % Defines points and point lists in .geo syntax to .geo file template.

    % Fetch.
    n = size(pt, 1);
    id = 0:(n-1);
    is_point = pt(:, 5) == 1;
    is_insitu = is_point & pt(:, 3) == 0;
    is_topo = logical(pt(:, 3));

    % Point definitions.
    tmp_gc = cell(size(pt, 1), 1);
    for i = 1:n
        tmp_gc{i} = sprintf(['Point(pt_id+%d) = {', ...
                             make_fspec(2, '%.17g'),', 0};'], ...
                             i-1, pt(i, 1:2));
    end

    % Helper.
    pt2str = @(x) sprintf(make_fspec(length(id(x)), 'pt_id+%d'), id(x));

    % Append point point list.
    if any(is_point)
        tmp_gc = [tmp_gc; ['pt_list() += {', pt2str(is_point), '};']];
    end

    % Append global topography describing point list.
    if any(is_topo)
        tmp_gc = [tmp_gc; ['topo_pt_list() += {', ...
                           pt2str(is_topo), '};']];
    end

    % Append global insitu point list.
    if any(is_insitu)
        tmp_gc = [tmp_gc; ['insitu_pt_list() += {',...
                           pt2str(is_insitu), '};']];
    end

    % Add physical point entities.
    if any(is_point)
        tmp_gc = [tmp_gc; ['Physical Point("point", 1) = {', ...
                           pt2str(is_point), '};']];
    end

    % Write geo code.
    gc = strjoin(tmp_gc, newline);
end

%% Geo code helper.

function [point_geo, args] = get_point(args)

    % Fetch info from given arguments.
    pt_all = fetch_point_info_from_args(args);

    % Ensure min set of sorted points for surface generation in Gmsh.
    [pt_all, args] = ensure_gmsh_surf(pt_all, args);

    % Write point info into .geo file template.
    point_geo = write_pt(pt_all);
end

function face_geo = get_face(args)

    % Write domain info into .geo file template.
    face_geo = [sprintf('domain_r = %.17g;', args.domain_r), newline, ...
                sprintf('Disk(fc_id) = {%.17g, %.17g, 0, domain_r};', ...
                       args.domain_c)];
end

function refine_geo = get_meshing(args)

    % Write refinement info into .geo file template.
    refine_geo = strjoin(repmat({'RefineMesh;'}, 1, args.ref), newline);
end
