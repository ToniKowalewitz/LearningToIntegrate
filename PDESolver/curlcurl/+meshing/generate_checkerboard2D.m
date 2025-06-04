function [mesh, varargout] = generate_checkerboard2D(varargin)
    % Creates a half-disc containing given rectangles and points.
    %
    % Current constraints:
    % - Half-sphere is centered at [0, 0]
    % - Blocks aren't allowed to intersect.
    % - Blocks aren't allowed to form cavities.
    % - Points are expected to be placed at y = 0.
    %
    % Coordinate system: Cartesian right-hand-side
    %                    i.e. z-axis oriented upwards
    %
    % Supported (tested) Gmsh version: 4.5.6, MSH file format version 2
    %
    % SYNTAX
    %   [mesh, varargout] = generate_checkerboard(varargin)
    %
    % OPTIONAL PARAMETER
    %   block      ... Matrix [n x 4] of n blocks definitions with columns:
    %                    1-2 origin point   x, y
    %                    3-4 block extents dx, dy
    %   point      ... Matrix [n x 2] of n point coordinates.
    %   domain_r   ... Scalar, denoting the sphere radius.
    %   size_at_pt ... Scalar, denoting the cellsize around given point.
    %   refinement ... Scalar, uniform (Gmsh) mesh refinements.
    %   keep_files ... Boolean, denoting if .geo and .msh files shouldn't be
    %                     deleted.
    %   marker     ... Vector: -1 <= x_i <= 2, of geometric entities to be
    %                  additionally exported from mesh:
    %                   -1 -> vector of physical entity marker
    %                    0 -> vector of point    entity marker
    %                    1 -> vector of line     entity marker
    %                    2 -> vector of face     entity marker
    %
    % OUTPUT PARAMETER
    %   mesh      ... Object from Mesh class.
    %   varargout ... Entity marker vectors, ordered by given input marker
    %                 sorting.

    % Define possible input keys and its properties checks.
    input_keys = {'block', 'point', 'domain_r', 'marker', ...
                  'size_at_pt', 'refinement', 'keep_files'};
    assertBlock = @(x) assert(ismatrix(x) && size(x, 2) == 4);
    assertScalar = @(x) isnumeric(x) && isscalar(x) && isreal(x);
    assertPoint = @(x) ismatrix(x) && size(x, 2) == 2;
    assertMarker = @(x) isvector(x) && (all(x >= -1) && all(x <= 2));

    % Create inputParser object and set possible inputs with defaults.
    parser_obj = inputParser();
    parser_obj.addParameter(input_keys{1}, [], assertBlock);
    parser_obj.addParameter(input_keys{2}, [], assertPoint);
    parser_obj.addParameter(input_keys{3}, 1, assertScalar);
    parser_obj.addParameter(input_keys{4}, [], assertMarker);
    % FIXME: default 1 might be quite tiny.
    parser_obj.addParameter(input_keys{5}, 1, assertScalar);
    parser_obj.addParameter(input_keys{6}, 0, assertScalar);
    parser_obj.addParameter(input_keys{7}, false, @islogical);

    % Exctract all properties from inputParser.
    parse(parser_obj, varargin{:});
    args = parser_obj.Results;

    % Prepare paths
    path = fileparts(mfilename('fullpath'));
    templatefile = strcat(path, filesep, 'geocode', filesep, 'checkerboard2D.in');
    geofile = strcat(path, filesep, 'checkerboard2D.geo');
    mshfile = strcat(path, filesep, 'checkerboard2D.msh');

    % Sanity checks.
    check_input_consistency(args);
    assert(args.domain_r/args.size_at_pt <= 1e6);

    % Build geo code for domain.
    domain_r = sprintf('%.17g', args.domain_r);

    % Build geo code for blocks.
    block_geo = get_block_geo(args);

    % Build geo code for points.
    point_geo = get_point_geo(args);

    % Build geo code for mesh refinement.
    size_at_pt = sprintf('%.17g', args.size_at_pt);
    meshing_geo_code = strjoin(repmat({'RefineMesh;'}, 1, args.refinement), ...
                               newline);

    % Build geo code from template
    geo_template = fileread(templatefile);
    geo_code = sprintf(geo_template, ...
        domain_r, ...
        block_geo, ...
        point_geo, ...
        size_at_pt, ...
        meshing_geo_code);

    % Write out geo code
    f = fopen(geofile, 'w');
    fprintf(f, geo_code);
    fclose(f);

    % Run gmsh to generate msh
    cmd = sprintf('gmsh %s -save -format msh2', geofile);
    util.run_sys_cmd(cmd);

    % Read mesh
    varargout = cell(size(args.marker));
    [mesh, varargout{:}] = io.read_msh(mshfile, args.marker);

    % Clean up.
    if ~args.keep_files
        delete(geofile, mshfile);
    end
end

%% Helper

function check_input_consistency(args)

    % Check points.
    if ~isempty(args.point)
        if any(args.point(:, 2))
            error('Points need to be placed at y = 0.');
        end
        if any(vecnorm(args.point.') > args.domain_r)
            error('Some points are lying outside the domain.');
        end
    end

    % Check blocks.
    for bb = 1:size(args.block, 1)
        pt = get_corner_pt(args.block(bb, :));
        % Below surface.
        assert(all(pt(:, 2) <= 0));
        % Within half-sphere.
        assert(all(vecnorm(pt, 2, 2) < args.domain_r));
        % No intersection with other blocks.
        % FIXME: may use Gilbert–Johnson–Keerthi algorithm.
    end
end

function pt = get_corner_pt(block)
    % Get list of corner points from block definition in Gmsh format.

    pt = [block(1)         , block(2)         ;
          block(1)         , block(2)+block(4);
          block(1)+block(3), block(2)         ;
          block(1)+block(3), block(2)+block(4)];
end

function fs = make_fspec(num, spec)
    % Create string of 'num' comma separated format specifier 'spec'.

    fs = strjoin(repmat({spec}, num, 1), ', ');
end

%% Geocode helper.

function point_geo = get_point_geo(args)
    % Get points and point list .geo syntax for .geo file template.

    pt = args.point;
    n_pt = size(pt, 1);
    point_geo = cell(n_pt+1, 1);
    for i = 1:n_pt
        point_geo{i} = sprintf('Point(point_id+%d) = {%.17g, %.17g, 0};', ...
                                i-1, pt(i, :));
    end
    point_geo{n_pt+1} = sprintf(['point_list() += {', ...
                                make_fspec(n_pt, 'point_id+%d'), '};'], ...
                                0:n_pt-1);
    point_geo = strjoin(point_geo, newline);
end

function block_geo = get_block_geo(args)
    % Get box(es) .geo syntax for .geo file template.

    bl = args.block;
    n_bl = size(bl, 1);
    block_geo = cell(n_bl+1, 1);
    for i = 1:n_bl
        block_geo{i} = sprintf('Rectangle(block_id+%d) = {%.17g, %.17g, 0, %.17g, %.17g};', ...
                                i-1, bl(i, :));
    end
    block_geo{n_bl+1} = sprintf(['block_list() += {', make_fspec(n_bl, 'block_id+%d'), '};'], ...
                                 0:n_bl-1);
    block_geo = strjoin(block_geo, newline);
end
