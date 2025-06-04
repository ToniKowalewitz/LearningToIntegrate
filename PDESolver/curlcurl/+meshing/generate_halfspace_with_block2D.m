function [mesh, varargout] = generate_halfspace_with_block2D(varargin)
    % Creates a 2D half-space containing a (dipped) block.
    %
    % For given points, variing in two coordinats, also a halfspace with
    % topography can be constructed.
    % The defined block is allowed to pierce through the surface of the
    % halfspace which will cut of the upper part of the block.
    %
    % Coordinate system: Cartesian right-hand-side
    %                    i.e. z-axis oriented upwards
    %
    % Supported (tested) Gmsh version: 4.5.6, MSH file format version 2
    %
    % SYNTAX
    %   [mesh, varargout] = generate_halfspace_with_block2D(varargin)
    %
    % OPTIONAL PARAMETER
    %   block_c    ... Vector [2 x 1] of dike centre.
    %   block_h    ... Scalar block height.
    %   block_w    ... Scalar block width.
    %   block_a    ... Scalar block dipping angle.
    %   point      ... Matrix [n x 2], coordinates of points in mesh.
    %   refinement ... Scalar, uniform (Gmsh) mesh refinements.
    %   domain_r   ... Scalar, denoting half of the domain width and
    %                  domain depth.
    %   size_at_pt ... Scalar, denoting the cell sizes at points.
    %   keep_files ... Boolean, denoting if .geo and .msh files shouldn't be
    %                  deleted.
    %   marker ... Vector: -1 <= x_i <= 3, of geometric entities to be
    %              additionally exported from mesh:
    %                 -1 -> vector of physical entity marker
    %                  0 -> vector of point    entity marker
    %                  1 -> vector of line     entity marker
    %                  2 -> vector of face     entity marker
    %                  3 -> vector of volume   entity marker
    %
    % OUTPUT PARAMETER
    %   mesh      ... Object from Mesh class.
    %   varargout ... Entity marker vectors, ordered by given input marker
    %                 sorting.

    % Define possible input keys and its properties checks.
    input_keys = {'block_c', 'block_h', 'block_w', 'block_a'...
                  'point', 'domain_r', 'marker', 'refinement', ...
                  'keep_files', 'size_at_pt'};
    assertParam = @(x) assert(isscalar(x) && x > 0);
    assertCenter = @(x) assert(isvector(x) && length(x) == 2);
    assertPoint = @(x) assert(ismatrix(x) && size(x, 2) == 2);
    assertMarker = @(x) isvector(x) && (all(x >= -1) && all(x <= 3));
    assertScalar = @(x) isscalar(x) && x > 0;

    % Create inputParser object and set possible inputs with defaults.
    parser_obj = inputParser();
    parser_obj.addParameter(input_keys{1}, [0, 0], assertCenter);
    parser_obj.addParameter(input_keys{2}, 10, assertParam);
    parser_obj.addParameter(input_keys{3}, 2, assertParam);
    parser_obj.addParameter(input_keys{4}, 0, @isscalar);
    parser_obj.addParameter(input_keys{5}, [], assertPoint);
    parser_obj.addParameter(input_keys{6}, 5000, assertParam);
    parser_obj.addParameter(input_keys{7}, [], assertMarker);
    parser_obj.addParameter(input_keys{8}, 0, @isscalar);
    parser_obj.addParameter(input_keys{9}, false, @islogical);
    parser_obj.addParameter(input_keys{10}, 10, assertScalar);

    % Exctract all properties from inputParser.
    parse(parser_obj, varargin{:});
    args = parser_obj.Results;

    % Sanity check.
    assert(isempty(args.point) || ...
           all(vecnorm(args.point, 2, 2) < args.domain_r));
    assert(args.domain_r/args.size_at_pt <= 1e6);

    % Check if block may intersect other boundary than surface or lies
    % above the surface.
    bl_range = sqrt(args.block_h^2 + args.block_w^2);
    if ~isempty(args.point)
        assert((args.block_c(2) - bl_range) > (min(args.point(:, 2)) - ...
                                               args.domain_r));
        assert((args.block_c(2) - bl_range) < max(args.point(:, 2)));
        assert((args.block_c(1) + bl_range) < (max(args.point(:, 1)) + ...
                                               args.domain_r));
        assert((args.block_c(1) - bl_range) > (min(args.point(:, 1)) - ...
                                               args.domain_r));
    end

    % Prepare paths
    path = fileparts(mfilename('fullpath'));
    templatefile = strcat(path, filesep, 'geocode', filesep, 'halfspace_with_block2D.in');
    geofile = strcat(path, filesep, 'halfspace_with_block2D.geo');
    mshfile = strcat(path, filesep, 'halfspace_with_block2D.msh');

    % Build geo code for points
    [pt, is_ele, is_surf] = add_domain_bnd(args.point, args.domain_r);
    points_geo_code = {};
    for i = 1:size(pt, 1)
        points_geo_code{i} = sprintf('Point(%d) = {%.17g, %.17g, 0};', ...
                                     i, pt(i, :));  %#ok<AGROW>
    end
    points_geo_code = strjoin(points_geo_code, newline);
    make_format_char = @(number) strjoin(repmat({'%d'}, number, 1), ', ');
    ele_list = sprintf(['{', make_format_char(size(args.point, 1)), '}'], ...
                       find(is_ele));
    surf_list = sprintf(['{', make_format_char(length(find(is_surf))), '}'], ...
                       find(is_surf));

    % Build geo code for meshing and refinement by splitting.
    meshing_geo_code = strjoin(repmat({'RefineMesh;'}, ...
                               1, args.refinement), newline);

    % Build geo code from template
    geo_template = fileread(templatefile);
    geo_code = sprintf(geo_template, ...
        sprintf('{%.17g, %.17g, 0}', args.block_c), ...
        num2str(args.block_w, '%.17g'), ...
        num2str(args.block_h, '%.17g'), ...
        num2str(args.block_a, '%.17g'), ...
        num2str(args.size_at_pt, '%.17g'), ...
        num2str(args.domain_r, '%.17g'), ...
        points_geo_code, ele_list, surf_list, meshing_geo_code);

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

function [pt, is_ele, is_surf] = add_domain_bnd(points, domain_r)
    % Add domain boundary points to point list.

    % Define domain boundaries.
    if ~isempty(points)
        pt_xmin = points(:, 1) == min(points(:, 1));
        pt_xmax = points(:, 1) == max(points(:, 1));
        y_xmin = points(pt_xmin, 2);
        y_xmax = points(pt_xmax, 2);
    else
        [y_xmin, y_xmax] = deal(0);
    end
    xmin = -domain_r;
    xmax =  domain_r;
    ymin = -domain_r;
    points_bnd = [xmin, y_xmin;
                  xmax, y_xmax;
                  xmin, ymin;
                  xmax, ymin];

    % Summarize and sort clockwise.
    pt = [points; points_bnd];
    is_ele = [true(size(points, 1), 1); false(size(points_bnd, 1), 1)];
    is_surf = [true(size(points, 1), 1); [true; true; false; false]];
    map = sort_clockwise(pt);
    pt = pt(map, :);
    is_ele = is_ele(map);
    is_surf = is_surf(map);
end

function map = sort_clockwise(pt)
    % Sort points by referring to angle w.r.t. center point.

    % FIXME: assumes that points form convex surface.
    pt_mean_x = mean(pt(:, 1));
    pt_mean_y = mean(pt(:, 2));
    pt_angle = atan2(pt(:, 2) - pt_mean_y, pt(:, 1) - pt_mean_x);
    [~, map] = sort(pt_angle, 'descend');
end
