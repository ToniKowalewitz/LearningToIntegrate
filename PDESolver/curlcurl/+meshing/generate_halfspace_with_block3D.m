function [mesh, varargout] = generate_halfspace_with_block3D(varargin)
    % Creates a 3D half-space containing a (dipped) block.
    %
    % The defined block is allowed to pierce through the surface of the
    % halfspace which will cut of the upper part of the block.
    %
    % Coordinate system: Cartesian right-hand-side
    %                    i.e. z-axis oriented upwards
    %
    % Supported (tested) Gmsh version: 4.5.6, MSH file format version 2
    %
    % SYNTAX
    %   [mesh, varargout] = generate_halfspace_with_block3D(varargin)
    %
    % OPTIONAL PARAMETER
    %   block_o    ... Vector [3 x 1] of block origin, i.e. the corner
    %                  point (xmin, ymin, zmin).
    %   block_dz   ... Scalar block height.
    %   block_dx   ... Scalar block width.
    %   block_dy   ... Scalar block length.
    %   block_a    ... Scalar block dipping angle w.r.t. x.
    %   point      ... Matrix [n x 3], coordinates of points in mesh.
    %   refinement ... Scalar, uniform (Gmsh) mesh refinements.
    %   domain_r   ... Scalar, denoting half of the domain width and
    %                      domain depth.
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
    input_keys = {'block_o', 'block_dx', 'block_dy', 'block_dz', ...
                  'block_a', 'point', 'domain_r', 'marker', ...
                  'refinement', 'keep_files', 'size_at_pt'};
    assertParam = @(x) assert(isscalar(x) && x > 0);
    assertCenter = @(x) assert(isvector(x) && length(x) == 3);
    assertPoint = @(x) assert(ismatrix(x) && size(x, 2) == 3);
    assertMarker = @(x) isvector(x) && (all(x >= -1) && all(x <= 3));
    assertScalar = @(x) isscalar(x) && x > 0;

    % Create inputParser object and set possible inputs with defaults.
    parser_obj = inputParser();
    parser_obj.addParameter(input_keys{1}, [-1, -1, -10], assertCenter);
    parser_obj.addParameter(input_keys{2}, 2, assertParam);
    parser_obj.addParameter(input_keys{3}, 2, assertParam);
    parser_obj.addParameter(input_keys{4}, 10, assertParam);
    parser_obj.addParameter(input_keys{5}, 0, @isscalar);
    parser_obj.addParameter(input_keys{6}, [], assertPoint);
    parser_obj.addParameter(input_keys{7}, 5000, assertParam);
    parser_obj.addParameter(input_keys{8}, [], assertMarker);
    parser_obj.addParameter(input_keys{9}, 0, @isscalar);
    parser_obj.addParameter(input_keys{10}, false, @islogical);
    parser_obj.addParameter(input_keys{11}, 10, assertScalar);

    % Exctract all properties from inputParser.
    parse(parser_obj, varargin{:});
    args = parser_obj.Results;

    % Sanity check.
    assert(isempty(args.point) || ...
           all(vecnorm(args.point, 2, 2) < args.domain_r));
    assert(args.domain_r/args.size_at_pt <= 1e6);

    % Check if block may intersect other boundary than surface or lies
    % above the surface.
    args.block_extent = [args.block_dx, args.block_dy, args.block_dz];
    args.block_centre = args.block_o + args.block_extent./2;
    block_w = round(norm(args.block_extent));
    bl_range = sqrt(args.block_dx^2 + args.block_dz^2);
    if ~isempty(args.point)
        assert((args.block_centre(3) - bl_range) > ...
                (min(args.point(:, 3)) - args.domain_r));
        assert((args.block_centre(3) - bl_range) < ...
                max(args.point(:, 3)));
        assert((args.block_centre(2) + args.block_dy) < ...
                (max(args.point(:, 2)) + args.domain_r));
        assert((args.block_centre(2) - args.block_dy) > ...
                (min(args.point(:, 2)) - args.domain_r));
        assert((args.block_centre(1) + bl_range) < ...
                (max(args.point(:, 1)) + args.domain_r));
        assert((args.block_centre(1) - bl_range) > ...
                (min(args.point(:, 1)) - args.domain_r));
    end

    % Prepare paths
    path = fileparts(mfilename('fullpath'));
    templatefile = strcat(path, filesep, 'geocode', filesep, 'halfspace_with_block3D.in');
    geofile = strcat(path, filesep, 'halfspace_with_block3D.geo');
    mshfile = strcat(path, filesep, 'halfspace_with_block3D.msh');

    % Build geo code for points.
    if ~isempty(args.point) && any(args.point(:, 3) ~= 0)
        warning(['No topography supported, setting z = 0 for all ', ...
                 'given points.']);
        args.point(:, 3) = 0;
    end
    points_geo_code = {};
    for i = 1:size(args.point, 1)
        points_geo_code{i} = sprintf('Point(%d) = {%.17g, %.17g, %.17g};', ...
                                     i, args.point(i, :));  %#ok<AGROW>
    end
    points_geo_code = strjoin(points_geo_code, newline);

    % Build geo code for domain.
    domain_geo_code = get_domain(args.point, args.domain_r);

    % Build geo code for block.
    block_geo_code = get_block(args);

    % Build geo code for meshing and refinement by splitting.
    meshing_geo_code = strjoin(repmat({'RefineMesh;'}, ...
                               1, args.refinement), newline);

    % Build geo code from template
    geo_template = fileread(templatefile);
    geo_code = sprintf(geo_template, ...
        num2str(args.size_at_pt, '%.17g'), ...
        num2str(args.domain_r, '%.17g'), ...
        num2str(block_w, '%.17g'), ...
        points_geo_code, ...
        domain_geo_code, ...
        block_geo_code, ...
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

function domain_geo = get_domain(points, domain_r)
    % Defines box by centre and extents.

    % Define domain boundaries.
    xmin = -domain_r;
    xmax =  domain_r;
    ymin = -domain_r;
    ymax =  domain_r;
    if isempty(points)
        zmax = 0;
    else
        zmax = max(points(:, 3));
    end
    zmin = -domain_r;

    domain_origin = [xmin, ymin, zmin];
    domain_extent = [xmax-xmin, ymax-ymin, zmax-zmin];

    domain_geo = sprintf(['Box(1) = {%d', repmat(', %.17g', 1, 5), '};'], ...
                                             domain_origin, domain_extent);
end

function block_geo = get_block(args)
    % Defines box by centre and extents.

    block_geo = sprintf(['Box(2) = {%d', repmat(', %.17g', 1, 5), '};'], ...
                                          args.block_o, args.block_extent);
    if args.block_a ~= 0
        % Define block rotation.
        rotate_geo = sprintf(['Rotate{{0, 1, 0}, {%.17g, %.17g, %.17g}, ', ...
                              '(%.17g * Pi)/180} {Volume{2};}'], ...
                                          args.block_centre, args.block_a);
        % Summarize.
        block_geo = {block_geo};
        block_geo = [block_geo, {rotate_geo}];
        block_geo = strjoin(block_geo, newline);
    end
end
