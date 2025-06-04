function [mesh, varargout] = read_tetgen(filename, output)
    % Read mesh info from tetgen output file set.
    %
    % SYNTAX
    %   [mesh, [varargout]] = read_tetgen(filename, output)
    %
    % INPUT PARAMETER
    %   filename ... Char, filename of tetgen output (with extenstion as
    %                mentioned below).
    %   output   ... Vector, elements denoting output entity dimension:
    %                   -2 raw data
    %                    0 vertex   markers
    %                    1 line     markers
    %                    2 face     markers
    %                    3 volume   markers
    %
    % OUTPUT PARAMETERS
    %   mesh      ... object from Mesh class
    %   varargout ... entity markers
    %                    -2 -> Struct, containing .mat file info
    %                     0 -> Vector of markers
    %                     1 -> Vector of markers
    %                     2 -> Vector of markers
    %                     3 -> Vector of markers
    %
    % REMARKS
    %
    %   Returned markers are zero for entities not marked.
    %   Cell markers is dense vector, other markers are sparse.
    %
    % FORMAT
    %   for filename there are assumed to exist files with extenstion:
    %   [required]
    %       .node ... contains ALL vertices with coordinates and markers
    %            header: num vtx, num coords(3), num attributes, bnd marker
    %             table: id, x, y, z, ..., maker
    %       .ele  ... contains ALL vertex2cell connectivities with markers
    %            header: num tet, num vtx per tet (4 or 10), region marker
    %             table: id, vtx1, vtx2, vtx3, vtx4, ..., maker
    %   [optional]
    %       .face ... contains ALL vertex2face connectivities with markers
    %            header: num face, bnd marker
    %             table: id, vtx1, vtx2, vtx3, maker
    %       .edge ... contains ALL vertex2edge connectivities with markers
    %            header: num edge, bnd marker
    %             table: id, vtx1, vtx2, maker
    %   markers
    %       By default TetGen boundary markers (which will be handed down
    %       each entity belonging to a PLC) have values 0 or 1.
    %       It is assumed that only markers differing from default are
    %       required.
    %
    %       For identifying
    %           nodes ... attribute marker (only the first)
    %           cells ... region marker
    %           faces ... boundary marker
    %           edges ... boundary marker
    %       are used.
    %
    %       Note, that user defined edge boundary markers are preserved
    %       afer meshing, whereas node boundary markers are overwritten by
    %       markers of superior entities.

    if nargin < 2
        output = [];
    end

    % Set required and optional file extentions.
    req_files = {'node', 'ele'};
    opt_files = {'face', 'edge'};

    %% Check input.

    assert(all(cellfun(@(x) exist([filename, '.', x], 'file'), ...
           req_files) == 2));

    %% Load from files.

    % Vertex and vtx2cell info.
    data = struct();
    for ii = 1:length(req_files)
        data.(req_files{ii}) = read_tetgen_file(req_files{ii});

        % Ensure that counting starts from 1.
        data.(req_files{ii}) = sanity_check(data, req_files{ii});
    end

    % Edge and face info.
    for ii = 1:length(opt_files)
        if exist([filename, '.', opt_files{ii}], 'file')
            data.(opt_files{ii}) = read_tetgen_file(opt_files{ii});

            % Ensure that counting starts from 1.
            data.(opt_files{ii}) = sanity_check(data, opt_files{ii});
        end
    end

    %% Build mesh object.

    mesh = meshing.Mesh(3, data.(req_files{1}).table(:, 2:4).', ...
                        uint32(data.(req_files{2}).table(:, 2:5).'));

    %% Extract mesh topology.

    % Prepare and check variable output
    n_outs = length(output);
    nargoutchk(n_outs+1, n_outs+1);
    varargout = cell(1, n_outs);
    for mm = 1:n_outs
        switch output(mm)
            case -2
            % Extract raw data.
            varargout{mm} = data;

            case {0}
            % Extract vertex ids.
            varargout{mm} = extract_entity_marker(mesh, ...
                data.(req_files{1}), output(mm));

            case {3}
            % Extract cell ids.
            varargout{mm} = extract_entity_marker(mesh, ...
                data.(req_files{2}), output(mm));

            case {2}
            % Extract face ids.
            if isfield(data, opt_files{1})
                varargout{mm} = map_entity_marker(mesh, ...
                    data.(opt_files{1}), output(mm));
            end

            case {1}
            % Extract line ids.
            if isfield(data, opt_files{2})
                varargout{mm} = map_entity_marker(mesh, ...
                    data.(opt_files{2}), output(mm));
            end
        end
    end

    function info = read_tetgen_file(extension)
        fullname = [filename, '.', extension];
        info = struct();

        % Read header line.
        fid = fopen(fullname, 'r');
        info.header = str2num(fgetl(fid)); %#ok<ST2NM>
        assert(isvector(info.header) && all(~isnan(info.header)));
        if info.header(2) == 10
            error('Tetrahedra definition with 10 vertices not supported.');
        end

        % Compatibility for MATLAB <= 2018b
        if verLessThan('matlab', '9.6')
            readmatrix_ = @(varargin) table2array(readtable(varargin{:}));
        else
            readmatrix_ = @readmatrix;
        end

        % Read file content.
        opts = detectImportOptions(fullname, ...
                                   'FileType', 'text', ...
                                   'CommentStyle', '#');
        info.table = readmatrix_(fullname, opts);
        fclose(fid);
    end
end

function em = extract_entity_marker(mesh, data, dim)
    % OUTPUT
    %   em ... Vector [mesh.num_entities(dim) x 1] of entity ids.
    %   wherein
    %       id ... entity marker set by user
    %              (TetGen default markers are ignored)

    % Fetch.
    head = data.header;
    num_entities = mesh.num_entities(dim);

    % Set vector.
    if dim ~= 3
        % Node.
        if head(3) == 0 % use attribute markers
            [i, s] = deal([]);
        else
            % Find non default marker.
            marker = data.table(:, 5); % use first attribute markers
            i = find(marker ~= 0 & marker ~= 1);
            s = marker(i);
        end
        em = sparse(i, 1, s, num_entities, 1);
    else
        % Cells.
        if head(end) == 0 % use boundary markers
            em = zeros(num_entities, 1);
        else
            em = data.table(:, end);
            % Expecting non default marker are set.
            if any(em == 0 | em == 1)
                warning('TetGen default cell marker obtained.');
            end
        end
    end
end

function em = map_entity_marker(mesh, data, dim)
    % OUTPUT
    %   em ... Sparse vector [mesh.num_entities(dim) x 1] of entity ids.
    %   wherein
    %       id ... entity marker set by user
    %              (TetGen default markers are ignored)

    % Fetch.
    head = data.header;
    if head(end) ~= 0
        marker = data.table(:, end);
        data = data.table(:, 2:end-1);
    end

    % Get connectivities.
    w = warning('off', 'Mesh:ExtraConnStored');  % Mute warning
    mesh.compute_connectivity(dim, 0);
    warning(w);  % Restore warning
    entities_to_vertices = mesh.get_connectivity(dim, 0);
    num_entities = mesh.num_entities(dim);

    if head(end) == 0
        % Skip.
        [marker_map, marker_vals] = deal([]);
    else
        % Find non default marker.
        i = find(marker ~= 0 & marker ~= 1);
        marker_vals = marker(i);
        nondef_data = data(i, :);

        % Check and fetch data.
        assert(size(nondef_data, 2) == dim+1);
        cur_entity_markers_vertices = sort(nondef_data, 2); % ensure ascending numbering

        % Compute entity indices of marked entities
        [~, marker_map] = ismember(cur_entity_markers_vertices, ...
                                   entities_to_vertices.', 'rows');
        assert(~any(~marker_map), ...
            'Some entity-to-vertex relations in file are not part of the mesh.');
    end
    % Build sparse vector representation
    em = sparse(marker_map, 1, marker_vals, num_entities, 1);
end

function info = sanity_check(info, field)

    info = info.(field);
    assert(info.header(1) == size(info.table, 1));

    % Check entity id.
    if info.table(1, 1) == 0 && info.table(end, 1) == info.header(1)-1
        info.table(:, 1) = info.table(:, 1) + 1;
    elseif info.table(1, 1) == 1 && info.table(end, 1) == info.header(1)
        % Skip.
    else
        error('Unknown data type.');
    end

    % Check entity definitions.
    cols = 2:2+info.header(2)-1;
    if ~strcmp(field, 'node') && any(any(info.table(:, cols) == 0))
        info.table(:, cols) = info.table(:, cols) + 1;
    end
end
