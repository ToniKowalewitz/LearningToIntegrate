function [mesh, varargout] = read_msh(filename, output)
    % Read out mesh information from Gmsh .msh file (ascii, format Ver. 2).
    %
    % SYNTAX
    %   [mesh, [varargout]] = read_msh(filename, [output])
    %
    % INPUT PARAMETER
    %   filename ... Char, filename of .msh file
    %   output   ... Vector, elements denoting output entity dimension:
    %                   -1 physical names for markers
    %                    0 vertex   markers
    %                    1 line     markers
    %                    2 face     markers
    %                    3 volume   markers
    %
    % OUTPUT PARAMETERS
    %   mesh      ... object from Mesh class
    %   varargout ... entity markers
    %                    -1 -> Cell {[marker].' , ['physical name'].'}
    %                     0 -> Vector of markers
    %                     1 -> Vector of markers
    %                     2 -> Vector of markers
    %                     3 -> Vector of markers
    %
    % REMARKS
    %
    %   The "markers" are associated with Gmsh "physical groups".
    %   Returned markers are zero for entities not marked.
    %   Cell markers is dense vector, other markers are sparse.
    %   Marker names cells for each requested marker entity including a
    %   vector of physical markers and a cell array  of respective physical
    %   names.
    %
    %   For Gmsh 4, run
    %
    %     gmsh -<dim> -f msh2 <file.geo>
    %
    %   to get the supported format version 2.

    if nargin < 2
        output = [];
    end

    % Parse file content
    [ele_content, vtx_content, phys_content] = read_sections(filename);
    topology_data = extract_topology_data(ele_content);

    % Extract mesh topology
    [dim, cells] = extract_topology(topology_data);
    if any(dim < output)
       error('Requesting marker exceeds dimensionality of current mesh.');
    elseif any(output < -1)
       error('Some outputs are not supported.');
    end

    % Extract vertex coordinates
    vertices = extract_vertex_coordinates(vtx_content);
    if dim < 3
        % Prune coordinates 3d-padding by (almost) zeros
        vertices = prune_z_component(vertices);
        assert(size(vertices, 2) == dim);
    end

    % Build mesh object
    mesh = meshing.Mesh(dim, vertices.', cells.');

    % Prepare and check variable output
    n_outs = length(output);
    nargoutchk(n_outs+1, n_outs+1);
    varargout = cell(1, n_outs);

    for mm = 1:n_outs
        switch output(mm)
            case -1
            % Extract marker names
            varargout{mm} = extract_marker_names(phys_content, ...
                                                 output(1:end ~= mm));

            case dim
            % Extract cell markers
            varargout{mm} = extract_cell_markers(dim, topology_data);

            case num2cell(0:dim-1)
            % Extract lower dim markers
            [entity_markers_values, entity_markers_vertices] = ...
                extract_entity_markers(dim, output(mm), topology_data);

            % Build return values
            varargout{mm} = compute_entity_markers(mesh, output(mm), ...
                            entity_markers_vertices, entity_markers_values);
        end
    end
end

function marker_names = extract_marker_names(phys_content_str, output)

    n_outs = length(output);
    marker_names = cell(1, n_outs);
    if isempty(phys_content_str)
        return;
    else
        % Separate string content
        fmt = @(str) textscan(str, '%f %f %q');
        phys_content_cell = cellfun(fmt, phys_content_str, ...
                                    'UniformOutput', false);
        phys_content_cell = vertcat(phys_content_cell{:});

        % Extract markers
        entity_marker = cell2mat(phys_content_cell(:, 1));
        physical_marker = cell2mat(phys_content_cell(:, 2));

        % Extract names
        physical_names = phys_content_cell(:, 3);
        physical_names = vertcat(physical_names{:});

        % Build return
        for ml = 1:n_outs
            marker_req = entity_marker == output(ml);
            marker_names{ml} = {physical_marker(marker_req), ...
                                physical_names(marker_req)};
        end
    end
end

function entity_markers = compute_entity_markers(mesh, ...
                            entity_dim, entity_markers_vertices, ...
                            entity_markers_values)
    assert(entity_dim < mesh.dim);
    if entity_dim == 0
        % Special case of vertex markers
        entity_markers_entities = entity_markers_vertices;
    else
        % Get entity-to-vertex connectivity
        w = warning('off', 'Mesh:ExtraConnStored');  % Mute warning
        mesh.compute_connectivity(entity_dim, 0);
        warning(w);  % Restore warning
        entities_to_vertices = mesh.get_connectivity(entity_dim, 0);

        % Compute entity indices of marked entities
        [~, entity_markers_entities] = ismember(entity_markers_vertices, ...
                                                entities_to_vertices.', 'rows');
        assert(~any(~entity_markers_entities), ...
            'Some entitie-to-vertex relations in file are not part of the mesh.');
    end

    % Build sparse vector representation
    num_entities = mesh.num_entities(entity_dim);
    entity_markers = sparse(entity_markers_entities, 1, entity_markers_values, ...
                            num_entities, 1);
end


function [ele_content, vtx_content, phys_content] = read_sections(filename)
    % Get whole content from file.
    assert(ischar(filename), 'expected a filename');
    file_content = fileread(filename);

    % Store each line of the input file in a seperate cell.
    file_content = textscan(file_content, '%s', 'delimiter', '\n', ...
        'whitespace', '');
    file_content = file_content{1};

    % Check file.
    assert(~isempty(file_content), 'empty file detected');

    % Go through content and search for leading keywords which separates
    % information blocks (see gmsh file format documentation).
    name_list = {'$MeshFormat', '$EndMeshFormat', ...
                 '$PhysicalNames', '$EndPhysicalNames', ...
                 '$Nodes', '$EndNodes', ...
                 '$Elements', '$EndElements'};
    [~, name_in_file] = ismember(file_content, name_list);
    [idx_in_file, ~, idx_in_name] = find(name_in_file);
    assert(all(ismember([1:2, 5:8], idx_in_name)), ...
        ['Could not observe relevant name tags in the provided ', ...
        'file. Check if input is a Gmsh .msh file (given in ASCII ', ...
        'format version 2']);
    idx_format_start = idx_in_file(idx_in_name == 1) + 1;
    idx_format_stop = idx_in_file(idx_in_name == 2) - 1;
    idx_node_start = idx_in_file(idx_in_name == 5) + 2;
    idx_node_stop = idx_in_file(idx_in_name == 6) - 1;
    idx_element_start = idx_in_file(idx_in_name == 7) + 2;
    idx_element_stop = idx_in_file(idx_in_name == 8) - 1;
    idx_phys_start = idx_in_file(idx_in_name == 3) + 2;
    idx_phys_stop = idx_in_file(idx_in_name == 4) - 1;

    % Check if file type is supported.
    format_content = file_content(idx_format_start:idx_format_stop);
    format_content = strsplit(format_content{:});
    if sscanf(format_content{1}, '%d') ~= 2 ...
        error(sprintf(['MSH file format version %s.x detected. ', ...
            'Only MSH file format version 2.x is supported yet.'], ...
            format_content{1})); %#ok
    end

    %% Get element information.

    % Information structure:
    % $
    % Number of ele.
    % ele. number, ele. type, number of ele. tags, [tag list], [vertex list]
    % $End
    %
    % Note:
    %   [tag list] as many columns/entries as given by number of ele. tags
    %   [vertex list] differs in length for edges and cells
    %
    % ele. type == 1  -> edge
    %           == 2  -> triangle
    %           == 4  -> tetrahedron
    %           == 15 -> point
    %
    %    1. tag == physical entity id
    %                -> phys. line / phys. surface / phys. volume
    %    2. tag == elementary geometrical entity id
    %                -> straight line / plane surface id / volume id
    ele_content_num = file_content(idx_element_start - 1);
    ele_content_num = str2double(ele_content_num{:});
    ele_content_tmp = file_content(idx_element_start:idx_element_stop);
    n_ele_cells = length(ele_content_tmp);
    ele_content = cell(n_ele_cells, 1);
    ele_content_tmp = char(ele_content_tmp);
    for ii = 1:n_ele_cells
       ele_content{ii} = sscanf(ele_content_tmp(ii,:), '%f').';
    end
    assert(size(ele_content, 1) == ele_content_num, ...
        'Mismatch between obtained and expected number of elements.');

    %% Get vertex information.

    % Information structure:
    % $
    % Number of vertices
    % vtx. number, x-coord, y-coord, z-coord
    % $End
    vtx_content_num = file_content(idx_node_start - 1);
    vtx_content_num = str2double(vtx_content_num{:});
    vtx_content = file_content(idx_node_start:idx_node_stop);
    assert(size(vtx_content, 1) == vtx_content_num, ...
        'Mismatch between obtained and expected number of vertices.');

    %% Get physical entity information.

    % Information structure:
    % $
    % Number of phys. names
    % phys. dimension, phys. number, phys. name
    % $End
    % phys. dimension == 0 -> point
    %                 == 1 -> line (edge)
    %                 == 2 -> face (triangle)
    %                 == 3 -> volume (tetrahedron)
    % phys. number == physical entity id (-> see above)
    if ~isempty(idx_phys_start)
        phys_content_num = file_content(idx_phys_start - 1);
        phys_content_num = str2double(phys_content_num{:});
        phys_content = file_content(idx_phys_start:idx_phys_stop);
        assert(size(phys_content, 1) == phys_content_num, ...
            'Mismatch between obtained and expected number of physical entities.');
    else
        phys_content = {};
    end
end


function vertex_coords = extract_vertex_coordinates(vtx_content)
    % Use a trick to speed up transforming the strings within the cells
    % into an array of numbers (as number of rows and cols don't change).
    cols = size(str2double(strsplit(vtx_content{1})), 2);
    rows = size(vtx_content, 1);
    vtx_content = sprintf('%s ', vtx_content{:});
    vertex_coords = reshape(sscanf(vtx_content, '%f'), cols, rows).';
    vertex_coords(:, 1) = [];
end


function vertices = prune_z_component(vertices)
    tol = 10;
    n_vtx = size(vertices, 1);
    nonzero_cols = ~(sum(abs(vertices), 1) < tol * eps * n_vtx);
    vertices = vertices(:, nonzero_cols);
end


function topology_data = extract_topology_data(ele_content)
    % FIXME: We could avoid cell arrays and preallocate mats accordingly

    % Get types of all entities
    ele_types = cellfun(@(x) x(2), ele_content);

    % Sanity check
    supported_ele_types = [1, 2, 4, 15];
    found_ele_types = unique(ele_types);
    assert(all(ismember(found_ele_types, supported_ele_types)), ...
        sprintf(['File contains unsupported element types.\n', ...
        'Currently supported:', ...
        '\n id \t type', ...
        '\n 1 \t edge', ...
        '\n 2 \t triangle', ...
        '\n 4 \t tetrahedron', ...
        '\n 15 \t point']));
    assert(any(ismember(found_ele_types, [2, 4])), ...
        ['File solely contains point and/or edge information. ', ...
        'Only 2D meshes (containing triangle-elements) or 3D meshes ', ...
        '(containing thetrahedra-elements) are supported.']);

    % Extract entities; store with index offset 1
    topology_data{0+1} = cell2mat(ele_content(ele_types == 15));
    topology_data{1+1} = cell2mat(ele_content(ele_types == 1));
    topology_data{2+1} = cell2mat(ele_content(ele_types == 2));
    topology_data{3+1} = cell2mat(ele_content(ele_types == 4));
end


function [dim, cells] = extract_topology(topology_data)
    dim = find(cellfun(@(data) size(data, 2) > 0, topology_data), ...
               1, 'last') - 1;
    cell_data = topology_data{dim+1};
    cells = cell_data(:, end-dim:end);
    cells = uint32(cells);
end


function cell_markers = extract_cell_markers(dim, topology_data)
    cell_data = topology_data{dim+1};
    cell_markers = cell_data(:, 4);
end


function [entity_markers_values, entity_markers_vertices] = ...
            extract_entity_markers(dim, entity_dim, topology_data)
    assert(entity_dim < dim);
    entity_data = topology_data{entity_dim+1};

    % Create dummy columns in the corner case with no data
    if size(entity_data, 2) == 0
        entity_data = reshape(entity_data, 0, 5+entity_dim);
    end

    entity_markers_values = entity_data(:, 4);
    entity_markers_vertices = sort(entity_data(:, end-entity_dim:end), 2);
end
