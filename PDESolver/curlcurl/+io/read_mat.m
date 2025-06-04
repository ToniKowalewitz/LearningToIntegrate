function [mesh, varargout] = read_mat(filename, output)
    % Read mesh info from .mat file (assuming FORMAT data structure).
    %
    % SYNTAX
    %   [mesh, [varargout]] = read_mat(filename, [output])
    %
    % INPUT PARAMETER
    %   filename ... Char, filename of .msh file
    %   output   ... Vector, elements denoting output entity dimension:
    %                   -2 raw data
    %                   -1 physical names for markers
    %                    0 vertex   markers
    %                    1 line     markers
    %                    2 face     markers
    %                    3 volume   markers
    %
    % OUTPUT PARAMETERS
    %   mesh      ... object from Mesh class
    %   varargout ... entity markers
    %                    -2 -> Struct, containing .mat file info
    %                    -1 -> Cell {[marker].' , ['physical name'].'}
    %                     0 -> Vector of markers
    %                     1 -> Vector of markers
    %                     2 -> Vector of markers
    %                     3 -> Vector of markers
    %
    % REMARKS
    %
    %   Returned markers are zero for entities not marked.
    %   Cell markers is dense vector, other markers are sparse.
    %   Marker names are cells for each requested marker entity.
    %
    % FORMAT
    %   data from filename is assumed to be structured as:
    %
    %   req. fields:
    %   data.vtx ... Matrix [n x dim] of vertice coordinates.
    %   data.ele ... Matrix [m x dim+1] of vertex2tetrahedra connectivity.
    %
    %   opt. fields:
    %   data.selection may containing
    %       .ele      ... cell entities
    %       .vtx      ... point entities
    %       .vtx2line ... line entities
    %       .vtx2face ... face entities
    %   wherein
    %       .name.id   ... Vector or Matrix of entity tags
    %   and arbitrary additional fields e.g.
    %            .meta ... String with meta info
    %            .rho  ... Scalar of domain resistivities
    %   e.g.
    %       data.selectio.tet.domain_1.id = [354; 355; 365]
    %       data.selectio.tet.domain_1.rho = 50
    %       data.selectio.tet.domain_1.meta = 'first layer beneeth surface'

    if nargin < 2
        output = [];
    end

    % Set parameter.
    req_vars = {'vtx', 'ele'};
    opt_vars = {'selection'};
    opt_vars_field = {'vtx', 'vtx2line', 'vtx2face', 'ele'};

    %% Check input.

    % Load file.
    data = load(filename);
    if ~all(isfield(data, req_vars))
        error('Incomplete mesh info, check for required fields.');
    end
    if ~isfield(data, opt_vars)
        warning('No entity marker available.');
        skip_marker = true;
    else
        skip_marker = false;
    end
    % FIXME: true, or are coords stored in three components anyway?
    dim = size(data.vtx, 2);
    assert(all(dim <= 3 & size(data.ele, 2) <= dim+1), ...
        'Check dimensions of required fields');

    %% Build mesh object.

    mesh = meshing.Mesh(dim, data.vtx.', data.ele.');

    %% Extract mesh topology.

    if skip_marker
        % Skip.
    else
        if any(dim < output)
           error('Requested marker exceed dimensionality of mesh.');
        elseif any(output < -2)
           error('Some outputs are not supported.');
        end

        % Extract marker names and ids.
        marker_names = extract_marker_names(output, data, opt_vars_field);

        % Prepare and check variable output
        n_outs = length(output);
        nargoutchk(n_outs+1, n_outs+1);
        varargout = cell(1, n_outs);
        for mm = 1:n_outs
            switch output(mm)
                case -2
                % Extract raw data.
                varargout{mm} = data;

                case -1
                % Extract marker names.
                varargout{mm} = marker_names;

                case {0, 3}
                % Extract vertex or cell ids.
                cur_info = marker_names{mm-1}{2};
                if ~isempty(cur_info)
                    varargout{mm} = extract_entity_marker(mesh, ...
                        data.(opt_vars{:}).(opt_vars_field{output(mm)+1}), ...
                        cur_info, output(mm));
                end

                case {1, 2}
                % Extract line or face ids.
                cur_info = marker_names{mm-1}{2};
                if ~isempty(cur_info)
                    varargout{mm} = map_entity_marker(mesh, ...
                        data.(opt_vars{:}).(opt_vars_field{output(mm)+1}), ...
                        cur_info, output(mm));
                end
            end
        end
    end
end

function mn = extract_marker_names(output, data, names)
    % OUTPUT
    %   mn ... Cell [l x 1] with l = length(output)-1
    % each cell containing a cell [1 x 2] with:
    %   {1} vector [n x 1] with n = 1 ... Scalar, id
    %   {2} cell [n x 2]   with n = 1 ... String, name

    % Extract only entity marker requests.
    outs = output >= 0;
    outs = output(outs);
    n_outs = length(outs);

    % Initialize marker_names cell.
    ds = data.selection;
    mn = cell(1, n_outs);
    for ii = 1:n_outs
        cur_string = names{outs(ii)+1};
        mn{ii} = cell(1, 2);
        mn{ii}{2} = create_marker_cell();
        mn{ii}{1} = (1:length(mn{ii}{2})).';
    end

    % Shortcut.
    function mc = create_marker_cell()
        % Extract field name string and connect with ids.
        if isfield(ds, cur_string)
            cur_struct = ds.(cur_string);
            field_names = fieldnames(cur_struct);
            mc = field_names;
        else
            mc = {};
        end
    end
end

function em = map_entity_marker(mesh, data, info, dim)
    % OUTPUT
    %   em ... Sparse vector [mesh.num_entities(dim) x 1] of entity ids.
    %   wherein
    %       id ... ongoing number w.r.t. line number in info

    % Get connectivities.
    w = warning('off', 'Mesh:ExtraConnStored');  % Mute warning
    mesh.compute_connectivity(dim, 0);
    warning(w);  % Restore warning
    entities_to_vertices = mesh.get_connectivity(dim, 0);
    num_entities = mesh.num_entities(dim);

    [marker_map, marker_vals] = deal([]);
    for ii = 1:size(info, 1)
        % Check and fetch data.
        if ~isfield(data.(info{ii}), 'id')
            error('Entity id missing.');
        end
        cur_data = data.(info{ii}).id;
        assert(size(cur_data, 2) == dim+1);
        cur_entity_markers_vertices = sort(cur_data, 2); % ensure ascending numbering
        cur_entity_markers_values = ii + zeros(size(cur_entity_markers_vertices, 1), 1);

        % Compute entity indices of marked entities
        [~, cur_entity_marker_map] = ismember(cur_entity_markers_vertices, ...
                                              entities_to_vertices.', 'rows');
        assert(~any(~cur_entity_marker_map), ...
            'Some entity-to-vertex relations in file are not part of the mesh.');
        marker_map = [marker_map; cur_entity_marker_map]; %#ok<AGROW>
        marker_vals = [marker_vals; cur_entity_markers_values]; %#ok<AGROW>
    end
    % Build sparse vector representation
    em = sparse(marker_map, 1, marker_vals, num_entities, 1);
end

function em = extract_entity_marker(mesh, data, info, dim)
    % OUTPUT
    %   em ... Vector [mesh.num_entities(dim) x 1] of entity ids.
    %   wherein
    %       id ... ongoing number w.r.t. line number in info

    % Get connectivities and init entity vector.
    em = zeros(mesh.num_entities(dim), 1);
    for ii = 1:size(info, 1)
        % Check and fetch data.
        if ~isfield(data.(info{ii}), 'id')
            error('Entity id missing.');
        end
        cur_data = data.(info{ii}).id;
        assert(size(cur_data, 2) == 1);
        em(cur_data) = ii;
    end
    % FIXME: implement more elegant.
    if dim == 0
        em = sparse(em);
    end
end
