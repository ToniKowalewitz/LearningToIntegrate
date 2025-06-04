classdef XDMF < handle

    methods (Access = public)

        function obj = XDMF(filename, encoding)
            % Create XDMF output handler for output of FE functions
            %
            % INPUT ARGS
            %   filename ... relative or absolute path of XDMF file (will be overwritten)
            %   encoding ... optional; 'binary' (default), 'xml', or 'hdf5'

            obj.filename = filename;

            if nargin < 2
                encoding = 'binary';
            end
            obj.encoding = lower(encoding);

            xdmf = create_xdmf_node();
            domain = add_domain_node(xdmf);
            obj.gc = add_grid_collection_node(domain, 'TimeSeries');

        end

        function write(obj, dofmap, x, t)
            % Write finite element function to XDMF
            %
            % INPUT PARAMETERS
            %   dofmap ... Struct, FE function space
            %   x      ... Vector [dofmap.dim, 1], DOFs of FE function
            %   t      ... Real scalar, time or index
            %
            % NOTES
            %   Currently it is only possible to use a matching
            %   dofmap in subsequent writes, which is enforced.
            %
            %   File contents may be cached; XDMF.flush() can be
            %   used when needed.

            obj.write_func(dofmap, x, t)

        end

        function flush(obj)
            % Dump possibly cached contents to disk

            xmlwrite(obj.filename, obj.gc.getOwnerDocument());

        end

        function delete(obj)

            obj.flush();

        end

    end

    properties (SetAccess = immutable)

        filename;
        encoding;

    end

    properties (SetAccess = private)

        counter_grid = 0;
        counter_file = 0;
        gc;
        dofmap_ptr = 0;
        embedding_get_cell_dofs;
        embedding_get_x;
        first_celldofs_node_;

    end

    methods (Access = private)

        function write_func(obj, dofmap, x, t)

            grid = add_grid_node(obj.gc, obj.get_next_grid_name());
            add_time_node(grid, t);
            attr = add_attribute_node(grid, dofmap.element, 'u');

            if obj.dofmap_ptr == 0
                obj.write_func_first(grid, attr, dofmap, x)
                obj.dofmap_ptr = util.getPr(dofmap);
            else
                assert(obj.dofmap_ptr == util.getPr(dofmap));
                obj.write_func_second(grid, attr, x)
            end

        end

        function write_func_first(obj, grid, attr, dofmap, x)

            % Prepare mesh topology and geometry XDMF nodes
            [num_vertices_per_cell, num_cells] = size(dofmap.mesh.cells);
            topology = add_topology_node(grid, num_cells, num_vertices_per_cell);
            geometry = add_geometry_node(grid, dofmap.mesh.dim);

            % Write out heavy data for mesh
            obj.add_data_item(topology, dofmap.mesh.cells-1);  % convert to zero-based indices
            obj.add_data_item(geometry, dofmap.mesh.vertex_coords);

            % Prepare operators for embedding into XDMF-supported function space
            [obj.embedding_get_cell_dofs, obj.embedding_get_x] = ...
                prepare_embedding(dofmap);

            % Write out heavy data for dofmap
            celldofs = obj.add_data_item(attr, obj.embedding_get_cell_dofs());
            celldofs.setAttribute('Name', 'celldofs');

            % Write out heavy data for function values
            obj.add_data_item(attr, obj.embedding_get_x(x));

            % Remember celldofs node (if appropriate)
            if strcmp(obj.encoding, 'binary')
                obj.first_celldofs_node_ = celldofs;
            end

        end

        function write_func_second(obj, grid, attr, x)

            % Write mesh pointer to the first mesh
            add_topology_geometry_ptr(grid, 'TimeSeries', ...
                                      obj.get_first_grid_name());

            % Write function dofmap
            switch obj.encoding
            case 'xml'
                % NB: None of the following pointers work in ParaView 5.8
                %add_celldofs_ptr1(attr, 'TimeSeries', obj.get_first_grid_name(), 'u');
                %add_celldofs_ptr2(attr, 'TimeSeries', obj.get_first_grid_name(), 'u');
                %add_celldofs_ptr3(attr, 'TimeSeries', obj.get_first_grid_name(), 'u');
                celldofs = obj.add_data_item(attr, obj.embedding_get_cell_dofs());
                celldofs.setAttribute('Name', 'celldofs');
            case 'binary'
                % NB: Notice cloneNode(true), one cannot reuse a node in XML tree
                attr.appendChild(obj.first_celldofs_node_.cloneNode(true));
            case 'hdf5'
                error('not implemented');
            end

            % Write function values
            obj.add_data_item(attr, obj.embedding_get_x(x));

        end

        function child = add_data_item(obj, parent, array)
            % Add 'DataItem' entry to XDMF XML and write out heavy data (array)

            switch obj.encoding
            case 'xml'

                % Create XML entry and dump data as XML text node
                child = add_data_item_xml(parent, array);

            case 'binary'

                % Construct filename for binary file
                [fpath, fname] = fileparts(obj.filename);
                counter = obj.get_next_file_counter();
                binfilename = fullfile(fpath, strcat(fname, counter, '.bin'));

                % Create XML entry and dump data to binary file
                child = add_data_item_binary(parent, array, binfilename);

            case 'hdf5'
                error('encoding "%s" not yet implemented', obj.encoding);
            otherwise
                error('unknown encoding "%s"', obj.encoding);
            end

        end

        function name = get_first_grid_name(obj)  %#ok<MANU>
            name = sprintf('grid%d', 0);
        end

        function name = get_next_grid_name(obj)
            name = sprintf('grid%d', obj.counter_grid);
            obj.counter_grid = obj.counter_grid + 1;
        end

        function counter = get_next_file_counter(obj)
            counter = sprintf('%04d', obj.counter_file);
            obj.counter_file = obj.counter_file + 1;
        end

    end

end


function child = add_data_item_xml(parent, array)

    doc = parent.getOwnerDocument();

    % Prepare data item node
    child = create_data_item_node(doc, array);
    child.setAttribute('Format', 'XML');

    % Prepare char array with heavy data (!)
    switch class(array)
    case 'double'
        fmt = strcat(repmat('%.17g ', 1, size(array, 1)), '\n');
    case 'uint32'
        fmt = strcat(repmat('%d ', 1, size(array, 1)), '\n');
    otherwise
        error('not implemented');
    end
    data = sprintf(fmt, array);

    % Add text node (with heavy data!) to data item node
    child.appendChild(doc.createTextNode(data));
    parent.appendChild(child);

end


function child = add_data_item_binary(parent, array, filename)

    doc = parent.getOwnerDocument();

    % Prepare data item node
    child = create_data_item_node(doc, array);
    child.setAttribute('Format', 'Binary');

    % Write out array to binary file
    [fid, msg] = fopen(filename, 'wb');
    if fid == -1
        error('Opening "%s" for writing failed with reason:\n%s', filename, msg);
    end
    count = fwrite(fid, array, class(array));
    assert(count == numel(array));
    fclose(fid);

    % Add text node (with filename) to data item node
    child.appendChild(doc.createTextNode(filename));
    parent.appendChild(child);

end


function child = create_data_item_node(doc, array)

    assert(ismatrix(array));
    [dim1, dim2] = size(array);

    switch class(array)
    case 'double'
        scalar_type = 'Float';
        scalar_size = '8';
    case 'uint32'
        scalar_type = 'UInt';
        scalar_size = '4';
    otherwise
        error('not implemented');
    end

    child = doc.createElement('DataItem');
    child.setAttribute('Dimensions', sprintf('%d %d', dim2, dim1));
    child.setAttribute('NumberType', scalar_type);
    child.setAttribute('Precision', scalar_size);

end


function xdmf = create_xdmf_node()

    doc = com.mathworks.xml.XMLUtils.createDocument('Xdmf');  %#ok<*JAPIMATHWORKS>

    xdmf = doc.getDocumentElement();
    xdmf.setAttribute('Version', '3.0');
    xdmf.setAttribute('xmlns:xi', 'http://www.w3.org/2001/XInclude');

end


function child = add_domain_node(parent)

    doc = parent.getOwnerDocument();

    child = doc.createElement('Domain');
    parent.appendChild(child);

end


function child = add_grid_collection_node(parent, name)

    doc = parent.getOwnerDocument();

    child = doc.createElement('Grid');
    child.setAttribute('GridType', 'Collection');
    child.setAttribute('CollectionType', 'Temporal');
    child.setAttribute('Name', name);
    parent.appendChild(child);

end


function child = add_time_node(parent, t)

    doc = parent.getOwnerDocument();

    child = doc.createElement('Time');
    child.setAttribute('Value', sprintf('%.17g', t));
    parent.appendChild(child);

end


function child = add_topology_geometry_ptr(parent, grid_collection_name, grid_name)

    doc = parent.getOwnerDocument();

    child = doc.createElement('xi:include');
    parent.appendChild(child);

    fmt = 'xpointer(//Grid[@Name="%s"]/Grid[@Name="%s"]/*[self::Topology or self::Geometry])';
    ptr = sprintf(fmt, grid_collection_name, grid_name);
    child.setAttribute('xpointer', ptr);

end


function child = add_grid_node(parent, name)

    doc = parent.getOwnerDocument();

    child = doc.createElement('Grid');
    child.setAttribute('Name', name);
    child.setAttribute('GridType', 'Uniform');
    parent.appendChild(child);

end


function child = add_topology_node(parent, num_cells, num_vertices_per_cell)

    doc = parent.getOwnerDocument();

    switch num_vertices_per_cell
    case 3
        topology_type = 'Triangle';
    case 4
        topology_type = 'Tetrahedron';
    otherwise
        error('not impl');
    end

    child = doc.createElement('Topology');
    child.setAttribute('NumberOfElements', num2str(num_cells));
    child.setAttribute('TopologyType', topology_type);
    child.setAttribute('NodesPerElement', num2str(num_vertices_per_cell));
    parent.appendChild(child);

end


function child = add_geometry_node(parent, dim)

    doc = parent.getOwnerDocument();

    switch dim
    case 2
        geometry_type = 'XY';
    case 3
        geometry_type = 'XYZ';
    otherwise
        error('not impl');
    end

    child = doc.createElement('Geometry');
    child.setAttribute('GeometryType', geometry_type);
    parent.appendChild(child);

end


function child = add_attribute_node(parent, element, name)

    doc = parent.getOwnerDocument();

    dim = element.simplex.dim;
    order = element.order;

    child = doc.createElement('Attribute');
    child.setAttribute('ItemType', 'FiniteElementFunction');
    child.setAttribute('ElementDegree', num2str(order));
    if dim == 2
        child.setAttribute('ElementCell', 'triangle');
    elseif dim == 3
        child.setAttribute('ElementCell', 'tetrahedron');
    else
        error('not implemented');
    end

    switch element.family
    case 'NedelecKind1'
        child.setAttribute('ElementFamily', 'DG');
        child.setAttribute('Center', 'Other');
        child.setAttribute('AttributeType', 'Vector');
    case 'Lagrange'
        % TODO: Want special code path for P1, avoiding FiniteElementFunction?
        child.setAttribute('ElementFamily', 'CG');
        child.setAttribute('Center', 'Other');
        child.setAttribute('AttributeType', 'Scalar');
    case 'DiscontinuousLagrange'
        % TODO: Want special code path for P0, avoiding FiniteElementFunction?
        child.setAttribute('ElementFamily', 'DG');
        child.setAttribute('Center', 'Other');
        child.setAttribute('AttributeType', 'Scalar');
    otherwise
        error('not implemented');
    end

    child.setAttribute('Name', name);
    parent.appendChild(child);

end


function child = add_celldofs_ptr1(parent, grid_collection_name, grid_name, attribute_name)  %#ok<DEFNU>

    doc = parent.getOwnerDocument();

    child = doc.createElement('xi:include');
    parent.appendChild(child);

    fmt = 'xpointer(//Grid[@Name="%s"]/Grid[@Name="%s"]/Attribute[@Name="%s"]/DataItem[@Name="celldofs"])';
    ptr = sprintf(fmt, grid_collection_name, grid_name, attribute_name);
    child.setAttribute('xpointer', ptr);

end


function child = add_celldofs_ptr2(parent, grid_collection_name, grid_name, attribute_name)  %#ok<DEFNU>

    doc = parent.getOwnerDocument();

    child = doc.createElement('DataItem');
    parent.appendChild(child);

    fmt = '/Xdmf/Domain/Grid[@Name="%s"]/Grid[@Name="%s"]/Attribute[@Name="%s"]/DataItem[@Name="celldofs"]';
    ptr = sprintf(fmt, grid_collection_name, grid_name, attribute_name);
    child.setAttribute('Reference', ptr);

end


function child = add_celldofs_ptr3(parent, grid_collection_name, grid_name, attribute_name)  %#ok<DEFNU>

    doc = parent.getOwnerDocument();

    child = doc.createElement('DataItem');
    parent.appendChild(child);

    fmt = '/Xdmf/Domain/Grid[@Name="%s"]/Grid[@Name="%s"]/Attribute[@Name="%s"]/DataItem[@Name="celldofs"]';
    ptr = sprintf(fmt, grid_collection_name, grid_name, attribute_name);
    child.setAttribute('Reference', 'XML');
    child.appendChild(doc.createTextNode(ptr));

end


function [get_cell_dofs, get_x] = prepare_embedding(dofmap)
    switch dofmap.element.family
    case 'NedelecKind1'
        [get_cell_dofs, get_x] = prepare_embedding_nedelec(dofmap);
    case {'Lagrange', 'DiscontinuousLagrange'}
        [get_cell_dofs, get_x] = prepare_embedding_lagrange(dofmap);
    otherwise
        error('Not implemented');
    end
end


function [get_cell_dofs, get_x] = prepare_embedding_nedelec(dofmap)
    % Embedding of Nedelec into (vector-valued) discontinuous Lagrange

    % Fetch some data
    dim = dofmap.element.simplex.dim;
    assert(dim == dofmap.mesh.dim);
    assert(strcmp(dofmap.element.family, 'NedelecKind1'));
    num_cells = dofmap.mesh.num_entities(dim);
    cells = dofmap.mesh.cells;
    cell_dofs = dofmap.cell_dofs;
    coords = dofmap.mesh.vertex_coords;

    % Target FE space is vector DG
    order = dofmap.element.order;
    element = fe.create_discontinuous_lagrange_element(dim, order);
    num_dofs_per_cell = element.fe_space_dim;
    num_vector_dofs_per_cell = dim*element.fe_space_dim;

    % Create local embedding matrix (and transpose for later convenience)
    M = create_transfer_matrix(dofmap.element, element);
    M = reshape(permute(M, [3, 2, 1]), [], dofmap.element.fe_space_dim);

    % NB: The following two function could operate on a single cell,
    %     thus saving a lot of memory for the prize of small slow down...

    function y = embed_nedelec_in_dg(x)
        % Compute discontinuous Lagrange function (given by y)
        % for the embedded Nedelec function (given by x)

        % Allocate vector in target space
        y = zeros(num_cells*num_vector_dofs_per_cell, 1);

        for c = 1:num_cells

            % Populate geometry for pullback
            B = coords(:, cells(1:dim, c)) - coords(:, cells(dim+1, c));

            % Get local dofs
            x_cell = x(cell_dofs(:, c));

            % Interpolate
            % FIXME: Put this complicated (reshape) logic into a (testable) function?!?
            temp = reshape(M*x_cell, num_dofs_per_cell, dim);
            cell_coeffs = reshape(temp/B, num_vector_dofs_per_cell, 1);

            % Store in global vector
            dof0 = num_vector_dofs_per_cell*(c-1);
            y(dof0+1:dof0+num_vector_dofs_per_cell) = cell_coeffs;

        end

        % Ensure row vector
        y = y(:).';

    end

    % How to permute our DG dofs into XDMF DG dofs
    perm = get_xdmf_dof_permutation(element);  % permutation for scalar element
    offsets = (1:num_dofs_per_cell:num_vector_dofs_per_cell) - 1;
    perm = reshape(perm(:) + offsets(:).', 1, num_vector_dofs_per_cell);  % permutation for all components

    % Check that the result is a permutation
    assert(all(sort(perm) == 1:num_vector_dofs_per_cell));

    function cell_dofs = get_cell_dofs_()
        % Forge DG dofmap

        % Create range and reshape
        num_dofs = num_cells*num_vector_dofs_per_cell;
        cell_dofs = uint32(0):uint32(num_dofs-1);  % zero-based indexing
        cell_dofs = reshape(cell_dofs, num_vector_dofs_per_cell, num_cells);

        % Permute from our dof ordering to XDMF dof ordering
        cell_dofs = cell_dofs(perm, :);

        % Check we have correct type
        assert(isa(cell_dofs, 'uint32'));
    end

    get_cell_dofs = @get_cell_dofs_;
    get_x = @embed_nedelec_in_dg;

end


function [get_cell_dofs, get_x] = prepare_embedding_lagrange(dofmap)
    % Embedding of our Lagrange into XDMF Lagrange

    % Lagrange and discontinuous Lagrange are directly represented in XDMF,
    % up to permutation of dofs and zero-based indices
    perm = get_xdmf_dof_permutation(dofmap.element);

    % Permute and convert to zero-based indices
    get_cell_dofs = @() dofmap.cell_dofs(perm, :) - 1;

    % Ensure row vector
    get_x = @(x) x(:).';

end


function perm = get_xdmf_dof_permutation(element)
    % Return permutation of dofs to get XDMF dof ordering
    %
    % The ordering is documented at [1] and is given
    % by FEniCS/UFC ordering [2]
    %
    % [1] https://www.xdmf.org/index.php/XDMF_Model_and_Format#Attribute
    % [2] https://fenicsproject.org/pub/documents/ufc/ufc-user-manual/ufc-user-manual.pdf#page=64

    % XDMF/Paraview only supports Lagrange and Raviart-Thomas,
    % which is not implemented yet
    assert(any(strcmp(element.family, {'Lagrange', 'DiscontinuousLagrange'})));

    % In this routine we work out the permutation for scalar element only
    assert(element.value_shape == 1);

    % HACK: Switch from discontinuous to Lagrange if needed;
    % DiscontinuousLagrange does not have proper entity dofs
    if element.order > 0 && strcmp(element.family, 'DiscontinuousLagrange')
        element = fe.create_lagrange_element(element.simplex.dim, element.order);
    end

    % Vertex dofs and edge dofs
    vdofs = element.get_entity_dofs(0);
    edofs = element.get_entity_dofs(1);

    % Function to compute edge index from its vertices
    v2e_ = element.simplex.get_connectivity(0, 1);
    edge = @(vertex1, vertex2) intersect(v2e_(:, vertex1), v2e_(:, vertex2));

    switch element.order
    case 0

        % Trivial permutation, only one dof
        perm = 1;

    case 1

        % Trivial permutation (vdofs is usually monotone)
        perm = vdofs;

    case 2
        switch element.simplex.dim
        case 1
            % Edge 12 is only edge
            perm = [vdofs(1), vdofs(2), ...
                    edofs(edge(1, 2))];
        case 2
            % Edge 23 is adjacent to vertex 1, etc.
            perm = [vdofs(1), vdofs(2), vdofs(3), ...
                    edofs(edge(2, 3)), edofs(edge(1, 3)), edofs(edge(1, 2))];
        case 3
            % Edge 34 is adjacent to vertices 1,2 (proceed lexicographically), etc.
            perm = [vdofs(1), vdofs(2), vdofs(3), vdofs(4), ...
                    edofs(edge(3, 4)), edofs(edge(2, 4)), edofs(edge(2, 3)), ...
                    edofs(edge(1, 4)), edofs(edge(1, 3)), edofs(edge(1, 2))];
        end

    otherwise
        error('element %s order %d not supported in XDMF', ...
              element.family, element.order);
    end

    % Check that the result is a permutation
    assert(all(sort(perm) == 1:element.fe_space_dim));

end


function M = create_transfer_matrix(origin_element, destination_element)
    % Cell-local interpolation/embedding operator between two elements

    o = origin_element;
    d = destination_element;

    dim = o.simplex.dim;
    assert(d.simplex.dim == dim);

    % Value size
    vs = dim;  % NB: Now assuming vector-valued elements

    basis_flat = @(x) reshape(o.tabulate_basis(x), o.fe_space_dim*vs, 1);
    M = d.evaluate_dual_basis(basis_flat);
    M = reshape(M, o.fe_space_dim, vs, d.fe_space_dim);

end
