function em = orient_wire(mesh, em, mv, mv_plus, mv_minus)
    % Partition marked edges according to orientation
    %
    % SYNTAX
    %
    %   em = orient_wire(mesh, em, mv, mv_plus, mv_minus)
    %
    % INPUT PARAMETERS
    %   mesh     ... Instance of Mesh class.
    %   em       ... Sparse vector marking mesh edges.
    %   mv       ... Marker value of marked edges.
    %   mv_plus  ... Desired marker value of positively oriented edges.
    %   mv_minus ... Desired marker value of negatively oriented edges.
    %
    % OUTPUT PARAMETER
    %   em       ... Sparse vector marking partitioned edges.
    %
    % NOTES
    %   Relation (em == mv) has to mark edges which form
    %   either a wire or a loop. If tree is found, i.e.,
    %   a vertex connected to three or more edges,
    %   error is raised. If marked edges consists of
    %   more than one disconnected components, error
    %   is raised as well.

    e2v = mesh.get_connectivity(1, 0);

    % Get edge vertex connectivity on the wire
    edge_vertices = e2v(:, em == mv);
    num_edges = size(edge_vertices, 2);
    assert(size(edge_vertices, 1) == 2);

    % No edges, return zeros
    if num_edges == 0
        em = sparse([], [], [], numel(em), 1, 0);
        return
    end

    % Sort wire vertex indices, remember sorting permutation
    [edge_vertices, perm] = sort(edge_vertices(:));
    invperm = invert_perm(perm);

    % Find endpoint (any if closed loop) vertex and corresponding edge
    v0 = find_first_vert(edge_vertices);
    e0 = vert_to_edge(num_edges, perm, v0);

    % Initialize orientations, orient first edge positively
    orient = zeros(num_edges, 1);
    orient(e0) = 1;

    % Loop until end of wire or loop around
    v2 = v0;
    while true

        % Move to adjacent vertex of the current edge
        v1 = get_other_edge_vert(num_edges, perm, invperm, v2);

        % Find same vertex on adjacent edge
        v2 = get_connected_vert(edge_vertices, v1);

        % We reached the end of the wire or looped around
        if v2 == 0 || v2 == v0
            break
        end

        % Compute edges indices
        [e1, lv1] = vert_to_edge(num_edges, perm, v1);
        [e2, lv2] = vert_to_edge(num_edges, perm, v2);

        % Figure out relative orientation of two edges
        relative_orientation = (-1)^(lv1 + lv2 + 1);

        % Extend orientation from already oriented to unoriented
        if orient(e1) ~= 0
            assert(orient(e2) == 0);
            orient(e2) = relative_orientation * orient(e1);
        elseif orient(e2) ~= 0
            orient(e1) = relative_orientation * orient(e2);
        else
            error('no edge oriented for some vertex! wire has multiple components?');
        end
    end

    % Something went wrong, not all edges oriented
    if nnz(orient) < num_edges
        error('not all edges oriented! wire has multiple components or is it tree?');
    end

    % Partition edge markers by orientations
    inds = find(em == mv);
    em = sparse([], [], [], numel(em), 1, numel(orient));
    em(inds(orient == +1)) = mv_plus;
    em(inds(orient == -1)) = mv_minus;
end


function inverse = invert_perm(perm)
    % Invert permutation

    inverse = zeros(numel(perm), 1);
    for j = 1:numel(perm)
        inverse(perm(j)) = j;
    end
end


function j = get_connected_vert(verts, j)
    % For a vertex given by index j (indexing sorted vertices)
    % find the second index of the same vertex (in the adjacent
    % edge)

    % Previous vertex
    if j > 1
        a = verts(j-1);
    else
        a = 0;
    end

    % Current vertex
    b = verts(j);

    % Following vertex
    if j < numel(verts)
        c = verts(j+1);
    else
        c = 0;
    end

    if a == b
        % Predecessor is connected
        assert(b ~= c);
        j = j-1;
    elseif b == c
        % Successor is connected
        j = j+1;
    else
        % No connected vertex
        j = 0;
    end
end


function [e, lv] = vert_to_edge(num_edges, perm, j)
    % Convert a vertex given by index j (indexing sorted vertices)
    % into (local_vertex, edge)-pair

    % Convert to flattened index into edge-vertex topology
    j = perm(j);

    % Convert to unflattened pair of (local_vertex, edge)
    [lv, e] = ind2sub([2, num_edges], j);
end


function j = get_other_edge_vert(num_edges, perm, invperm, j)
    % For a vertex given by index j (indexing sorted vertices)
    % find the other vertex on the same edge given by
    % index j (indexing sorted vertices)

    % Convert to flattened index into edge-vertex topology
    j = perm(j);

    % Convert to unflattened pair of (local_vertex, edge)
    [lv, e] = ind2sub([2, num_edges], j);

    % Pick the other vertex on the edge
    if lv == 1
        lv = 2;
    else
        lv = 1;
    end

    % Convert back to flattened index into edge-vertex topology
    j = sub2ind([2, num_edges], lv, e);

    % Convert to index into sorted vertices
    j = invperm(j);
end


function j = find_first_vert(verts)
    % Find endpoint vertex of wire (or any if closed loop)
    %
    % Input is sorted (!) list of vertex indices defining
    % an open-ended wire or a closed loop. Returns an index
    % to vertex which is endpoint or any if loop. Raises
    % error if vertex of degree > 2 is found.

    % Chech first two vertices
    for j = 1
        b = verts(j);
        c = verts(j+1);

        % Success first two vertices different
        if b < c
            return
        end
    end

    % Check all of successive 3 vertices
    for j = 2:numel(verts)-1
        a = verts(j-1);
        b = verts(j);
        c = verts(j+1);

        % Success, all 3 vertices different
        if a < b && b < c
            return
        end

        % Error, all 3 vertices same
        if a == b && b == c
            error('vertex of degree >= 3');
        end
    end

    % Check last two vertices
    for j = numel(verts)
        a = verts(j-1);
        b = verts(j);

        % Success last two vertices different
        if a < b
            return
        end
    end

    % Closed loop, return arbitrary
    j = 1;
end
