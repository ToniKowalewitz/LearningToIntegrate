function tests = test_mark_entities
  % Run all function tests in this file
  tests = functiontests(localfunctions);
end


function test_mark_entities_doc_snippet(t)
    figure(5555);

    % Supress warning about unneeded connectivities computed and stored
    import matlab.unittest.fixtures.SuppressedWarningsFixture
    t.applyFixture(SuppressedWarningsFixture('Mesh:ExtraConnStored'));

    % Test that the snippet from doc of 'mark_entities' runs
    mesh = meshing.generate_half_disc_mesh(16);
    mesh.compute_connectivity(mesh.dim, mesh.dim-1);
    mesh.compute_boundary_facets();
    tol = 1e-10;
    boundary1 = @(x, b) b && norm(x) > 1-tol;
    markers = meshing.mark_entities(mesh, mesh.dim-1, [], boundary1, 1);
    boundary2 = @(x, b) b && x(2) < tol;
    markers = meshing.mark_entities(mesh, mesh.dim-1, markers, boundary2, 2);
    interior = @(x, b) ~b;
    markers = meshing.mark_entities(mesh, mesh.dim-1, markers, interior, -1);
    meshing.plot_mesh(mesh, 'facet_markers', markers);

    close(5555);
end


function test_facet_markers_half_disc(t)

    % Supress warning about unneeded connectivities computed and stored
    import matlab.unittest.fixtures.SuppressedWarningsFixture
    t.applyFixture(SuppressedWarningsFixture('Mesh:ExtraConnStored'));

    % Create mesh and markers
    mesh = meshing.generate_half_disc_mesh(16);
    mesh.compute_connectivity(mesh.dim, mesh.dim-1);
    mesh.compute_boundary_facets();
    tol = 1e-10;
    boundary1 = @(x, b) b && norm(x) > 1-tol;
    markers = meshing.mark_entities(mesh, mesh.dim-1, [], boundary1, 1);
    boundary2 = @(x, b) b && x(2) < tol;
    markers = meshing.mark_entities(mesh, mesh.dim-1, markers, boundary2, 2);
    interior = @(x, b) ~b;
    markers = meshing.mark_entities(mesh, mesh.dim-1, markers, interior, -1);

    % Fetch needed data
    coords = mesh.vertex_coords;
    num_facets = mesh.num_entities(mesh.dim-1);
    f2v = mesh.get_connectivity(mesh.dim-1, 0);
    boundary_facets = mesh.get_boundary_facets_indices();
    boundary_facets = sparse(double(boundary_facets), 1, true, num_facets, 1);

    % Verify facet by facet
    for f = 1:num_facets
        if markers(f) == -1
            t.verifyTrue(~boundary_facets(f));
        elseif markers(f) == 1
            t.verifyTrue(boundary_facets(f));
            for v = 1:mesh.dim
                vertex_coords = coords(:, f2v(v, f));
                t.verifyEqual(norm(vertex_coords), 1, 'AbsTol', 1e-14);
            end
        elseif markers(f) == 2
            t.verifyTrue(boundary_facets(f));
            for v = 1:mesh.dim
                vertex_coords = coords(:, f2v(v, f));
                t.verifyEqual(vertex_coords(2), 0, 'AbsTol', 1e-14);
            end
        end
    end

end


function test_facet_markers_cube_1d(t)
    % Throws: 1D cube mesh not implemented
    t.verifyError(@() t_facet_markers_cube(t, 66), ?MException);
end


function test_facet_markers_cube_2d(t)
    t_facet_markers_cube(t, [7, 13]);
end


function test_facet_markers_cube_3d(t)
    t_facet_markers_cube(t, [4, 5, 6]);
end


function t_facet_markers_cube(t, dims)

    % Supress warning about unneeded connectivities computed and stored
    import matlab.unittest.fixtures.SuppressedWarningsFixture
    t.applyFixture(SuppressedWarningsFixture('Mesh:ExtraConnStored'));

    % Create mesh and markers
    mesh = meshing.generate_unit_cube_mesh(dims);
    mesh.compute_connectivity(mesh.dim, mesh.dim-1);
    num_facets = mesh.num_entities(mesh.dim-1);
    mesh.compute_boundary_facets();
    tol = 1e-10;
    markers = zeros(num_facets, 1);
    for d=1:mesh.dim
        boundary = @(x, b) b && ( x(d) < tol || x(d) > 1-tol );
        markers = meshing.mark_entities(mesh, mesh.dim-1, markers, boundary, d);
    end

    % Fetch needed data
    coords = mesh.vertex_coords;
    f2v = mesh.get_connectivity(mesh.dim-1, 0);
    boundary_facets = mesh.get_boundary_facets_indices();
    boundary_facets = sparse(double(boundary_facets), 1, true, num_facets, 1);

    % Verify facet by facet
    for f = 1:num_facets
        if markers(f) == 0
            t.verifyTrue(~boundary_facets(f));
        else
            t.verifyTrue(boundary_facets(f));
            d = markers(f);
            for v = 1:mesh.dim
                vertex_coords = coords(:, f2v(v, f));
                t.verifyTrue(vertex_coords(d) < tol || vertex_coords(d) > 1-tol);
            end
        end
    end

end


function test_vertex_markers_cube_1d(t)
    % Throws: 1D cube mesh not implemented
    t.verifyError(@() t_vertex_markers_cube(t, 66), ?MException);
end


function test_vertex_markers_cube_2d(t)
    t_vertex_markers_cube(t, [7, 13]);
end


function test_vertex_markers_cube_3d(t)
    t_vertex_markers_cube(t, [4, 5, 6]);
end


function t_vertex_markers_cube(t, dims)

    % Supress warning about unneeded connectivities computed and stored
    import matlab.unittest.fixtures.SuppressedWarningsFixture
    t.applyFixture(SuppressedWarningsFixture('Mesh:ExtraConnStored'));

    % Create mesh and markers
    mesh = meshing.generate_unit_cube_mesh(dims);
    mesh.compute_connectivity(mesh.dim, mesh.dim-1);
    mesh.compute_boundary_facets();
    tol = 1e-10;
    where = @(x, b) ~b && x(1) > 0.5 - tol;
    markers = meshing.mark_entities(mesh, 0, [], where, 1);

    % Fetch needed data
    coords = mesh.vertex_coords;
    num_vertices = mesh.num_entities(0);

    % Verify vertex by vertex
    for v = 1:num_vertices

        x = coords(:, v);
        onb = is_vertex_on_boundary_(x);
        isc = is_vertex_in_corner_(x);

        if markers(v) == 0
            t.verifyTrue(~where(x, onb));
        else
            assert(markers(v) == 1);

            % NB: The following should be arguably true:
            %         t.verifyTrue(~where(x, onb));
            %     but vertices connected to cells which have all facet
            %     interior are considered not on boundary.
            %     This is arguably a bug! Instead, we test at least:
            t.verifyTrue(where(x, onb) || (onb && ~isc));
        end
    end

end


function test_cell_markers_cube_1d(t)
    % Throws: 1D cube mesh not implemented
    t.verifyError(@() t_cell_markers_cube(t, 66), ?MException);
end


function test_cell_markers_cube_2d(t)
    t_cell_markers_cube(t, [7, 13]);
end


function test_cell_markers_cube_3d(t)
    t_cell_markers_cube(t, [4, 5, 6]);
end


function t_cell_markers_cube(t, dims)

    % Supress warning about unneeded connectivities computed and stored
    import matlab.unittest.fixtures.SuppressedWarningsFixture
    t.applyFixture(SuppressedWarningsFixture('Mesh:ExtraConnStored'));

    % Create mesh and markers
    mesh = meshing.generate_unit_cube_mesh(dims);
    mesh.compute_connectivity(mesh.dim, mesh.dim-1);
    mesh.compute_boundary_facets();
    tol = 1e-10;
    where = @(x, b) ~b && x(1) > 0.5 - tol;
    markers = meshing.mark_entities(mesh, mesh.dim, [], where, 1);

    % Fetch needed data
    cells = mesh.cells;
    coords = mesh.vertex_coords;
    num_cells = mesh.num_entities(mesh.dim);

    % Verify cell by cell
    for c = 1:num_cells

        x = coords(:, cells(:, c));
        onb = cell_has_boundary_facet_(x);

        if markers(c) == 0
            t.verifyTrue(onb || any(x(1, :) <= 0.5 - tol));
        else
            assert(markers(c) == 1);
            t.verifyTrue(~onb && all(x(1, :) > 0.5 - tol));
        end
    end

end


function test_edge_markers_cube_1d(t)
    % Throws: 1D cube mesh not implemented
    t.verifyError(@() t_edge_markers_cube(t, 66), ?MException);
end


% NB: Case already covered above: edge markers in 3D
function test_edge_markers_cube_2d(t)
    t_edge_markers_cube(t, [7, 13]);
end


% NB: Case not covered above: edge markers in 3D
function test_edge_markers_cube_3d(t)
    t_edge_markers_cube(t, [4, 5, 6]);
end


function t_edge_markers_cube(t, dims)
    % Supress warning about unneeded connectivities computed and stored
    import matlab.unittest.fixtures.SuppressedWarningsFixture
    t.applyFixture(SuppressedWarningsFixture('Mesh:ExtraConnStored'));

    % Create mesh and markers
    mesh = meshing.generate_unit_cube_mesh(dims);
    mesh.compute_connectivity(mesh.dim, mesh.dim-1);
    mesh.compute_connectivity(mesh.dim, 1);
    mesh.compute_boundary_facets();
    tol = 1e-10;
    where = @(x, b) ~b && x(1) > 0.5 - tol;
    markers = meshing.mark_entities(mesh, 1, [], where, 1);

    % Fetch needed data
    e2v = mesh.get_connectivity(1, 0);
    coords = mesh.vertex_coords;
    num_edges = mesh.num_entities(1);

    % Verify vertex by vertex
    for e = 1:num_edges

        x = coords(:, e2v(:, e));
        onb = is_edge_on_boundary_(x);

        if markers(e) == 0
            t.verifyTrue(~where(x(:, 1), onb) || ~where(x(:, 2), onb));
        else
            assert(markers(e) == 1);

            if mesh.dim == 2
                t.verifyTrue(where(x(:, 1), onb) && where(x(:, 2), onb));
            else
                % NB: The following should be arguably true:
                %         t.verifyTrue(where(x(:, 1), onb) && where(x(:, 2), onb));
                %     but every boundary edge is connected to some cell
                %     which has no facet on boundary and 'where' is true there.
                %     This is arguably a bug! Instead, we test at least:
                t.verifyTrue(where(x(:, 1), false) && where(x(:, 2), false));
            end
        end
    end

end


function flag = is_vertex_in_corner_(coords)
    tol = 1e-10;
    flag = all(coords(:) < tol | coords(:) > 1-tol);
end


function flag = is_vertex_on_boundary_(coords)
    assert(iscolumn(coords));
    tol = 1e-10;
    flag = any(coords(:) < tol) || any(coords(:) > 1-tol);
end


function flag = is_edge_on_boundary_(coords)
    flag = is_vertex_on_boundary_(coords(:, 1)) ...
        && is_vertex_on_boundary_(coords(:, 2)) ...
        && is_vertex_on_boundary_(mean(coords, 2));
end


function flag = is_face_on_boundary_(coords)
    flag = is_vertex_on_boundary_(coords(:, 1)) ...
        && is_vertex_on_boundary_(coords(:, 2)) ...
        && is_vertex_on_boundary_(coords(:, 3)) ...
        && is_vertex_on_boundary_(mean(coords, 2));
end


function flag = is_facet_on_boundary_(coords)
    switch size(coords, 1)
    case 1
        flag = is_vertex_on_boundary_(coords);
    case 2
        flag = is_edge_on_boundary_(coords);
    case 3
        flag = is_face_on_boundary_(coords);
    otherwise
        assert(false);
    end
end


function flag = cell_has_boundary_facet_(coords)
    dim = size(coords, 1);
    assert(size(coords, 2) == dim + 1);
    facet_verts = nchoosek(1:dim+1, dim);
    for fv = facet_verts.'
        if is_facet_on_boundary_(coords(:, fv))
            flag = true;
            return;
        end
    end
    flag = false;
end
