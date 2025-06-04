function tests = test_orient_wire
  % Run all function tests in this file
  tests = functiontests(localfunctions);
end


function test_orient_wire_2d_wire(t)

    import matlab.unittest.fixtures.SuppressedWarningsFixture
    t.applyFixture(SuppressedWarningsFixture('Mesh:ExtraConnStored'));

    mesh = meshing.generate_half_disc_mesh(16);
    mesh.compute_connectivity(mesh.dim, mesh.dim-1);
    mesh.compute_connectivity(mesh.dim, 1);
    mesh.compute_boundary_facets();
    tol = 1e-10;

    boundary = @(x, b) b && x(2) < tol;
    markers = meshing.mark_entities(mesh, 1, [], boundary, 1);
    markers = meshing.orient_wire(mesh, markers, 1, 1, -1);

    e2v = mesh.get_connectivity(1, 0);
    coords = mesh.vertex_coords;

    e0 = find(markers, 1);
    x0 = coords(1, e2v(:, e0));
    orientation0 = markers(e0) * sign(x0(2)-x0(1));
    assert(orientation0 == 1 || orientation0 == -1);

    for e = find(markers).'
        x = coords(1, e2v(:, e));
        orientation = markers(e) * sign(x(2)-x(1));
        t.verifyEqual(orientation, orientation0);
    end

    %meshing.plot_mesh(mesh, 'edge_markers', markers);
end


function test_orient_wire_2d_loop(t)

    import matlab.unittest.fixtures.SuppressedWarningsFixture
    t.applyFixture(SuppressedWarningsFixture('Mesh:ExtraConnStored'));

    mesh = meshing.generate_half_disc_mesh(16);
    mesh.compute_connectivity(mesh.dim, mesh.dim-1);
    mesh.compute_connectivity(mesh.dim, 1);
    mesh.compute_boundary_facets();
    boundary = @(x, b) b;
    markers = meshing.mark_entities(mesh, 1, [], boundary, 1);
    markers = meshing.orient_wire(mesh, markers, 1, 1, -1);

    e2v = mesh.get_connectivity(1, 0);
    coords = mesh.vertex_coords;

    % Auxiliary boundary marker saying whether y == 0 or not
    tol = 1e-10;
    up = @(x, b) b && norm(x) > 1-tol;
    down = @(x, b) b && x(2) < tol;
    up_down = meshing.mark_entities(mesh, 1, [], up, -1);
    up_down = meshing.mark_entities(mesh, 1, up_down, down, 1);

    e0 = find(markers, 1);
    x0 = coords(1, e2v(:, e0));
    orientation0 = markers(e0) * sign(x0(2)-x0(1)) * up_down(e0);
    assert(orientation0 == 1 || orientation0 == -1);

    for e = find(markers).'
        x = coords(1, e2v(:, e));
        orientation = markers(e) * sign(x(2)-x(1)) * up_down(e);
        t.verifyEqual(orientation, orientation0);
    end

    %meshing.plot_mesh(mesh, 'edge_markers', markers);
end


function test_orient_wire_cube_1d(t)
    % Throws: 1D cube mesh not implemented
    t.verifyError(@() t_orient_wire_cube_wire(t, 66), ?MException);
    t.verifyError(@() t_orient_wire_cube_loop(t, 66), ?MException);
    t.verifyError(@() t_orient_wire_cube_empty(t, 66), ?MException);
    t.verifyError(@() t_orient_wire_cube_tree(t, 66), ?MException);
end


function test_orient_wire_cube_2d(t)
    t_orient_wire_cube_wire(t, [7, 13]);
    t_orient_wire_cube_loop(t, [7, 13]);
    t_orient_wire_cube_empty(t, [7, 13]);
    t_orient_wire_cube_tree(t, [7, 13]);
end


function test_orient_wire_cube_3d(t)
    t_orient_wire_cube_wire(t, [4, 5, 6]);
    t_orient_wire_cube_loop(t, [4, 5, 6]);
    t_orient_wire_cube_empty(t, [4, 5, 6])
    t_orient_wire_cube_tree(t, [2, 3, 1])
end


function t_orient_wire_cube_wire(t, dims)

    import matlab.unittest.fixtures.SuppressedWarningsFixture
    t.applyFixture(SuppressedWarningsFixture('Mesh:ExtraConnStored'));

    mesh = meshing.generate_unit_cube_mesh(dims);
    mesh.compute_connectivity(mesh.dim, mesh.dim-1);
    mesh.compute_connectivity(mesh.dim, 1);
    mesh.compute_boundary_facets();
    tol = 1e-10;

    boundary = @(x, b) b && all(x(2:end) < tol);
    markers = meshing.mark_entities(mesh, 1, [], boundary, 1);
    markers = meshing.orient_wire(mesh, markers, 1, 1, -1);

    e2v = mesh.get_connectivity(1, 0);
    coords = mesh.vertex_coords;

    e0 = find(markers, 1);
    x0 = coords(1, e2v(:, e0));
    orientation0 = markers(e0) * sign(x0(2)-x0(1));
    assert(orientation0 == 1 || orientation0 == -1);

    for e = find(markers).'
        x = coords(1, e2v(:, e));
        orientation = markers(e) * sign(x(2)-x(1));
        t.verifyEqual(orientation, orientation0);
    end

    %meshing.plot_mesh(mesh, 'edge_markers', markers);
end


function t_orient_wire_cube_loop(t, dims)

    import matlab.unittest.fixtures.SuppressedWarningsFixture
    t.applyFixture(SuppressedWarningsFixture('Mesh:ExtraConnStored'));

    mesh = meshing.generate_unit_cube_mesh(dims);
    mesh.compute_connectivity(mesh.dim, mesh.dim-1);
    mesh.compute_connectivity(mesh.dim, 1);
    mesh.compute_boundary_facets();
    tol = 1e-10;
    if mesh.dim == 2
        boundary = @(x, b) b;
        markers = meshing.mark_entities(mesh, 1, [], boundary, 1);
    elseif mesh.dim == 3
        boundary = @(x, b) b && x(3) > 1-tol && x(1)*(1-x(1)) < tol;
        markers = meshing.mark_entities(mesh, 1, [], boundary, 1);
        boundary = @(x, b) b && x(3) > 1-tol && x(2)*(1-x(2)) < tol;
        markers = meshing.mark_entities(mesh, 1, markers, boundary, 1);
    end
    markers = meshing.orient_wire(mesh, markers, 1, 1, -1);

    e2v = mesh.get_connectivity(1, 0);
    coords = mesh.vertex_coords;

    % Auxiliary boundary marker
    top = @(x, b) b && x(2) > 1-tol;
    left = @(x, b) b && x(1) < tol;
    bottom = @(x, b) b && x(2) < tol;
    right = @(x, b) b && x(1) > 1-tol;
    boundary = meshing.mark_entities(mesh, 1, [], top, -1);
    boundary = meshing.mark_entities(mesh, 1, boundary, left, -1);
    boundary = meshing.mark_entities(mesh, 1, boundary, bottom, 1);
    boundary = meshing.mark_entities(mesh, 1, boundary, right, 1);

    e0 = find(markers, 1);
    x0 = coords(:, e2v(:, e0));
    orientation0 = markers(e0) * sign(sum(x0(:, 2)-x0(:, 1))) * boundary(e0);
    assert(orientation0 == 1 || orientation0 == -1);

    for e = find(markers).'
        x = coords(:, e2v(:, e));
        orientation = markers(e) * sign(sum(x(:, 2)-x(:, 1))) * boundary(e);
        t.verifyEqual(orientation, orientation0);
    end

    %meshing.plot_mesh(mesh, 'edge_markers', markers);
end


function t_orient_wire_cube_empty(t, dims)

    import matlab.unittest.fixtures.SuppressedWarningsFixture
    t.applyFixture(SuppressedWarningsFixture('Mesh:ExtraConnStored'));

    mesh = meshing.generate_unit_cube_mesh(dims);
    mesh.compute_connectivity(mesh.dim, mesh.dim-1);
    mesh.compute_connectivity(mesh.dim, 1);
    mesh.compute_boundary_facets();
    boundary = @(x, b) false;
    markers = meshing.mark_entities(mesh, 1, [], boundary, 1);
    markers = meshing.orient_wire(mesh, markers, 1, 1, -1);

    t.verifyEqual(norm(markers), 0);

    %meshing.plot_mesh(mesh, 'edge_markers', markers);
end


function t_orient_wire_cube_tree(t, dims)

    import matlab.unittest.fixtures.SuppressedWarningsFixture
    t.applyFixture(SuppressedWarningsFixture('Mesh:ExtraConnStored'));

    mesh = meshing.generate_unit_cube_mesh(dims);
    mesh.compute_connectivity(mesh.dim, mesh.dim-1);
    mesh.compute_connectivity(mesh.dim, 1);
    mesh.compute_boundary_facets();

    % Empty edge markers
    markers = sparse([], [], [], mesh.num_entities(1), 1, 32);

    % Construct sorted edge vertices
    e2v = mesh.get_connectivity(1, 0);
    [verts, perm] = sort(e2v(:));

    % Loop over sorted edge vertices
    for j = 1:numel(verts)-2

        % Pick only vertices of degree >= 3
        if verts(j) == verts(j+1) && verts(j+1) == verts(j+2)

            % Edge indices connected to one vertex (of degree >= 3)
            edges = 1 + floor((perm(j:j+2)-1)/2);

            % Mark tree on such edges
            markers(:) = 0;      %#ok<SPRIX>
            markers(edges) = 1;  %#ok<SPRIX>

            % Check that error is raised
            fail = @() meshing.orient_wire(mesh, markers, 1, 1, -1);
            t.verifyError(fail, ?MException);

        end
    end

    %meshing.plot_mesh(mesh, 'edge_markers', markers);
end
