function tests = test_mumps
  % Run all function tests in this file
  tests = functiontests(localfunctions);
end


function test_mumps_simple(t)
    % Simple example of using MUMPS in matlab

    % Assemble system
    A = assemble_matrix_(t);
    b = ones(size(A, 1), 1);

    % Solve
    solver = solving.MUMPS(A, struct('SYM', 1));
    x = solver.solve(b);
    clear('solver');

    % Check solution
    err = norm(b-A*x, 'inf');
    t.verifyLessThanOrEqual(err, sqrt(eps));
end


function A = assemble_matrix_(t)
    mesh = meshing.generate_half_disc_mesh(16);
    element = fe.create_lagrange_element(mesh.dim, 2);
    import matlab.unittest.fixtures.SuppressedWarningsFixture
    t.applyFixture(SuppressedWarningsFixture('Mesh:ExtraConnStored'));
    for d = element.get_dof_entity_dims()
        mesh.compute_connectivity(mesh.dim, d);
        mesh.clear_connectivity(d, 0);
    end
    dofmap = assembling.build_dofmap(mesh, element);
    A = assembling.assemble_laplace(dofmap, 1, -1);
end
