function tests = test_hsl_mi20
  % Run all function tests in this file
  tests = functiontests(localfunctions);
end


function test_hsl_mi20_works(t)
    % Check that HSL_MI20 runs and gives solution

    % Supress warning
    import matlab.unittest.fixtures.SuppressedWarningsFixture
    t.applyFixture(SuppressedWarningsFixture('Mesh:ExtraConnStored'));

    % Assemble system
    mesh = meshing.generate_half_disc_mesh(16);
    mesh.compute_connectivity(mesh.dim, mesh.dim-1);
    element = fe.create_lagrange_element(mesh.dim, 2);
    dofmap = assembling.build_dofmap(mesh, element);
    A = assembling.assemble_laplace(dofmap, 1, -1);
    b = ones(size(A, 1), 1);

    % Solve by 10 V-cycles of HSL_MI20
    c = hsl_mi20_control();
    c.v_iterations = 10;
    amg = solving.HSLMI20(A, c);
    x = amg.solve(b);

    % Compate with exact solution
    t.verifyEqual(x, A\b, 'RelTol', 1e-7);
end
