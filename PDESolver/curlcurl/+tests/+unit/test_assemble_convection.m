function tests = test_assemble_convection
  % Run all function tests in this file
  tests = functiontests(localfunctions);
end


function test_assemble_convection_reference(t)
    mesh = meshing.generate_unit_cube_mesh([1, 1]);

    import matlab.unittest.fixtures.SuppressedWarningsFixture
    t.applyFixture(SuppressedWarningsFixture('Mesh:ExtraConnStored'));
    mesh.compute_connectivity(mesh.dim, mesh.dim-1);

    element = fe.create_lagrange_element(mesh.dim, 1);
    dofmap = assembling.build_dofmap(mesh, element);
    wind = @(x, c) [0, 1.0];
    degree_rise = 0;
    M = assembling.assemble_convection(dofmap, wind, degree_rise);
    M = full(M);

    % Reference value obtained by FEniCS code:
    %   from dolfin import *
    %   mesh = UnitSquareMesh(1, 1, 'left')
    %   V = FunctionSpace(mesh, 'P', 1)
    %   u, v = TrialFunction(V), TestFunction(V)
    %   a = inner(as_vector([0, 1]), grad(u)*v)*dx
    %   print(assemble(a).array())
    M_ref = [
        1/6, -1/6,   0,    0;
        1/6, -1/6, 1/6, -1/6;
        1/6, -1/6, 1/6, -1/6;
          0,    0, 1/6, -1/6;
    ];

    % Compare our value with reference
    perm = [4, 2, 3, 1];
    t.verifyEqual(M(perm, perm), M_ref, 'RelTol', 1e-14);
end
