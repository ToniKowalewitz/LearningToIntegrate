function tests = test_rectangle_mesh
  % Run all function tests in this file
  tests = functiontests(localfunctions);
end


function test_docstring_code_snippet(t)  %#ok<INUSD>

    mesh = meshing.generate_rectangle_mesh(0, 2, 0, 4);
    figure();
    meshing.plot_mesh(mesh);

    mesh = meshing.generate_rectangle_mesh(0, 2, 0, 4, 'Hmax', 0.1);
    figure();
    meshing.plot_mesh(mesh);

    close();
    close();

end
