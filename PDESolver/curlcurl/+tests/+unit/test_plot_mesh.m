function tests = test_plot_mesh
    % Run all function tests in this file
    tests = functiontests(localfunctions);
end


function test_plot_mesh_doc_snippet_2d(t)
    figure(5577);

    % Supress warning about unneeded connectivities computed and stored
    import matlab.unittest.fixtures.SuppressedWarningsFixture
    t.applyFixture(SuppressedWarningsFixture('Mesh:ExtraConnStored'));

    % Test that the snippet from doc of 'plot_mesh' runs
    mesh = meshing.generate_unit_cube_mesh([3, 3]);
    mesh.compute_connectivity(1, 0);
    cm = mod(1:mesh.num_entities(mesh.dim), 3);
    em = mod(1:mesh.num_entities(1), 4);
    vm = mod(1:mesh.num_entities(0), 5);
    meshing.plot_mesh(mesh, 'cell_markers', cm, 'edge_markers', em, ...
                      'vertex_markers', vm);

    close(5577);
end


function test_plot_mesh_doc_snippet_3d_1(t)
    figure(5678);

    % Supress warning about unneeded connectivities computed and stored
    import matlab.unittest.fixtures.SuppressedWarningsFixture
    t.applyFixture(SuppressedWarningsFixture('Mesh:ExtraConnStored'));

    % Test that the snippet from doc of 'plot_mesh' runs
    mesh = meshing.generate_unit_cube_mesh([3, 3, 3]);
    mesh.compute_connectivity(1, 0);
    cm = mod(1:mesh.num_entities(mesh.dim), 3);
    em = mod(1:mesh.num_entities(1), 4);
    vm = mod(1:mesh.num_entities(0), 5);
    pl = @(x) 3*x(1,:) + 1*x(2,:) + 2*x(3,:) - 2;
    meshing.plot_mesh(mesh, 'cell_markers', cm, 'edge_markers', em, ...
                      'vertex_markers', vm, 'plane', pl);

    close(5678);
end


function test_plot_mesh_doc_snippet_3d_2(t)
    figure(5679);

    % Supress warning about unneeded connectivities computed and stored
    import matlab.unittest.fixtures.SuppressedWarningsFixture
    t.applyFixture(SuppressedWarningsFixture('Mesh:ExtraConnStored'));

    % Test that the snippet from doc of 'plot_mesh' runs
    mesh = meshing.generate_unit_cube_mesh([3, 3, 3]);
    pl = @(x) x(2,:) - 1/3;
    meshing.plot_mesh(mesh, 'plane', pl, 'preferred_side', +1);

    close(5679);
end


function test_plot_mesh_3d_unstructured_intersection(~)
    figure(5680);

    domain_r = 20;
    pt = -domain_r/2:2:domain_r/2;
    pt = [pt; 0*pt; 0*pt].';
    pl = @(x) x(2, :);
    [X, Z] = ndgrid([-domain_r, domain_r], [-domain_r, 5]);

    [mesh, pm] = meshing.generate_mesh3D('point', pt, ...
                                         'marker', 0, ...
                                         'domain_r', domain_r);

    set(gcf, 'Units', 'normalized', 'OuterPosition', [.25, 0, .5, 1]);
    subplot(3, 1, 1);
        draw_subplot(1);
    subplot(3, 1, 2);
        draw_subplot(0);
    subplot(3, 1, 3);
        draw_subplot(-1);

    close(5680);

    function draw_subplot(side)

        meshing.plot_mesh(mesh, 'vertex_marker', pm, ...
                                'plane', pl, ...
                                'preferred_side', side);
        hold on
            surf(X, 0*X, Z, 'EdgeColor', 'Red', 'FaceColor', 'Red', ...
                 'DisplayName', 'cutting plane');
            plot3(0, -2*domain_r/3, -domain_r/2, '.r', ...
                  'MarkerSize', 24, ...
                  'DisplayName', 'view point');
        hold off
        view([0, -.25, 1]);
        title(sprintf('side: %d', side));
    end
end
