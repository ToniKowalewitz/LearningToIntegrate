function verify_mesh_topology(t, mesh)
    % Check has sane topology data and is connected

    % Supress warning about unneeded connectivities computed and stored
    import matlab.unittest.fixtures.SuppressedWarningsFixture
    t.applyFixture(SuppressedWarningsFixture('Mesh:ExtraConnStored'));

    % Check that each cell has distinct vertices
    t.verifyTrue(all(all(unique(mesh.cells, 'rows') == mesh.cells)));

    % Get cell-facet connectivity
    mesh.compute_connectivity(mesh.dim, mesh.dim-1);
    mesh.clear_connectivity(mesh.dim-1, 0);
    c2f = mesh.get_connectivity(mesh.dim, mesh.dim-1);

    % Check that mesh is connected, i.e., each cell contains
    % at least one facet which is connected to two cells
    facet_counts = reshape(grouptransform(c2f(:), c2f(:), @numel), size(c2f));
    t.verifyTrue(all(max(facet_counts) == 2));

    % Sanity check
    t.verifyEqual(unique(facet_counts), [1; 2]);
end
