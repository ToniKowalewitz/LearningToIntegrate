function coords = get_lagrange_dofs_coordinates(dofmap)
    % Provide global coordinates of all DOFs
    %
    % SYNTAX
    %   coords = get_lagrange_dofs_coordinates(dofmap)
    %
    % INPUT PARAMETER
    %   dofmap ... Struct, containing cell-to-dof mapping as well as the
    %              underlying mesh and FE element objects
    %
    % OUTPUT PARAMETER
    %   coords ... Matrix [dofmap.mesh.dim, dofmap.dim], of DOF
    %              coordinates

    assert(strcmp(dofmap.element.family, 'Lagrange'), ...
           'Expected Lagrange element');

    % Allocate result
    coords = zeros(dofmap.mesh.dim, dofmap.dim);

    % Interpolate coordinate-function (exploit nodality of element)
    for i = 1:dofmap.mesh.dim
        coords(i, :) = assembling.interpolate(dofmap, @(x, c) x(i));
    end
end
