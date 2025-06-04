function dofmap = build_dofmap(mesh, element)
  % Build dofmap for element and mesh.
  %
  % FIXME: Fix the docstring

  assert(mesh.dim == element.simplex.dim, ...
         'Nonmatching mesh and element dimensions');

  % Build cell dofs - mapping of cells to global dofs.
  % cell_dofs(i, j) is global dof number for cell j and local dof
  % index i.
  cell_dofs = zeros(element.fe_space_dim, size(mesh.cells, 2), 'uint32');
  offset = 0;

  for dim = 0:mesh.dim
    local_dofs = element.get_entity_dofs(dim);
    dofs_per_entity = size(local_dofs, 1);

    if dofs_per_entity == 0
        continue
    end

    connectivity = mesh.get_connectivity(mesh.dim, dim);

    for k = 1:dofs_per_entity
      k1 = k - dofs_per_entity + offset;
      cell_dofs(local_dofs(k, :), :) = dofs_per_entity*connectivity + k1;
    end

    offset = offset + dofs_per_entity*mesh.num_entities(dim);
  end

  % Build dofmap structure and store what it depends on
  % FIXME: I am not sure about copy/reference semantics here!!!
  dofmap = struct;
  dofmap.cell_dofs = cell_dofs;
  dofmap.dim = offset;
  dofmap.mesh = mesh;
  dofmap.element = element;
end
