% Bulk impedance
sigma = @(x, c) 1 + x(1)^2 + x(2)^2;

% Electrodes configuration
impedances = [1/2, 1/3, 1/3, 1/4, 1/4, 1/3, 1/3, 1/2];
source_currents = [10, 30, 30, 60, -60, -30, -30, -10];
num_electrodes = numel(impedances);
phi_electrodes = linspace(0, pi, num_electrodes+2);
phi_electrodes = phi_electrodes(2:end-1);
x_electrodes = [cos(phi_electrodes); sin(phi_electrodes)];

% Sanity check
assert(numel(source_currents) == num_electrodes);
assert(sum(source_currents) == 0);
assert(all(real(impedances) > 0));
assert(all(isreal(impedances)), 'TODO: implement/test complex scenario');

% Build mesh and needed connectivities
n = 32;
mesh = meshing.generate_half_disc_mesh(n);
mesh.compute_connectivity(mesh.dim, mesh.dim-1);
mesh.compute_connectivity(mesh.dim-1, mesh.dim);
mesh.compute_boundary_facets();

% Build boundary facet markers
tol = pi/(3*n);
markers = sparse(mesh.num_entities(mesh.dim-1), 1);
for j = 1:num_electrodes
    boundary = @(x, b) b && vecnorm(x - x_electrodes(:, j)) <= tol;
    markers = meshing.mark_entities(mesh, mesh.dim-1, markers, boundary, j);
end

% Build Lagrange function space of order 1
element = fe.create_lagrange_element(mesh.dim, 1);
dofmap = assembling.build_dofmap(mesh, element);

% Assemble Laplacian matrix
sigma_p0 = interpolate_to_p0(mesh, sigma);
assert(all(real(sigma_p0) > 0));
assert(all(isreal(sigma_p0)), 'TODO: implement/test complex scenario');
A11 = assembling.assemble_laplace(dofmap, sigma_p0, 0);

% Build vector of impedances (indexed by cells adjacent to respective boundary facet)
f2c = mesh.get_connectivity(mesh.dim-1, mesh.dim);
num_cells = mesh.num_entities(mesh.dim);
z = spalloc(num_cells, 1, nnz(markers));
for j = 1:num_electrodes
    cells_electrodes = f2c(:, markers == j);

    % boundary facets have just one adjacent cell
    assert(all(cells_electrodes(2, :) == 0));
    cells_electrodes = cells_electrodes(1, :);

    z(cells_electrodes) = impedances(j);  %#ok<SPRIX>
end

% Pop some plots
figure();
subplot(2, 1, 1); meshing.plot_mesh(mesh, 'facet_markers', markers);
subplot(2, 1, 2); meshing.plot_mesh(mesh, 'cell_markers', full(z));

% Add electrode impedance term to 11-block
g = @(x, n, c) 1/z(c);
degree_rise = 0;
A11 = assembling.assemble_robin_matrix(A11, dofmap, markers > 0, true, ...
                                       g, degree_rise);

% Assemble 12-block
g = @(x, n, c) -1/z(c);  % NB the minus sign!
degree_rise = 0;
A12 = [];
for j = 1:num_electrodes
    b = assembling.assemble_neumann_source(dofmap, markers, j, ...
                                           g, degree_rise);
    A12 = [A12, sparse(b(:))];  %#ok<AGROW>
    clear('b');
end

% Assemble 22-block
electrode_areas = zeros(num_electrodes, 1);
f2v = mesh.get_connectivity(mesh.dim-1, 0);
for j = 1:num_electrodes
    facets = find(markers == j);
    for f = facets.'
        verts = mesh.vertex_coords(:, f2v(:, f));
        electrode_areas(j) = electrode_areas(j) + facet_area(verts);
    end
end
A22_diag = electrode_areas(:) ./ impedances(:);
A22 = spdiags(A22_diag, 0, num_electrodes, num_electrodes);

% Assemble blocks together
A = [A11, A12; A12', A22];

% Plot sparsity pattern and report about constant nullspace
figure(); spy(A);
fprintf('Before BC: ||A*ones(n, 1)|| = %g\n', norm(sum(A, 2)));

% Asseble rhs
b = zeros(size(A, 1), 1);
b(dofmap.dim+1:end) = source_currents;

% Add additional BC to remove constant nullspace
dofs = 1; vals = 0;
[A, b] = assembling.apply_dirichlet_bc(A, b, dofs, vals);
fprintf('After BC:  ||A*ones(n, 1)|| = %g\n', norm(sum(A, 2)));

% Solve
x = A\b;
x1 = x(1:dofmap.dim);
x2 = x(dofmap.dim+1:end);

% Report solution at electrodes
fprintf('Electrode %d: potential = %f\n', [1:num_electrodes; x2.']);
fprintf('Sum of electrode potentials = %g\n', sum(x2));
assert(sum(x2) < 1e-10);

% Interpolate to vertices and plot
vertex_values = assembling.interpolate_vertex_values(dofmap, x);
figure();
subplot(1, 2, 1); meshing.plot_vertex_values(mesh, vertex_values);
subplot(1, 2, 2); meshing.plot_vertex_values(mesh, vertex_values);
hold('on'); plot3(x_electrodes(1, :), x_electrodes(2, :), x2, 'ro'); hold('off');


function V = facet_area(vertices)
    [dim, num_vertices] = size(vertices);
    assert(dim == num_vertices);
    switch dim
    case 2
        V = norm(vertices(:, 1) - vertices(:, 2));
    case 3
        % Heron's formula
        a = norm(vertices(:, 2) - vertices(:, 1));
        b = norm(vertices(:, 3) - vertices(:, 1));
        c = norm(vertices(:, 3) - vertices(:, 2));
        s = 0.5*(a+b+c);
        V = sqrt(s*(s-a)*(s-b)*(s-c));
    end
end


function u = interpolate_to_p0(mesh, func)
    element = fe.create_p0_element(mesh.dim);
    dofmap = assembling.build_dofmap(mesh, element);
    u = assembling.interpolate(dofmap, func);
end
