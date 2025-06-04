% Build mesh and needed connectivities
n = 32;
mesh = meshing.generate_half_disc_mesh(n);
mesh.compute_connectivity(mesh.dim, mesh.dim-1);

% Build boundary facet markers
mesh.compute_boundary_facets();
marker_value = 616;
tol = 1e-10;
boundary = @(x, b) b && norm(x) > 1-tol;
markers = meshing.mark_entities(mesh, mesh.dim-1, [], boundary, marker_value);
mesh.clear_connectivity(mesh.dim-1, 0);

% Build Lagrange function space
element = fe.create_lagrange_element(2, 2);
dofmap = assembling.build_dofmap(mesh, element);

% Assemble Dirac source
mesh.init_geometric_queries();
x0 = [0; 0];
b = assembling.assemble_point_sources(dofmap, x0);
mesh.clear_geometric_queries();

% Assemble Laplacian matrix
A = assembling.assemble_laplace(dofmap, 1, 0);

% Add Robin term to the matrix
g = @(x, n, c) 1/pi;
degree_rise = 0;
A = assembling.assemble_robin_matrix(A, dofmap, markers, marker_value, ...
                                     g, degree_rise);

% Clear unneeded connectivity
mesh.clear_connectivity(mesh.dim, mesh.dim-1);
mesh.clear_boundary_facets();

% Solve
x = A\b;

% Interpolate to vertices and plot
vertex_values = assembling.interpolate_vertex_values(dofmap, x);
meshing.plot_vertex_values(mesh, vertex_values);

% Compare with manufactured solution
u = @(x) 1-1/pi*log(norm(x));
degree_rise = 1;
[err_L2, err_H1] = assembling.assemble_error_h1(dofmap, x, u, degree_rise);
fprintf('Err L2: %f, Err H1: %f \n', err_L2, err_H1);
