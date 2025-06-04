function [err_L2, err_H1] = assemble_error_h1(dofmap, x_h, u_ex, degree_rise)
    % Assemble ||u_h-u_x||_2 and ||grad(u_h-u_x)||_2
    %
    % Assemble ||u_h-u_x||_2 and ||grad(u_h-u_x)||_2
    % where u_h is given by FE and exansion coefficients x_h.
    % Parameter degree_rise, non-negative integer, specifies
    % how much should be quadrature degree increased
    % to account for non-polynomiality of u_ex. If function
    % u_ex lived in the FE space, then zero degree_rise
    % would give exact quadrature.
    %
    % Exact solution u_ex has to admit symbolic computation
    % of its grad.


    assert(strcmp(dofmap.element.mapping, 'affine'));

    dim = dofmap.mesh.dim;

    % Pick quadrature rule
    quad_degree = 2*dofmap.element.order + degree_rise;
    [x, w] = fe.make_quadrature_rule(dim, quad_degree);
    num_quad_points = size(x, 1);

    % Tabulate basis and basis grads at quadrature points
    basis_grads = zeros(dofmap.element.fe_space_dim, dim, num_quad_points);
    basis       = zeros(dofmap.element.fe_space_dim,      num_quad_points);
    for k = 1:num_quad_points
        basis_grads(:, :, k) = dofmap.element.tabulate_basis_grad(x(k, :));
        basis      (:,    k) = dofmap.element.tabulate_basis     (x(k, :));
    end

    % Compute symbolic grad of exact solution
    p = sym('p', [dim, 1], 'real');
    u_ex_grad_sym = gradient(u_ex(p), p).';
    u_ex_grad = matlabFunction(u_ex_grad_sym, 'Vars', {p});

    % Fetch some data and preallocate temporaries
    cells = dofmap.mesh.cells;
    coords = dofmap.mesh.vertex_coords;
    cell_dofs = dofmap.cell_dofs;
    num_cells = size(dofmap.mesh.cells, 2);
    origin = zeros(dim, 1);
    jac = zeros(dim, dim);
    jac_inv = zeros(dim, dim);
    detJ = zeros(1, 1, 'double');
    detJinv = zeros(1, 1, 'double');
    x_h_dofs = zeros(1, dofmap.element.fe_space_dim);
    X = zeros(dim, 1);
    temp = zeros(dim, 1);

    % Loop over cells
    err_L2 = 0;
    err_H1 = 0;
    for c = 1:num_cells

        % Compute geometric quantities
        origin(:) = coords(:, cells(dim+1, c));
        jac(:, :) = coords(:, cells(1:dim, c)) - coords(:, cells(dim+1, c));
        jac_inv(:, :) = inv(jac);
        detJ(1) = abs(det(jac));
        detJinv(1) = 1/det(jac);

        % Fetch relevant dofs from x_h
        x_h_dofs(:) = x_h(cell_dofs(:, c));

        % Loop over quadrature points and assembly cell integrals
        for k = 1:num_quad_points
            X(:) = jac*x(k, :).' + origin;
            temp(1) = x_h_dofs*basis(:, k) - u_ex(X);
            err_L2(1) = err_L2 + detJ*w(k)*temp(1)'*temp(1);
            temp(:) = x_h_dofs*basis_grads(:, :, k)*jac_inv - u_ex_grad(X);
            err_H1(1) = err_H1 + detJ*w(k)*(temp'*temp);
        end

    end

    % Take square roots
    assert(abs(imag(err_L2)) < eps*real(err_L2), 'Expected real error');
    assert(abs(imag(err_H1)) < eps*real(err_H1), 'Expected real error');
    err_L2 = sqrt(abs(err_L2));
    err_H1 = sqrt(abs(err_H1));

end
