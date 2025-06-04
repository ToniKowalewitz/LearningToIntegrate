function [L, S, M] = assemble_laplace(dofmap, sigma, lambda)
    % Assemble diffusion-reaction operator
    %
    %   S_{ij} = \int \sigma \nabla \phi_i \cdot \nabla \phi_j
    %   M_{ij} = \int \lambda \phi_i \phi_j
    %   L = S - M
    %
    %   \phi_i          - H^1-conforming basis functions
    %   \sigma, \lambda - scalar or P_0 coefficients indexed by cell index
    %
    % INPUT PARAMETER
    %   dofmap ... Struct, containing mesh and FE element objects
    %              (the latter contain H^1-conforming basis functions)
    %              as well as the cell-2-DOF mapping.
    %   sigma  ... Vector of P_0 element values, indexed by cell index.
    %              Alternatively can be given as scalar (applying to all).
    %   lambda ... See sigma.
    %
    % OUTPUT PARAMETER
    %   L ... Sparse matrix [dofmap.dim x dofmap.dim]
    %   S ... Sparse matrix [dofmap.dim x dofmap.dim]
    %   M ... Sparse matrix [dofmap.dim x dofmap.dim]

    assert(strcmp(dofmap.element.mapping, 'affine'));

    % Fetch data
    local_element_dim = dofmap.element.fe_space_dim;
    dim = dofmap.mesh.dim;
    num_cells = size(dofmap.mesh.cells, 2);
    cells = dofmap.mesh.cells;
    coords = dofmap.mesh.vertex_coords;
    cell_dofs = dofmap.cell_dofs;

    % Pick quadrature rule
    assemble_mass = ~( isscalar(lambda) && lambda == 0 );
    if assemble_mass
        quad_degree = 2*dofmap.element.order;
    else
        quad_degree = 2*(dofmap.element.order - 1);
    end
    [x, w] = fe.make_quadrature_rule(dim, quad_degree);
    num_quad_points = size(x, 1);

    % Wrap constant values as p0 function
    sigma_cell  = wrap_as_p0(sigma,  num_cells);
    lambda_cell = wrap_as_p0(lambda, num_cells);

    % Tabulate basis and basis gradients at quadrature points
    basis_grads = zeros(local_element_dim, dim, num_quad_points);
    basis       = zeros(local_element_dim,      num_quad_points);
    for k = 1:num_quad_points
        basis_grads(:, :, k) = dofmap.element.tabulate_basis_grad(x(k, :));
        basis      (:,    k) = dofmap.element.tabulate_basis     (x(k, :));
    end

    % Preallocate temporaries
    [L, S, M, temp_S, temp_M] = deal(zeros(local_element_dim, local_element_dim));
    jac = zeros(dim, dim);
    jac_inv = zeros(dim, dim);
    temp = zeros(local_element_dim, dim);
    detJ = zeros(1, 1, 'double');

    % Preallocate assembly data
    nnz = num_cells*local_element_dim^2;
    I = zeros(nnz, 1, 'double');
    J = zeros(nnz, 1, 'double');
    [V_L, V_S, V_M] = deal(zeros(nnz, 1, 'double'));
    offsets = uint32(1:local_element_dim^2);

    % Loop over cells
    for c = 1:num_cells

        % Compute geometric quantities
        jac(:, :) = coords(:, cells(1:dim, c)) - coords(:, cells(dim+1, c));
        detJ(1) = abs(det(jac));
        jac_inv(:, :) = inv(jac);

        % Zero cell matrices
        [L(:, :), S(:, :), M(:, :)] = deal(0);

        % Loop over quadrature points and assembly cell integrals
        for k = 1:num_quad_points
            temp(:, :) = basis_grads(:, :, k)*jac_inv;
            temp_S(:, :) = sigma_cell(c)*detJ*w(k)*(temp*temp.');
            temp_M(:, :) = lambda_cell(c)*detJ*w(k)*basis(:, k)*basis(:, k).';
            S(:, :) = S + temp_S;
            M(:, :) = M + temp_M;
            L(:, :) = L + temp_S;
            if assemble_mass
                L(:, :) = L - temp_M;
            end
        end

        % Compute global dof indices and store cell matrix
        [J(offsets), I(offsets)] = meshgrid(cell_dofs(:, c));
        V_L(offsets) = L;
        if nargout >= 2
            V_S(offsets) = S;
        end
        if nargout == 3
            V_M(offsets) = M;
        end
        offsets(:) = offsets + local_element_dim^2;

    end

    % Assemble sparse matrices
    L = sparse(I, J, V_L, dofmap.dim, dofmap.dim);
    if nargout >= 2
        S = sparse(I, J, V_S, dofmap.dim, dofmap.dim);
    end
    if nargout == 3
        M = sparse(I, J, V_M, dofmap.dim, dofmap.dim);
    end
end


function g = wrap_as_p0(f, num_cells)
    if isscalar(f)
        g = @(i) f;
    else
        assert(numel(f) == num_cells);
        g = f;
    end
end
