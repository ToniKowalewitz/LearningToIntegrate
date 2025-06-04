function func = create_observation(electrode_coords, Mtx, Mrx, ...
                                   dofmap, facet_marker, bc_facet_tag, ...
                                   solver_type)
    % Prepare function for computing DC resistivity observations and Jacobian
    %
    % INPUT PARAMETER
    %   electrode_coords ... matrix [dim, num_electrodes]
    %   Mtx              ... matrix [num_electrodes, num_observations]
    %   Mrx              ... matrix [num_electrodes, num_observations]
    %   dofmap           ... struct, represents H^1-conforming function space
    %   facet_marker     ... (sparse) vector, marks facets on underlying mesh
    %   bc_facet_tag     ... scalar, facet marker value representing zero BC
    %   solver_type      ... optional argument, character vector
    %
    % OUTPUT PARAMETER
    %   func ... function of signature:
    %                d = compute_observation(sigma);
    %                [d, J] = compute_observation(sigma);

    if nargin < 7
        solver_type = 'backslash';
    end

    bc_dofs = build_dirichlet_dofs(dofmap, facet_marker, bc_facet_tag);
    F = assemble_electrodes(dofmap, bc_dofs, electrode_coords);
    assemble_forward_operator = create_fwd_assembler(dofmap, bc_dofs);
    switch solver_type
    case 'backslash'
        solver_fwd = create_fwd_solver_direct(F);
    case 'cg/amg'
        solver_fwd = create_fwd_solver_iterative(F);
    end
    [assemble_observation_only, assemble_observation_and_jacobian] = ...
        create_observation_assembler(dofmap, Mtx, Mrx);

    function [d, J] = compute_observation(sigma)
        % Compute observations and possibly Jacobian
        %
        % SYNOPSIS
        %   d = compute_observation(sigma);
        %   [d, J] = compute_observation(sigma);

        % TODO: Can we optimize for repeated assembly with known sparsity?
        % NB: Probably not worth with Matlab's sparse format
        A = assemble_forward_operator(sigma);

        % TODO: Solve into preallocated solution vector?
        U = solver_fwd(A);

        if nargout == 1
            d = assemble_observation_only(U, A);
        else
            [d, J] = assemble_observation_and_jacobian(U, sigma);
        end
    end

    func = @compute_observation;
end


function func = create_fwd_solver_direct(F)  %#ok<*MSNU,DEFNU>
    function U = solve(A)
        % Perform forward solve directly

        fprintf('Forward solve directly: ');
        t = tic();

        U = A\full(F);
        assert(~issparse(U));

        fprintf('%f seconds\n', toc(t));
    end

    func = @solve;
end


function func = create_fwd_solver_iterative(F)  %#ok<*MSNU,DEFNU>
    function U = solve(A)
        % Perform forward solve using CG/AMG

        fprintf('Forward solve CG/AMG...\n');
        t = tic();

        amg = solving.HSLMI20(A);
        U = zeros(size(F));
        for j = 1:size(F, 2)
            U(:, j) = pcg(A, F(:,j), 1e-8, 100, @amg.solve);
        end
        clear amg;

        assert(~issparse(U));

        fprintf('Forward solve CG/AMG: %f seconds\n', toc(t));
    end

    func = @solve;
end


function bc_dofs = build_dirichlet_dofs(dofmap, facet_markers, facet_tag)
    fprintf('Build Dirichlet dofs: ');
    t = tic();

    value = @(x) 0;

    w = warning('off', 'Mesh:ExtraConnStored');  % Mute warning
    dofmap.mesh.compute_connectivity(dofmap.mesh.dim, dofmap.mesh.dim-1);
    warning(w);  % Restore warning
    dofmap.mesh.clear_connectivity(dofmap.mesh.dim-1, 0);
    dofmap.mesh.compute_boundary_facets();

    [bc_dofs, ~] = assembling.build_dirichlet_dofs(dofmap, facet_markers, {facet_tag}, {value});

    fprintf('%f seconds\n', toc(t));
end


function F = assemble_electrodes(dofmap, bc_dofs, electrode_coords)
    fprintf('Assemble electrodes: ');
    t = tic();

    dofmap.mesh.init_geometric_queries();
    F = assembling.assemble_point_sources(dofmap, electrode_coords);
    dofmap.mesh.clear_geometric_queries();

    [~, F] = assembling.apply_dirichlet_bc([], F, bc_dofs, 0);

    assert(size(F, 2) == size(electrode_coords, 2));
    assert(issparse(F));

    fprintf('%f seconds\n', toc(t));
end


function func = create_fwd_assembler(dofmap, bc_dofs)
    function A = assemble_forward_operator(sigma)
        fprintf('Assembling forward operator: ');
        t = tic();

        A = assembling.assemble_laplace(dofmap, sigma, 0);
        [A, ~] = assembling.apply_dirichlet_bc(A, [], bc_dofs, 0);

        fprintf('%f seconds\n', toc(t));
    end

    func = @assemble_forward_operator;
end


function [f1, f2] = create_observation_assembler(dofmap, Mtx, Mrx)

    function d = assemble_observation_only(U, A)
        fprintf('Assembling observation only: ');
        t = tic();

        d = dot((A*U)*Mrx, U*Mtx, 1);
        d = d(:);

        assert(~issparse(d));

        fprintf('%f seconds\n', toc(t));
    end

    function [d, J] = assemble_observation_and_jacobian(U, sigma)
        fprintf('Assembling observation and sensitivity: ');
        t = tic();

        % Assemble sensitivity matrix
        J = assembling.assemble_laplace_sensitivity(dofmap, U, Mrx, Mtx);

        % Assemble observables
        d = J*sigma;

        % Scale by sigma and invert sign
        diag_sigma = spdiags(sigma, 0, numel(sigma), numel(sigma));
        J = -J*diag_sigma;

        assert(~issparse(d));
        assert(~issparse(J));

        fprintf('%f seconds\n', toc(t));
    end

    f1 = @assemble_observation_only;
    f2 = @assemble_observation_and_jacobian;
end
