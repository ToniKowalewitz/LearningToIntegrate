function check_convergence(elements, assemble_mat, assemble_rhs, assemble_err, mesh_resolutions, u, boundary)
% CHECK_CONVERGENCE  Solve on sequence of meshes and compare with exact solution
%
% SYNTAX
%
%   check_convergence(elements, assemble_func, mesh_resolutions, u)
%
% INPUT PARAMETERS
%
%   elements        ... array with finite elements
%   assemble_mat    ... function taking dofmap and returning system matrix
%   assemble_rhs    ... function taking dofmap and returning system right-hand side
%   assemble_err    ... function taking dofmap and x and returning L2 and energy errors
%   mesh_resolution ... cells with array of integers (indexed by FE order)
%   u               ... exact solution
%   boundary        ... function taking coordinate and evaluating true on boundary

    % Loop over FE elements
    for element = elements

        dim = element.simplex.dim;
        order = element.order;

        % Variables for storing results
        num_dofs{order} = [];  %#ok<AGROW>
        errs_L2{order} = [];  %#ok<AGROW>
        errs_en{order} = [];  %#ok<AGROW>

        % Loop over prescribed mesh resolutions
        for n = mesh_resolutions{order}
            fprintf('*** %s order %d ***\n', element.family, order);
            fprintf('*** Unit cube mesh %d^%d ***\n', n, dim);

            % Build mesh
            ns(1:dim) = n;
            tic; mesh = meshing.generate_unit_cube_mesh(ns);
            fprintf('Generate unit cube mesh in time %.3f sec\n', toc);

            % Build additional mesh connectivities (for building dofmap and BCs)
            tic;
            w = warning('off', 'Mesh:ExtraConnStored');  % Mute warning
            for d = element.get_dof_entity_dims()
                mesh.compute_connectivity(mesh.dim, d);
                mesh.clear_connectivity(d, 0);
            end
            mesh.compute_connectivity(mesh.dim, mesh.dim-1);
            mesh.clear_connectivity(mesh.dim-1, 0);
            warning(w);  % Restore warning
            mesh.compute_boundary_facets();
            fprintf('Build mesh connectivity in time %.3f sec\n', toc);

            % Prepare geometric queries
            % FIXME: Do only if necessary?
            mesh.init_geometric_queries();

            % Build dofmap
            tic; dofmap = assembling.build_dofmap(mesh, element);
            fprintf('Build dofmap in time %.3f sec\n', toc);

            % Clear connectivities not needed any more
            for d = 1:mesh.dim-1
                mesh.clear_connectivity(mesh.dim, d);
            end

            % Assemble operator
            tic; A = assemble_mat(dofmap);
            fprintf('Assemble operator in time %.3f sec\n', toc);

            % Assemble rhs
            tic; b = assemble_rhs(dofmap);
            fprintf('Assemble right-hand side in time %.3f sec\n', toc);

            % Compute boundary dof indices and dof values
            tic; [bdofs, bvals] = assembling.build_dirichlet_dofs(dofmap, {boundary}, {u});
            fprintf('Build Dirichlet map in time %.3f sec\n', toc);

            % Apply boundary conditions to system matrix and rhs
            tic; [A, b] = assembling.apply_dirichlet_bc(A, b, bdofs, bvals);
            fprintf('Apply Dirichlet bc in time %.3f sec\n', toc);

            % Solve linear system
            tic; x = A\b;
            fprintf('Solve in time %.3f sec\n', toc);

            % Assemble errors and store them
            tic; [err_L2, err_en] = assemble_err(dofmap, x);
            fprintf('Assemble errors in %.3f sec\n', toc);
            fprintf('L2 error: %f, en error: %f\n', err_L2, err_en);
            errs_L2{order}(end+1) = err_L2;
            errs_en{order}(end+1) = err_en;
            num_dofs{order}(end+1) = dofmap.dim;

            % Interpolate exact solution on space (just for fun; coverage)
            tic; x_ex = assembling.interpolate(dofmap, @(x, c) u(x));  %#ok<NASGU>
            fprintf('Interpolate exact solution in time %.3f sec\n', toc);

            % Interpolate FE solution at vertices and plot in 2d
            if strcmp(element.family, 'Lagrange')
                tic; vertex_values = assembling.interpolate_vertex_values(dofmap, x);
                fprintf('Interpolate vertex values in time %.3f sec\n', toc);
                if dim == 2
                    figure(999);
                    meshing.plot_vertex_values(mesh, vertex_values);
                end
            end

            fprintf('\n');
        end
        fprintf('\n');
    end

    % Plot convergence rates
    figure;
    subplot(2, 1, 1); plot_rates('L^2', num_dofs, errs_L2);
    subplot(2, 1, 2); plot_rates('energy', num_dofs, errs_en);

    % Report convergence rates, LS fit and last rate
    report_rates(dim, 'L^2', num_dofs, errs_L2);
    report_rates(dim, 'energy', num_dofs, errs_en);
end


function plot_rates(errtype, num_dofs, errs)
    for data = [num_dofs; errs]
        loglog(data{1}, data{2}, '-o'); hold on
    end
    xlabel('Number dofs'), ylabel(strcat(errtype, ' error'))
end


function report_rates(dim, errtype, num_dofs, errs)
    fprintf('%s-error vs (#dofs)^(1/%d): LS fit, last rate\n', errtype, dim);
    for order = 1:size(num_dofs, 2)
      fit = polyfit(log10(num_dofs{order}), dim*log10(errs{order}), 1);
      last_rate = dim * log(errs{order}(end)/errs{order}(end-1)) / log(num_dofs{order}(end)/num_dofs{order}(end-1));
      fprintf('    order %d, LS %1.2f, last %1.2f\n', order, fit(1), last_rate);
    end
end
