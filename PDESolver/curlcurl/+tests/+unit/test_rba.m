function tests = test_rba
    % Run all function tests in this file
    tests = functiontests(localfunctions);
end


function test_rba_scalar(t)

    function verify_rba_(ns, xs, expected_reduction, last_good_iteration)

        % First iteration is always good (compare to infinity)
        err0 = +inf;

        % Loop over number of poles
        for n = ns

            % Compute maximal difference between RBA and exp(-x) at gridpoints
            r = rba.create_rba_scalar(n);
            err = vecnorm(r(xs) - exp(-xs), inf);

            %fprintf('n=%d, err=%f, H=%f\n', n, err/err0, expected_reduction);

            if n <= last_good_iteration
                % Good iterations: error reduction in every iteration
                t.verifyLessThanOrEqual(err/err0, expected_reduction);
            else
                % Numerical noise: we reached double precision
                t.verifyGreaterThan(err/err0, expected_reduction);
            end

            err0 = err;
        end
    end

    [a, b] = deal(0, 10);
    step = (b-a)/1024;
    xs = a:step:b;
    ns = 2:18;

    % Theoretical reduction rate
    halphen_const = 9.289025491920818918755449435951;
    relaxation_factor = 1.02;
    H = relaxation_factor / halphen_const;

    % Test real case
    verify_rba_(ns, xs, H, 13);

    % Test complex case
    verify_rba_(ns, xs+1i, 0.25, inf);

end


function test_rba_matrix(t)

    % Supress warning about unneeded connectivities computed and stored
    import matlab.unittest.fixtures.SuppressedWarningsFixture
    t.applyFixture(SuppressedWarningsFixture('Mesh:ExtraConnStored'));

    % Build function space
    mesh = meshing.generate_unit_cube_mesh([32, 32]);
    element = fe.create_nedelec_element(mesh.dim, 1);
    for d = element.get_dof_entity_dims()
        mesh.compute_connectivity(mesh.dim, d);
        mesh.clear_connectivity(d, 0);
    end
    dofmap = assembling.build_dofmap(mesh, element);

    % Assemble system
    [M, A] = assembling.assemble_hcurl_operators(dofmap);
    b = ones(dofmap.dim, 1);

    % Apply BCs
    mesh.compute_boundary_facets();
    bc_dofs = assembling.build_dirichlet_dofs(dofmap, {@(x) true}, {@(x) [0, 0]});
    [A, b] = assembling.apply_dirichlet_bc(A, b, bc_dofs, 0);
    M = assembling.apply_dirichlet_bc(M, [], bc_dofs, 0);

    function x = solve_rba_(A, M, x0, t, num_poles)
        r = rba.create_rba_matrix(num_poles, A, M, x0, @(A, b) A\b);
        x = r(t);
    end

    x0 = M\b;
    t_end = 0.5;

    % Compute very precise solution (large number of poles)
    x_good = solve_rba_(A, M, x0, t_end, 15);

    % Theoretical reduction rate
    halphen_const = 9.289025491920818918755449435951;
    relaxation_factor = 1.6;
    H = relaxation_factor / halphen_const;

    % Check that error decreases by factor H adding a pole
    err0 = +inf;
    for num_poles = 2:10
        x = solve_rba_(A, M, x0, t_end, num_poles);
        err = norm(x_good - x, inf);

        t.verifyLessThanOrEqual(err/err0, H);

        err0 = err;
    end

    % Just run
    for num_poles = 11:18
        x = solve_rba_(A, M, x0, t_end, num_poles);  %#ok<NASGU>
    end

end
