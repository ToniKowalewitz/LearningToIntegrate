function func = create_observation(dofmap, g, electrode_coords)
    % Prepare function for computing MT impedance/tipper observations and Jacobian
    %
    % INPUT PARAMETERS
    %   dofmap           ... struct, represents H(curl)-conforming function space
    %   g                ... cell array [num_polarizations, num_frequencies]
    %                        of functions representing boundary values
    %   electrode_coords ... matrix [dim, num_electrodes]
    %
    % OUTPUT PARAMETER
    %   func ... function of signature:
    %                d = compute_observation(lambda1, lambda2, sigma);
    %                [d, J] = compute_observation(lambda1, lambda2, sigma);
    %
    %            The coefficients correspond to PDE:
    %
    %              curl curl Ej + ( lambda1(j) * sigma + lambda2(j) ) * Ej = 0
    %
    %            with lambda1(j), lambda2(j) complex constants, sigma
    %            P0-coefficient, and j goes over 1:num_frequencies.
    %
    %     INPUT PARAMETERS
    %
    %       lambda1 ... Vector [num_frequencies, 1] of complex numbers
    %       lambda2 ... Vector [num_frequencies, 1] of complex numbers
    %       sigma   ... Vector [num_cells,       1] of complex numbers
    %
    %     OUTPUT PARAMETERS
    %
    %       d ... Vector [num_observations, 1] of complex numbers;
    %             in 3D format is:
    %               [
    %                Zxx11,          Zxx21,          ..., Zxx<num_electrodes>1,
    %                Zyx11,          Zyx21,          ..., Zyx<num_electrodes>1,
    %                Tx11,           Tx21,           ..., Tx<num_electrodes>1,
    %                Zxy11,          Zxy21,          ..., Zxy<num_electrodes>1,
    %                Zyy11,          Zyy21,          ..., Zyy<num_electrodes>1,
    %                Ty11,           Ty21,           ..., Ty<num_electrodes>1,
    %
    %                Zxx12,          Zxx22,          ..., Zxx<num_electrodes>2,
    %                Zyx12,          Zyx22,          ..., Zyx<num_electrodes>2,
    %                Tx12,           Tx22,           ..., Tx<num_electrodes>2,
    %                Zxy12,          Zxy22,          ..., Zxy<num_electrodes>2,
    %                Zyy12,          Zyy22,          ..., Zyy<num_electrodes>2,
    %                Ty12,           Ty22,           ..., Ty<num_electrodes>2,
    %                .               .                    .
    %                .               .                    .
    %                .               .                    .
    %                Zxx1<num_freq>, Zxx2<num_freq>, ..., Zxx<num_electrodes><num_freq>,
    %                Zyx1<num_freq>, Zyx2<num_freq>, ..., Zyx<num_electrodes><num_freq>,
    %                Tx1<num_freq>,  Tx2<num_freq>,  ..., Tx<num_electrodes><num_freq>,
    %                Zxy1<num_freq>, Zxy2<num_freq>, ..., Zxy<num_electrodes><num_freq>,
    %                Zyy1<num_freq>, Zyy2<num_freq>, ..., Zyy<num_electrodes><num_freq>,
    %                Ty1<num_freq>,  Ty2<num_freq>,  ..., Ty<num_electrodes><num_freq>,
    %               ]
    %
    %             in 2D format is:
    %               [
    %                Zxy11,          Zxy21,          ..., Zxy<num_electrodes>1,
    %                Zxy12,          Zxy22,          ..., Zxy<num_electrodes>2,
    %                .               .                    .
    %                .               .                    .
    %                .               .                    .
    %                Zxy1<num_freq>, Zxy2<num_freq>, ..., Zxy<num_electrodes><num_freq>,
    %               ]
    %
    %       J ... (optional; call without J is much cheaper)
    %             derivative of d w.r.t. cell values of sigma, i.e.,
    %
    %               J(i, j) = partial d(i) / partial sigma(j)
    %
    %             where j goes over 1, 2, ..., num_cells.
    %
    % NOTES
    %
    %   It should be num_polarizations==2 in 3D and num_polarizations==1 in 2D,
    %   because the implemented observation in 2D is currently limited to Zxy.
    %   The 3D observation is full [Zxx, Zyx, Tx, Zxy, Zyy, Ty].

    assert(size(g, 1) == get_num_polarizations(dofmap.mesh.dim));

    % Precompute BC
    % TODO: It is expensive; could consider memoization in g?
    [bc_dofs, bc_vals] = build_dirichlet_bc(dofmap, g);

    % Preassemble terms independent of sigma
    Q = assemble_electrodes(dofmap, bc_dofs, electrode_coords);
    assemble_forward_operator_phase_2 = assemble_forward_operator_phase_1(dofmap, bc_dofs);

    % Init other functions
    solver_fwd = create_fwd_solver_direct(Q);
    [assemble_observation_only, assemble_observation_and_jacobian] = ...
        create_observation_assembler(dofmap, Q);

    % Extract dimensions
    n_freq = size(g, 2);
    n_ele = size(electrode_coords, 2);
    n_cells = dofmap.mesh.num_entities(dofmap.mesh.dim);
    n_comp_per_freq = get_num_impedance_components_per_frequency(dofmap.mesh.dim, n_ele);

    % Assign return value
    func = @compute_observation;


    function [d, J] = compute_observation(lambda1, lambda2, sigma)

        assert(numel(lambda1) == n_freq);
        assert(numel(lambda2) == n_freq);
        assert(all(size(sigma) == [n_cells, 1]));

        % Preassemble terms independent of frequency
        assemble_forward_operator_phase_3 = assemble_forward_operator_phase_2(sigma);

        % Preallocate return value(s)
        d = zeros(n_freq*n_comp_per_freq, 1);
        if nargout >= 2
            J = zeros(n_freq*n_comp_per_freq, n_cells);
        end

        % Loop over frequencies
        for j = 1:n_freq

            [A, b] = assemble_forward_operator_phase_3(lambda1(j), lambda2(j), bc_vals(:, :, j));
            if nargout == 1
                E = solver_fwd(A, b);
                clear('A', 'b');
                d_ = assemble_observation_only(E);
                d_ = transform_to_impedance(dofmap.mesh.dim, lambda1(j), d_);
                d(n_comp_per_freq*(j-1)+1:n_comp_per_freq*j, 1) = d_;
            else
                [E, q] = solver_fwd(A, b);
                clear('A', 'b');
                [d_, J_] = assemble_observation_and_jacobian(lambda1(j), sigma, E, q);
                [d_, J_] = transform_to_impedance(dofmap.mesh.dim, lambda1(j), d_, J_);
                d(n_comp_per_freq*(j-1)+1:n_comp_per_freq*j, 1) = d_;
                J(n_comp_per_freq*(j-1)+1:n_comp_per_freq*j, :) = J_;
            end

        end

    end
end


function func = create_fwd_solver_direct(Q)
    function [E, q] = solve(A, b)
        % Perform forward solve directly

        fprintf('Forward solve directly (factorization): ');
        t = tic();
        lu = decomposition(A, 'lu');
        fprintf('%f seconds\n', toc(t));

        fprintf('Forward solve directly (primal solve): ');
        t = tic();
        E = lu\b;
        assert(~issparse(E));
        fprintf('%f seconds\n', toc(t));

        if nargout < 2
            return
        end

        fprintf('Forward solve directly (dual solve): ');
        t = tic();
        q = lu\full(Q);
        assert(~issparse(q));
        fprintf('%f seconds\n', toc(t));
    end

    func = @solve;
end


function [bc_dofs, bc_vals] = build_dirichlet_bc(dofmap, g)
    fprintf('Build Dirichlet BC: ');
    t = tic();

    w = warning('off', 'Mesh:ExtraConnStored');  % Mute warning
    dofmap.mesh.compute_connectivity(dofmap.mesh.dim, dofmap.mesh.dim-1);
    warning(w);  % Restore warning
    dofmap.mesh.clear_connectivity(dofmap.mesh.dim-1, 0);
    dofmap.mesh.compute_boundary_facets();

    z_ = zeros(1, dofmap.mesh.dim);
    bc_dofs = assembling.build_dirichlet_dofs(dofmap, {@(x) true}, {@(x) z_});
    clear('z_');

    bc_vals = zeros(numel(bc_dofs), numel(g));
    for ij = 1:numel(g)
        [~, bc_vals(:, ij)] = assembling.build_dirichlet_dofs(dofmap, {@(x) true}, g(ij));
    end
    bc_vals = reshape(bc_vals, [numel(bc_dofs), size(g)]);

    fprintf('%f seconds\n', toc(t));
end


function Q = assemble_electrodes(dofmap, bc_dofs, electrode_coords)
    fprintf('Assemble electrodes: ');
    t = tic();

    num_electrodes = size(electrode_coords, 2);

    dofmap.mesh.init_geometric_queries();
    Q1 = assembling.assemble_point_sources(dofmap, electrode_coords);
    Q2 = assembling.assemble_point_sources(dofmap, electrode_coords, 'curl');
    dofmap.mesh.clear_geometric_queries();

    % NB: The following layout (sorting of columns in Q) has to play well
    %     with transform_to_impedance() below

    switch dofmap.mesh.dim
    case 2
        Q = sparse([], [], [], dofmap.dim, 2*num_electrodes);
        assert(size(Q1, 2) == 2*num_electrodes);
        assert(size(Q2, 2) == num_electrodes);
        Q(:, 1:2:end) = Q1(:, 1:2:end);  % Ex
        Q(:, 2:2:end) = Q2(:, :);        % curlE
    case 3
        Q = sparse([], [], [], dofmap.dim, 5*num_electrodes);
        assert(size(Q1, 2) == 3*num_electrodes);
        assert(size(Q2, 2) == 3*num_electrodes);
        colsl = sort([1:5:5*num_electrodes, 2:5:5*num_electrodes]);
        colsr = sort([1:3:3*num_electrodes, 2:3:3*num_electrodes]);
        Q(:, colsl) = Q1(:, colsr);  % Ex, Ey
        colsl = sort([3:5:5*num_electrodes, 4:5:5*num_electrodes, 5:5:5*num_electrodes]);
        Q(:, colsl) = Q2(:, :);  % curlEx, curlEy, curlEz
    otherwise
        error('not implemented');
    end
    clear('Q1', 'Q2', 'colsl', 'colsr');

    % Apply BCs
    [~, Q] = assembling.apply_dirichlet_bc([], Q, bc_dofs, 0);

    assert(issparse(Q));

    fprintf('%f seconds\n', toc(t));
end


function func = assemble_forward_operator_phase_1(dofmap, bc_dofs)
    fprintf('Assembling forward operator (phase 1): ');
    t = tic();

    % Preassemble parts of the operator independent of sigma
    A0 = assembling.assemble_curl_curl(dofmap, 0, 1);
    A2 = assembling.assemble_curl_curl(dofmap, -1, 1/0);

    function func = assemble_forward_operator_phase_2(sigma)
        fprintf('Assembling forward operator (phase 2): ');
        t = tic();

        % Preassemble the part with sigma
        A1 = assembling.assemble_curl_curl(dofmap, -sigma, 1/0);

        function [A, b] = assemble_forward_operator_phase_3(lambda1, lambda2, bc_vals)
            fprintf('Assembling forward operator (phase 3): ');
            t = tic();

            % Assemble the final forward system
            A = A0 + lambda1*A1 + lambda2*A2;
            b = zeros(dofmap.dim, size(bc_vals, 2));
            [A, b] = assembling.apply_dirichlet_bc(A, b, bc_dofs, bc_vals);

            fprintf('%f seconds\n', toc(t));
        end

        func = @assemble_forward_operator_phase_3;

        fprintf('%f seconds\n', toc(t));
    end

    func = @assemble_forward_operator_phase_2;

    fprintf('%f seconds\n', toc(t));
end


function [f1, f2] = create_observation_assembler(dofmap, Q)
    % First (linear) stage of computing MT observations/derivatives
    %
    % SYNOPSIS
    %
    %   [assemble_d, assemble_d_J] = create_observation_assembler(dofmap, Q);
    %   d = assemble_d(E);
    %   [d, J] = assemble_d_J(lambda1, sigma, E, q);

    function d = assemble_observation_only(E)
        fprintf('Assembling observation only: ');
        t = tic();

        d = reshape(Q.'*E, [], 1);
        assert(~issparse(d));

        fprintf('%f seconds\n', toc(t));
    end

    function [d, J] = assemble_observation_and_jacobian(lambda1, sigma, E, q)
        fprintf('Assembling observation and sensitivity: ');
        t = tic();

        % Assemble sensitivity matrix
        J = -lambda1 * assembling.assemble_mass_sensitivity(dofmap, E, q);

        % Assemble observables
        d = reshape(Q.'*E, [], 1);

        % Scale by sigma
        % FIXME: Transformation should happen elsewhere
        % FIXME: Which one is faster?
        %diag_sigma = spdiags(sigma, 0, numel(sigma), numel(sigma));
        %J = J*diag_sigma;
        J = J .* (sigma(:).');

        assert(~issparse(d));
        assert(~issparse(J));

        fprintf('%f seconds\n', toc(t));
    end

    f1 = @assemble_observation_only;
    f2 = @assemble_observation_and_jacobian;
end


function n = get_num_polarizations(dim)
    switch dim
    case 2
        n = 1;  % Currently implemented only H-polarization (Zxy)
    case 3
        n = 2;
    end
end


function n = get_num_impedance_components_per_frequency(dim, num_electrodes)
    switch dim
    case 2
        n = num_electrodes;  % Zxy
    case 3
        n = 6*num_electrodes;  % Zxx, Zyx, Tx, Zxy, Zyy, Ty
    end
end


function [d, J] = transform_to_impedance(dim, lambda1, d, J)
    % Transform linear observations/derivatives to impedace/tipper
    %
    % SYNOPSIS
    %
    %   ZT = transform_to_impedance(dim, lambda1, d);
    %   [ZT, dZT] = transform_to_impedance(dim, lambda1, d, J);
    %
    % NOTES
    %
    %   See the implementation for ordering of entries in d, J, ZT, dZT

    fprintf('Transform to impedance: ');
    t = tic();

    assert(nargin-2 == nargout);

    % NB: The minus sign given by Fourier transform convention;
    %     has to match sign in the underlying curl-curl equation
    k = -1/lambda1;

    switch dim
    case 2

        n = size(d, 1) / 2;
        assert(is_integer(n));

        [E, H] = deal(zeros(n, 1));
        E(:, 1) =   d(1:2:end);  % E1x
        H(:, 1) = k*d(2:2:end);  % H1y

        Zxy = E./H;
        d = reshape(Zxy, [], 1);

        if nargout == 1
            fprintf('%f seconds\n', toc(t));
            return
        end

        m = size(J, 2);

        [dEi, dHi] = deal(zeros(m, 1));
        dZxy = zeros(1*n, m);
        for i = 0:n-1
            dEi(:, 1) =   J(2*i+1, :);  % dE1x
            dHi(:, 1) = k*J(2*i+2, :);  % dH1y
            dZxy(i+1, :) = transpose((dEi-E(i+1, :).*dHi./H(i+1, :))./H(i+1, :));
        end
        J = dZxy;

    case 3

        [compute_Z, compute_dZ] = generate_code();
        n = size(d, 1) / 10;
        assert(is_integer(n));

        B = zeros(n, 6);
        H = zeros(n, 4);
        B(:, 1) =   d(    1:5: 5*n);  % E1x
        B(:, 2) =   d(    2:5: 5*n);  % E1y
        B(:, 3) = k*d(    5:5: 5*n);  % H1z
        B(:, 4) =   d(5*n+1:5:10*n);  % E2x
        B(:, 5) =   d(5*n+2:5:10*n);  % E2y
        B(:, 6) = k*d(5*n+5:5:10*n);  % H2z
        H(:, 1) = k*d(    3:5: 5*n);  % H1x
        H(:, 2) = k*d(    4:5: 5*n);  % H1y
        H(:, 3) = k*d(5*n+3:5:10*n);  % H2x
        H(:, 4) = k*d(5*n+4:5:10*n);  % H2y

        d = compute_Z(B, H);
        d = reshape(d, [], 1);

        if nargout == 1
            fprintf('%f seconds\n', toc(t));
            return
        end

        m = size(J, 2);

        dBi = zeros(m, 6);
        dHi = zeros(m, 4);
        J_ = zeros(6*n, m);
        for i = 0:n-1
            dBi(:, 1) =   J(5*i    +1, :);  % dE1x
            dBi(:, 2) =   J(5*i    +2, :);  % dE1y
            dBi(:, 3) = k*J(5*i    +5, :);  % dH1z
            dBi(:, 4) =   J(5*(n+i)+1, :);  % dE2x
            dBi(:, 5) =   J(5*(n+i)+2, :);  % dE2y
            dBi(:, 6) = k*J(5*(n+i)+5, :);  % dH2z
            dHi(:, 1) = k*J(5*i    +3, :);  % dH1x
            dHi(:, 2) = k*J(5*i    +4, :);  % dH1y
            dHi(:, 3) = k*J(5*(n+i)+3, :);  % dH2x
            dHi(:, 4) = k*J(5*(n+i)+4, :);  % dH2y
            J_(i+1:n:end, :) = transpose(compute_dZ(B(i+1, :), dBi, H(i+1, :), dHi));
        end
        J = J_;

    end

    fprintf('%f seconds\n', toc(t));
end


function [Z, dZ] = generate_code()

    % Input args of generated code (note that the generated
    % code is vectorized along the singleton dimension!)
    B_ = sym('E', [1, 6]);    % [E1x, E1y, H1z, E2x, E2y, H2z]
    H_ = sym('H', [1, 4]);    % [H1x, H1y, H2x, H2y]
    dB_ = sym('dE', [1, 6]);  % [dE1x, dE1y, dH1z, dE2x, dE2y, dH2z]
    dH_ = sym('dH', [1, 4]);  % [dH1x, dH1y, dH2x, dH2y]

    % Reshape args as needed matrices
    B = reshape(B_, [3, 2]);
    H = reshape(H_, [2, 2]);
    dB = reshape(dB_, [3, 2]);
    dH = reshape(dH_, [2, 2]);

    % Compute impedance/tipper (stacked) and their derivatives
    % and reshape back (to allow for vectorization)
    Z_ = reshape(B/H, [], 6);
    dZ_ = reshape((dB-B/H*dH)/H, [], 6);

    % Generate code
    Z = matlabFunction(Z_, 'Vars', {B_, H_});
    dZ = matlabFunction(dZ_, 'Vars', {B_, dB_, H_, dH_});

end


function flag = is_integer(num)
    flag = floor(num) == ceil(num);
end
