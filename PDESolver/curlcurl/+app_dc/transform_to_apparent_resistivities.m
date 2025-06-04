function [M1, M2] = transform_to_apparent_resistivities(electrode_coords, M1, M2)
    % Transform DC potential/current measurements into apparent resistivities
    %
    % INPUT/OUTPUT ARGUMENTS
    %   electrode_coords ... Matrix [dim, num_electrodes]
    %   M1, M2           ... Matrix [num_electrodes, num_observations]
    %
    % REMARKS
    %   Input M1 and M2 are supposed to have at most one +1
    %   and at most -1 in each column. These corresponds to
    %   dipole or pole measurements/sources. The output M1
    %   and M2 correspond to apparent resistivity for 2D
    %   or 3D half-space.
    %
    %   Note that pole source/measurement will not work in 2D,
    %   because an electrode at infinity corresponds to infinite
    %   fundamental solution (log(r)).

    [num_electrodes, num_observations] = size(M1);
    assert(all(size(M1) == size(M2)));
    assert(size(electrode_coords, 2) == num_electrodes);

    dim = size(electrode_coords, 1);
    switch dim
    case 2
        % FIXME: Wrong code-path for 2.5D?
        % FIXME: Should we have different config factor for points which are not on boundary??
        config_factor = @config_factor_2d_halfspace;
    case 3
        % FIXME: Should we have different config factor for points which are not on boundary??
        config_factor = @config_factor_3d_halfspace;
    end

    for k = 1:num_observations
        [a, b] = extract_two_electrode_indices(M1(:, k));
        [m, n] = extract_two_electrode_indices(M2(:, k));
        xa = extract_electrode_coords(electrode_coords, a);
        xb = extract_electrode_coords(electrode_coords, b);
        xm = extract_electrode_coords(electrode_coords, m);
        xn = extract_electrode_coords(electrode_coords, n);
        scale = sqrt(abs(config_factor(xa, xb, xm, xn)));
        M1(:, k) = scale*M1(:, k);
        M2(:, k) = scale*M2(:, k);
    end
end


function [a, b] = extract_two_electrode_indices(electrodes)
   a = find(electrodes == 1);
   b = find(electrodes == -1);
   if isempty(a) && isempty(b) || numel(a) > 1 || numel(b) > 1 || nnz(electrodes) > numel(a) + numel(b)
       error('Don''t know how to handle electrode setup %s', mat2str(electrodes));
   end
end


function coord = extract_electrode_coords(electrode_coords, index)
    if isempty(index)
        coord = inf;
    else
        coord = electrode_coords(:, index);
    end
end


function k = config_factor_3d_halfspace(xa, xb, xm, xn)
    % NB: If one of xa, xb is inf and/or one of xm, xn is inf,
    %     this formula works just fine
    k = 2*pi/( ...
        + 1/vecnorm(xa-xm) - 1/vecnorm(xb-xm) ...
        - 1/vecnorm(xa-xn) + 1/vecnorm(xb-xn) ...
    );
end


function k = config_factor_2d_halfspace(xa, xb, xm, xn)
    if isinfvec(xa)
        assert(~isinfvec(xb));
        if isinfvec(xm)
            assert(~isinfvec(xn));
            k = +config_factor_2d_halfspace_pole_pole(xb, xn);
        elseif isinfvec(xn)
            k = -config_factor_2d_halfspace_pole_pole(xb, xm);
        else
            k = -config_factor_2d_halfspace_pole_dipole(xb, xm, xn);
        end
    elseif isinfvec(xb)
        if isinfvec(xm)
            assert(~isinfvec(xn));
            k = -config_factor_2d_halfspace_pole_pole(xa, xn);
        elseif isinfvec(xn)
            k = +config_factor_2d_halfspace_pole_pole(xa, xm);
        else
            k = +config_factor_2d_halfspace_pole_dipole(xa, xm, xn);
        end
    else
        if isinfvec(xm)
            assert(~isinfvec(xn));
            k = -config_factor_2d_halfspace_pole_dipole(xn, xa, xb);
        elseif isinfvec(xn)
            k = +config_factor_2d_halfspace_pole_dipole(xm, xa, xb);
        else
            k = +config_factor_2d_halfspace_dipole_dipole(xa, xb, xm, xn);
        end
    end
end


function flag = isinfvec(x)
    flag = any(isinf(x));
end


function k = config_factor_2d_halfspace_dipole_dipole(xa, xb, xm, xn)
    k = pi/( ...
        - log(vecnorm(xa-xm)) + log(vecnorm(xb-xm)) ...
        + log(vecnorm(xa-xn)) - log(vecnorm(xb-xn)) ...
    );
end


function k = config_factor_2d_halfspace_pole_dipole(xa, xm, xn)
    k = pi/( ...
        - log(vecnorm(xa-xm)) ...
        + log(vecnorm(xa-xn)) ...
    );
end


function k = config_factor_2d_halfspace_pole_pole(xa, xm)
    k = pi/( ...
        - log(vecnorm(xa-xm)) ...
    );
end
