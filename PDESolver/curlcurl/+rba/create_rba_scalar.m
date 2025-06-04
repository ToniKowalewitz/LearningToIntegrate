function r = create_rba_scalar(n)
    % Return rational best approximation of exp(-x) on (0, +INF)
    %
    % INPUT PARAMETER
    %   n ... number of poles; higher the number higher the precision
    %
    % OUTPUT PARAMETER
    %   r ... a function handle; r(x) approximates exp(x) uniformly
    %         on (0, +INF); r is vectorized, it takes array of any
    %         shape, computes r(x) entry-wise, and returns array of
    %         the same shape

    % Get poles and coefficients for complex case and real case
    [c_poles, c_coeffs, c_absterm] = rba.get_poles(n, 'include_conjugate', true);
    [r_poles, r_coeffs, r_absterm] = rba.get_poles(n, 'include_conjugate', false);

    % Sanity check
    assert(iscolumn(c_poles) && iscolumn(c_coeffs));
    assert(iscolumn(r_poles) && iscolumn(r_coeffs));
    assert(isscalar(c_absterm) && isreal(c_absterm));
    assert(isscalar(r_absterm) && isreal(r_absterm));

    function y = r_(x)
        % Reshape to row vector
        shp = size(x);
        x = x(:).';

        % Compute RBA values for all xs
        if isreal(x)
            y = sum(real(r_coeffs ./ (x - r_poles)), 1) + r_absterm;
        else
            y = sum(c_coeffs ./ (x - c_poles), 1) + c_absterm;
        end

        % Reshape back
        y = reshape(y, shp);
    end

    % Assign return value
    r = @r_;
end
