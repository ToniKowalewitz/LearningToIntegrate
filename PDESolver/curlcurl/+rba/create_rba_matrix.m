function r = create_rba_matrix(n, A, M, y0, solver)
    % Return rational best approximation of exp(-t*inv(M)*A)*y0
    %
    % Vector-valued expression
    %
    %   y(t) = exp(-t*inv(M)*A)*y0
    %
    % with t > 0, A and M HPD matrices, y0 vector, is a solution
    % of ordinary differential system
    %
    %   M*y' + A*y = 0    on (0, +Inf),
    %         y(0) = y0.
    %
    % This routine returns rational best approximation of y(t).
    %
    % INPUT PARAMETER
    %   n      ... number of poles; higher the number higher the precision
    %   A      ... HPD matrix
    %   M      ... HPD matrix
    %   y0     ... vector, initial condition
    %   solver ... function handle of signature
    %
    %                  x = solver(T, b)
    %
    %              which returns (approximate) solution of equation Tx=b;
    %              in the simplest case:
    %
    %                  solver = @(T, b) T\b
    %
    % OUTPUT PARAMETER
    %   y ... a function handle; y(t) approximates exp(-t*inv(M)*A)*y0

    real_case = isreal(A) && isreal(M);
    % TODO: Should we check symmetry?

    [poles, coeffs, coeff0] = rba.get_poles(n, 'include_conjugate', ~real_case);
    n = numel(poles);

    if real_case
        % We use half of the paired poles so need to extract real part below
        real_ = @real;
    else
        % No-op
        real_ = @(x) x;
    end

    % Precompute constant, repeatedly used vector
    b = M*y0;

    function y = r_(t)
        assert(ismatrix(A) && ismatrix(b) && size(A, 1) == size(b, 1));

        % Pole at infinity
        y = coeff0*y0;

        % Finite poles (solve shifted systems)
        for i = 1:n
            A_shift = t*A - poles(i)*M;
            y = y + real_( coeffs(i) * solver(A_shift, b) );
        end
    end

    r = @r_;
end
