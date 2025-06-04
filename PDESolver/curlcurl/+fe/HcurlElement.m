classdef HcurlElement < fe.FiniteElement
    % HcurlElement class describing H(curl)-conforming finite element.

    methods (Access = public)

        function obj = HcurlElement(family, mapping, order, simplex, ...
                                    value_shape, entity_dofs, ...
                                    tabulate_basis_handle, dual_basis_handle)
            % Constructor.

            % Call superclass constructor
            obj = obj@fe.FiniteElement(family, mapping, order, simplex, ...
                                       value_shape, entity_dofs, ...
                                       tabulate_basis_handle, dual_basis_handle);

            % Compute basis curls
            obj.compute_basis_curl(tabulate_basis_handle);

        end

        function CurlPhi = tabulate_basis_curl(obj, point)
            % Return value of basis members curls at point.

            % Take linear combination of original basis curls for nodality
            CurlPhi = obj.V_*obj.tabulate_basis_curl_(point);
        end

    end

    methods (Access = private, Sealed = true)

        function compute_basis_curl(obj, tabulate_basis_handle)
            % Symbolically compute curls of original basis

            % Define curl in 2d and 3d
            % NB: avoiding matlab's sym/curl as it does not operate on matrices
            if obj.simplex.dim == 2
                curl_ = @(f, x) diff(f(:, 2), x(1)) - diff(f(:, 1), x(2));
            elseif obj.simplex.dim == 3
                curl_ = @(f, x) [
                    diff(f(:, 3), x(2)) - diff(f(:, 2), x(3)), ...
                    diff(f(:, 1), x(3)) - diff(f(:, 3), x(1)), ...
                    diff(f(:, 2), x(1)) - diff(f(:, 1), x(2)), ...
                ];
            end

            % Compute basis curls symbolically
            x = sym('x', [1, obj.simplex.dim], 'real');
            tabulate_basis_curl_sym = curl_(tabulate_basis_handle(x), x);

            % Convert symbolic expression into to fast function
            tabulate_basis_curl_fun = matlabFunction(tabulate_basis_curl_sym, 'Vars', {x});

            % Store the result
            obj.tabulate_basis_curl_ = tabulate_basis_curl_fun;
        end

    end

    properties (Access = private)
        tabulate_basis_curl_   % Handle with curls of original basis
    end

end
