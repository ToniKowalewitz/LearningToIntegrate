classdef HdivElement < fe.FiniteElement
    % HdivElement class describing H(div)-conforming finite element.

    methods (Access = public)

        function obj = HdivElement(family, mapping, order, simplex, ...
                                   value_shape, entity_dofs, ...
                                   tabulate_basis_handle, dual_basis_handle)
            % Constructor.

            % Call superclass constructor
            obj = obj@fe.FiniteElement(family, mapping, order, simplex, ...
                                       value_shape, entity_dofs, ...
                                       tabulate_basis_handle, dual_basis_handle);

            % Compute basis divs
            obj.compute_basis_div(tabulate_basis_handle);

        end

        function DivPhi = tabulate_basis_div(obj, point)
            % Return value of basis members divs at point.

            % Take linear combination of original basis divs for nodality
            DivPhi = obj.V_*obj.tabulate_basis_div_(point);
        end

    end

    methods (Access = private, Sealed = true)

        function compute_basis_div(obj, tabulate_basis_handle)
            % Symbolically compute divs of original basis

            % Define divergence in 2d and 3d
            % NB: avoiding matlab's sym/divergence as it does not operate on matrices
            if obj.simplex.dim == 2
                div_ = @(f, x) diff(f(:, 1), x(1)) + diff(f(:, 2), x(2));
            elseif obj.simplex.dim == 3
                div_ = @(f, x) diff(f(:, 1), x(1)) + diff(f(:, 2), x(2)) + diff(f(:, 3), x(3));
            end

            % Compute basis divergence symbolically
            x = sym('x', [1, obj.simplex.dim], 'real');
            tabulate_basis_div_sym = div_(tabulate_basis_handle(x), x);

            % Convert symbolic expression into to fast function
            tabulate_basis_div_fun = matlabFunction(tabulate_basis_div_sym, 'Vars', {x});

            % Store the result
            obj.tabulate_basis_div_ = tabulate_basis_div_fun;
        end

    end

    properties (Access = private)
        tabulate_basis_div_   % Handle with divs of original basis
    end

end
