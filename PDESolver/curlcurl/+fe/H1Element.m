classdef H1Element < fe.FiniteElement
    % H1Element class describing H^1-conforming finite element.

    methods (Access = public)

        function obj = H1Element(family, mapping, order, simplex, ...
                                 value_shape, entity_dofs, ...
                                 tabulate_basis_handle, dual_basis_handle)
            % Constructor.

            % Call superclass constructor
            obj = obj@fe.FiniteElement(family, mapping, order, simplex, ...
                                       value_shape, entity_dofs, ...
                                       tabulate_basis_handle, dual_basis_handle);

            % Compute basis gradients
            obj.compute_basis_grad(tabulate_basis_handle);

        end

        function GradPhi = tabulate_basis_grad(obj, point)
            % Return value of basis members grads at point.

            % Take linear combination of original basis grads for nodality
            GradPhi = obj.V_*obj.tabulate_basis_grad_(point);
        end

    end

    methods (Access = private, Sealed = true)

        function compute_basis_grad(obj, tabulate_basis_handle)
            % Symbolically compute gradients of original basis

            % Compute basis gradients symbolically
            x = sym('x', [1, obj.simplex.dim], 'real');
            tabulate_basis_grad_sym = jacobian(tabulate_basis_handle(x), x);

            % Convert symbolic expression into to fast function
            tabulate_basis_grad_fun = matlabFunction(tabulate_basis_grad_sym, 'Vars', {x});

            % Store the result
            obj.tabulate_basis_grad_ = tabulate_basis_grad_fun;
        end

    end

    properties (Access = private)
        tabulate_basis_grad_   % Handle with grads of original basis
    end

end
