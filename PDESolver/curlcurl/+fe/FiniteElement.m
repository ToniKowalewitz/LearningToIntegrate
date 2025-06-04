classdef FiniteElement < handle
    % FiniteElement class describing a finite element.

    properties (GetAccess = public, SetAccess = private)
        family        % chararray describing element family
        mapping       % chararray describing transformation
        order         % order of element
        simplex       % reference simplex
        value_shape   % value shape of FE functions
        fe_space_dim  % cardinality of basis
        facet_dofs    % dof indices for each facet
    end

    methods (Access = public)

        function obj = FiniteElement(family, mapping, order, simplex, ...
                                     value_shape, entity_dofs, ...
                                     tabulate_basis_handle, dual_basis_handle)
            % Constructor.

            % Stored passed data
            obj.family = family;
            obj.mapping = mapping;
            obj.order = order;
            obj.simplex = simplex;
            obj.value_shape = value_shape;
            obj.entity_dofs_ = entity_dofs;

            % Compute element dimension and facet dofs
            obj.compute_topology();

            % Compute nodal basis
            obj.compute_nodal_basis(tabulate_basis_handle, dual_basis_handle);

        end

        function Phi = tabulate_basis(obj, point)
            % Return value of basis members at point.

            % Take linear combination of original basis for nodality
            Phi = obj.V_*obj.tabulate_basis_(point);
        end

        function dofs = evaluate_dual_basis(obj, v)
            % Return value of dual basis members on function.
            %
            % Provides interpolation operator, e.g., evaluation
            % at lattice points for Lagrange element.

            dofs = obj.dual_basis_(v);
        end

        function dofs = get_entity_dofs(obj, dim)
            % Return dof indices for given entities of given dimension.
            %
            % Result is a matrix with column indices indexing entities
            % of the given dimension. Hence number of rows is number
            % of dofs per entity.

            dofs = obj.entity_dofs_{dim+1};
        end

        function dims = get_dof_entity_dims(obj)
            % Returns vector with dimensions of entities supporting dofs.

            dims = [];
            for d = 0:obj.simplex.dim
                if not(isempty(obj.get_entity_dofs(d)))
                    dims = horzcat(dims, d);  %#ok<AGROW>
                end
            end
        end

    end

    methods (Access = private, Sealed = true)

        function compute_topology(obj)  %#ok<*PROP>
            %% Compute dimension
            fe_space_dim = 0;
            for edim = 0:obj.simplex.dim
                fe_space_dim = fe_space_dim + numel(obj.get_entity_dofs(edim));
            end

            % Sannity check (check dofs are range 1:fe_space_dim)
            dofs = zeros(0, 1, 'uint32');
            for edim = 0:obj.simplex.dim
                edofs = obj.get_entity_dofs(edim);
                dofs = vertcat(dofs, edofs(:));  %#ok<AGROW>
            end
            assert(all(sort(dofs) == (1:fe_space_dim).'));

            %% Compute facet dofs
            num_facets = obj.simplex.num_entities(obj.simplex.dim-1);
            facet_dofs = cell(num_facets, 1);

            % Iterate over entity dimensions up to facet dimension
            for edim = 0:obj.simplex.dim-1
                edofs = obj.get_entity_dofs(edim);

                % Get connections of facet to dof entities
                facet_entity_conn = obj.simplex.get_connectivity(obj.simplex.dim-1, edim);

                % Iterate over facets and process dofs which map to the facet
                for f = 1:num_facets
                    fdofs = edofs(:, facet_entity_conn(:, f));
                    fdofs = reshape(fdofs, 1, []);
                    facet_dofs{f} = [facet_dofs{f}, fdofs];
                end
            end
            facet_dofs = cell2mat(facet_dofs);

            % Store computed data
            obj.fe_space_dim = fe_space_dim;
            obj.facet_dofs = facet_dofs;
        end

        function compute_nodal_basis(obj, tabulate_basis_handle, dual_basis_handle)
            % Initialize the object given handles to function tabulating
            % basis member at point and dual basis members on function.
            % Given dual basis is kept and primal basis is modified so that
            % the element is nodal, i.e., action of dual basis on primal
            % basis gives identity matrix.

            % Compute action of dual basis on original space basis and invert it
            V = dual_basis_handle(tabulate_basis_handle);
            V = reshape(V, [obj.fe_space_dim, obj.fe_space_dim]);
            V = inv(V);

            % Store all data
            obj.tabulate_basis_ = tabulate_basis_handle;
            obj.dual_basis_ = dual_basis_handle;
            obj.V_ = V;
        end

    end

    properties (Access = private)
        entity_dofs_           % dof indices for each entity
        tabulate_basis_        % handle with original basis
        dual_basis_            % handle with original dual basis
    end

    properties (GetAccess = protected, SetAccess = private)
        V_                     % matrix giving modified (nodal) basis
    end

end
