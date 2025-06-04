function element = create_raviart_thomas_element(dim, order)
    % Create instance of Raviart Thomas element.
    %
    % SYNTAX
    %   element = create_raviart_thomas_element(dim, order)
    %
    % INPUT/OUTPUT PARAMETERS
    %
    %   dim     ... Dimension of the reference cell.
    %   order   ... Order of the element.
    %   element ... Instance of the HdivElement class.
    %
    % REMARKS
    %
    %   Reference simplex of dimension dim is given by
    %   vertex coordinates
    %
    %     [eye(dim), zeros(dim, 1)]
    %
    %   and lexicographic entity-vertex numbering, e.g.,
    %   edges of reference tetrahedron are given
    %   by vertices
    %
    %     [1,2],[1,3],[1,4],[2,3],[2,4],[3,4],
    %
    %   faces of reference tetrahedron are given
    %   by vertices
    %
    %     [1,2,3],[1,2,4],[1,3,4],[2,3,4],
    %
    %   which results in face-edge connectivity
    %
    %     [1,2,4],[1,3,5],[2,3,6],[4,5,6].
    %
    %   FIXME: We should move these assumptions from documentation
    %   to actual code and document it better.

    reference_simplex = fe.ReferenceSimplex(dim);

    switch dim
    case 2
        switch order

        case 1
            % Triangle, Raviart Thomas order 1
            entity_dofs = {
                reshape([],  0, 3);  % vertex dofs
                reshape(1:3, 1, 3);  % edge dofs
                reshape([],  0, 1);  % cell dofs
            };
            tabulate_basis_handle = @tabulate_basis_triangle_order_1;
            dual_basis_handle = @dual_basis_triangle_order_1;

        case 2
            % Triangle, Raviart Thomas order 2
            entity_dofs = {
                reshape([],  0, 3);  % vertex dofs
                reshape(1:6, 2, 3);  % edge dofs
                reshape(7:8, 2, 1);  % cell dofs
            };
            tabulate_basis_handle = @tabulate_basis_triangle_order_2;
            dual_basis_handle = @dual_basis_triangle_order_2;

        otherwise
          error('Order %d not implemented', order);
        end

    case 3
        switch order

        case 1
            % Tetrahedron, Raviart Thomas order 1
            entity_dofs = {
                reshape([],  0, 4);  % vertex dofs
                reshape([],  0, 6);  % edge dofs
                reshape(1:4, 1, 4);  % face dofs
                reshape([],  0, 1);  % cell dofs
            };
            tabulate_basis_handle = @tabulate_basis_tetrahedron_order_1;
            dual_basis_handle = @dual_basis_tetrahedron_order_1;

        case 2
            % Tetrahedron, Raviart Thomas order 2
            entity_dofs = {
                reshape([],    0, 4);  % vertex dofs
                reshape([],    0, 6);  % edge dofs
                reshape( 1:12, 3, 4);  % face dofs
                reshape(13:15, 3, 1);  % cell dofs
            };
            tabulate_basis_handle = @tabulate_basis_tetrahedron_order_2;
            dual_basis_handle = @dual_basis_tetrahedron_order_2;

        otherwise
            error('Order %d not implemented', order);
        end

    otherwise
        error('Dimension %d not implemented', dim);
    end

    % Functions are vector-valued
    value_shape = dim;

    % Build nodal element from its defining data
    % NB numbering in entity_dofs and facet_dofs must match
    % ordering in dual basis. This is not fool-proof implementation...
    % On the other hand the numbering in basis does not matter,
    % only important is that the basis spans the right space.
    element = fe.HdivElement('RaviartThomas', 'contravariant', order, ...
                             reference_simplex, value_shape, entity_dofs,  ...
                             tabulate_basis_handle, dual_basis_handle);
end


function Phi = tabulate_basis_triangle_order_1(p)
    assert(isrow(p));
    Phi = [
        p(1)-1, p(2)  ;
        p(1)  , p(2)-1;
        p(1)  , p(2)  ;
    ];
end


function Phi = tabulate_basis_triangle_order_2(p)  %#ok
    error('not implemented');
    assert(isrow(p));  %#ok
    Phi = [
    ];
end


function dofs = dual_basis_triangle_order_1(v)
    dof = @(p, w) v(p)*w;
    dofs = [
        % normal evaluations at edge midpoint times edge length
        dof([0.5, 0.5], [-1;-1]);  % e1
        dof([0.5, 0.0], [ 0;-1]);  % e2
        dof([0.0, 0.5], [ 1; 0]);  % e3
    ];
end


function dofs = dual_basis_triangle_order_2(v)
    t1 = 1/3;

    % Gauss points on [0,1]
    s1 = 0.21132486540518713; s2 = 0.78867513459481287;

    % [q1,q1], [q1,q2], [q2,q1] is degree 2 quadrature on
    % triangle [1,0],[0,1],[0,0] with weights 1/6;
    % see [Strang, Fix]
    q1 = 1/6; q2 = 2/3;

    dof = @(p, w) v(p)*w;
    dofs = [
        % normal evaluations at edge Gauss points
        % times edge length
        dof([s2, s1], [-1;-1]);  % e1
        dof([s1, s2], [-1;-1]);  % e1
        dof([s2,  0], [ 0;-1]);  % e2
        dof([s1,  0], [ 0;-1]);  % e2
        dof([ 0, s2], [ 1; 0]);  % e3
        dof([ 0, s1], [ 1; 0]);  % e3

        % directional evaluations at cell quad points in
        % direction of first and second edge (by lexicographic
        % numbering) times edge length;
        t1*(dof([q1, q1], [-1; 1]) + dof([q1, q2], [-1; 1]) + dof([q2, q1], [-1; 1]));  % c1
        t1*(dof([q1, q1], [-1; 0]) + dof([q1, q2], [-1; 0]) + dof([q2, q1], [-1; 0]));  % c1
    ];
end


function Phi = tabulate_basis_tetrahedron_order_1(p)
    assert(isrow(p));
    Phi = [
        p(1)-1, p(2)  , p(3)  ;
        p(1)  , p(2)-1, p(3)  ;
        p(1)  , p(2)  , p(3)-1;
        p(1)  , p(2)  , p(3)  ;
    ];
end


function Phi = tabulate_basis_tetrahedron_order_2(p)  %#ok
    error('not implemented');
    assert(isrow(p));  %#ok
    Phi = [
    ];
end


function dofs = dual_basis_tetrahedron_order_1(v)
    t1 = 1/3;
    dof = @(p, w) v(p)*w;
    dofs = [
        % normal evaluations at face midpoint times face area
        dof([t1, t1, t1], [ 0.5; 0.5; 0.5]);  % f1
        dof([t1, t1,  0], [   0;   0; 0.5]);  % f2
        dof([t1,  0, t1], [   0;-0.5;   0]);  % f3
        dof([ 0, t1, t1], [ 0.5;   0;   0]);  % f4
    ];
end


function dofs = dual_basis_tetrahedron_order_2(v)  %#ok
    error('not implemented');
    dof = @(p, w) v(p)*w;  %#ok
    dofs = [
    ];
end
