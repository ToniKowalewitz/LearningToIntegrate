function element = create_nedelec_element(dim, order)
    % Create instance of Nedelec element.
    %
    % SYNTAX
    %   element = create_nedelec_element(dim, order)
    %
    % INPUT/OUTPUT PARAMETERS
    %
    %   dim     ... Dimension of the reference cell.
    %   order   ... Order of the element.
    %   element ... Instance of the HcurlElement class.
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
            % Triangle, Nedelec order 1
            entity_dofs = {
                reshape([],  0, 3);  % vertex dofs
                reshape(1:3, 1, 3);  % edge dofs
                reshape([],  0, 1);  % cell dofs
            };
            tabulate_basis_handle = @tabulate_basis_triangle_order_1;
            dual_basis_handle = @dual_basis_triangle_order_1;

        case 2
            % Triangle, Nedelec order 2
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
            % Tetrahedron, Nedelec order 1
            entity_dofs = {
                reshape([],  0, 4);  % vertex dofs
                reshape(1:6, 1, 6);  % edge dofs
                reshape([],  0, 4);  % face dofs
                reshape([],  0, 1);  % cell dofs
            };
            tabulate_basis_handle = @tabulate_basis_tetrahedron_order_1;
            dual_basis_handle = @dual_basis_tetrahedron_order_1;

        case 2
            % Tetrahedron, Nedelec order 2
            entity_dofs = {
                reshape([],    0, 4);  % vertex dofs
                reshape( 1:12, 2, 6);  % edge dofs
                reshape(13:20, 2, 4);  % face dofs
                reshape([],    0, 1);  % cell dofs
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
    element = fe.HcurlElement('NedelecKind1', 'covariant', order, ...
                              reference_simplex, value_shape, entity_dofs,  ...
                              tabulate_basis_handle, dual_basis_handle);
end


function Phi = tabulate_basis_triangle_order_1(p)
    assert(isrow(p));
    Phi = [
        -p(2), p(1);
            1,    0;
            0,    1;
    ];
end


function Phi = tabulate_basis_triangle_order_2(p)
    assert(isrow(p));
    x = p(1); y = p(2);
    Phi = [
          1,    0;
          0,    1;
          x,    0;
          0,    x;
          y,    0;
          0,    y;
        -x*y, x*x;
        -y*y, x*y;
    ];
end


function dofs = dual_basis_triangle_order_1(v)
    dof = @(p, w) v(p)*w;
    dofs = [
        % tangential evaluations at edge midpoint (in direction of
        % lower-to-higher-vertex-index) times edge length
        dof([0.5, 0.5], [-1; 1]);  % e1
        dof([0.5, 0.0], [-1; 0]);  % e2
        dof([0.0, 0.5], [ 0;-1]);  % e3
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
        % directional evaluations at edge Gauss points
        % (in direction of lower-to-higher-vertex-index)
        % times edge length; [Monk (5.49)]
        dof([s2, s1], [-1; 1]);  % e1
        dof([s1, s2], [-1; 1]);  % e1
        dof([s2,  0], [-1; 0]);  % e2
        dof([s1,  0], [-1; 0]);  % e2
        dof([ 0, s2], [ 0;-1]);  % e3
        dof([ 0, s1], [ 0;-1]);  % e3

        % directional evaluations at cell quad points in
        % direction of first and second edge (by lexicographic
        % numbering) times edge length;
        t1*(dof([q1, q1], [-1; 1]) + dof([q1, q2], [-1; 1]) + dof([q2, q1], [-1; 1]));  % c1
        t1*(dof([q1, q1], [-1; 0]) + dof([q1, q2], [-1; 0]) + dof([q2, q1], [-1; 0]));  % c1
    ];
end


function Phi = tabulate_basis_tetrahedron_order_1(p)
    assert(isrow(p));
    x = p(1);
    y = p(2);
    z = p(3);
    Phi = [
           -y,     x,     0;  % e1
           -z,     0,     x;  % e2
        y+z-1,    -x,    -x;  % e3
            0,    -z,     y;  % e4
           -y, x+z-1,    -y;  % e5
           -z,    -z, x+y-1;  % e6
    ];
end


function Phi = tabulate_basis_tetrahedron_order_2(p)
    assert(isrow(p));
    x = p(1); y = p(2); z = p(3);
    w = x + y + z - 1;
    Phi = [
              -x*y,        x*x,          0;  % e1
              -y*y,        x*y,          0;  % e1
              -x*z,          0,        x*x;  % e2
              -z*z,          0,        x*z;  % e2
         x*(y+z-1),       -x*x,       -x*x;  % e3
        -w*(y+z-1),        x*w,        x*w;  % e3
                 0,       -y*z,        y*y;  % e4
                 0,       -z*z,        y*z;  % e4
              -y*y,  y*(x+z-1),       -y*y;  % e5
               y*w, -w*(x+z-1),        y*w;  % e5
              -z*z,       -z*z,  z*(x+y-1);  % e6
               z*w,        z*w, -w*(x+y-1);  % e6
              -y*z,        x*z,          0;  % f1
              -y*z,          0,        x*y;  % f1
               y*w,       -x*w,          0;  % f2
         y*(y+z-1),       -x*y,       -x*y;  % f2
               z*w,          0,       -x*w;  % f3
         z*(y+z-1),       -x*z,       -x*z;  % f3
                 0,        z*w,       -y*w;  % f4
              -y*z,  z*(x+z-1),       -y*z;  % f4
    ];
end


function dofs = dual_basis_tetrahedron_order_1(v)
    dof = @(p, w) v(p)*w;
    dofs = [
        % tangential evaluations at edge midpoint (in direction of
        % lower-to-higher-vertex-index) times edge length
        dof([0.5, 0.5, 0.0], [-1; 1; 0]);  % e1
        dof([0.5, 0.0, 0.5], [-1; 0; 1]);  % e2
        dof([0.5, 0.0, 0.0], [-1; 0; 0]);  % e3
        dof([0.0, 0.5, 0.5], [ 0;-1; 1]);  % e4
        dof([0.0, 0.5, 0.0], [ 0;-1; 0]);  % e5
        dof([0.0, 0.0, 0.5], [ 0; 0;-1]);  % e6
    ];
end


function dofs = dual_basis_tetrahedron_order_2(v)
    t1 = 1/3;

    % Gauss points on [0,1]
    s1 = 0.21132486540518713; s2 = 0.78867513459481287;

    % [q1,q1], [q1,q2], [q2,q1] is degree 2 quadrature on
    % triangle [1,0],[0,1],[0,0] with weights 1/6;
    % see [Strang, Fix]
    q1 = 1/6; q2 = 2/3;

    dof = @(p, w) v(p)*w;
    dofs = [
        % directional evaluations at edge Gauss points
        % (in direction of lower-to-higher-vertex-index)
        % times edge length; [Monk (5.49)]
        dof([s2, s1,  0], [-1; 1; 0]);  % e1
        dof([s1, s2,  0], [-1; 1; 0]);  % e1
        dof([s2,  0, s1], [-1; 0; 1]);  % e2
        dof([s1,  0, s2], [-1; 0; 1]);  % e2
        dof([s2,  0,  0], [-1; 0; 0]);  % e3
        dof([s1,  0,  0], [-1; 0; 0]);  % e3
        dof([ 0, s2, s1], [ 0;-1; 1]);  % e4
        dof([ 0, s1, s2], [ 0;-1; 1]);  % e4
        dof([ 0, s2,  0], [ 0;-1; 0]);  % e5
        dof([ 0, s1,  0], [ 0;-1; 0]);  % e5
        dof([ 0,  0, s2], [ 0; 0;-1]);  % e6
        dof([ 0,  0, s1], [ 0; 0;-1]);  % e6

        % tangential evaluations at facet quad points in
        % direction of first and second edge (by lexicographic
        % numbering) times edge length;
        t1*(dof([q1, q1, q2], [-1; 1; 0]) + dof([q1, q2, q1], [-1; 1; 0]) + dof([q2, q1, q1], [-1; 1; 0]));  % f1
        t1*(dof([q1, q1, q2], [-1; 0; 1]) + dof([q1, q2, q1], [-1; 0; 1]) + dof([q2, q1, q1], [-1; 0; 1]));  % f1
        t1*(dof([q1, q1,  0], [-1; 1; 0]) + dof([q1, q2,  0], [-1; 1; 0]) + dof([q2, q1,  0], [-1; 1; 0]));  % f2
        t1*(dof([q1, q1,  0], [-1; 0; 0]) + dof([q1, q2,  0], [-1; 0; 0]) + dof([q2, q1,  0], [-1; 0; 0]));  % f2
        t1*(dof([q1,  0, q1], [-1; 0; 1]) + dof([q1,  0, q2], [-1; 0; 1]) + dof([q2,  0, q1], [-1; 0; 1]));  % f3
        t1*(dof([q1,  0, q1], [-1; 0; 0]) + dof([q1,  0, q2], [-1; 0; 0]) + dof([q2,  0, q1], [-1; 0; 0]));  % f3
        t1*(dof([ 0, q1, q1], [ 0;-1; 1]) + dof([ 0, q1, q2], [ 0;-1; 1]) + dof([ 0, q2, q1], [ 0;-1; 1]));  % f4
        t1*(dof([ 0, q1, q1], [ 0;-1; 0]) + dof([ 0, q1, q2], [ 0;-1; 0]) + dof([ 0, q2, q1], [ 0;-1; 0]));  % f4
    ];
end
