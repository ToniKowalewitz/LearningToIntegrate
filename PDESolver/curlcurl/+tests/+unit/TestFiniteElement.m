classdef TestFiniteElement < matlab.unittest.TestCase


    properties (TestParameter)
        factory = {
            %#ok<*NBRAK>
            {nedelec,  2, 1, [2,]};
            {nedelec,  2, 2, [2,]};
            {nedelec,  3, 1, [3,]};
            {nedelec,  3, 2, [3,]};
            {rt,       2, 1, [2,]};
            {rt,       3, 1, [3,]};
            {lagrange, 1, 1, [1,]};
            {lagrange, 1, 2, [1,]};
            {lagrange, 1, 3, [1,]};
            {lagrange, 2, 1, [1,]};
            {lagrange, 2, 2, [1,]};
            {lagrange, 2, 3, [1,]};
            {lagrange, 2, 4, [1,]};
            {lagrange, 3, 1, [1,]};
            {lagrange, 3, 2, [1,]};
            {lagrange, 3, 3, [1,]};
            {lagrange, 3, 4, [1,]};
            {lagrange, 3, 5, [1,]};
            {dg      , 1, 0, [1,]};
            {dg      , 1, 1, [1,]};
            {dg      , 1, 2, [1,]};
            {dg      , 1, 3, [1,]};
            {dg      , 2, 0, [1,]};
            {dg      , 2, 1, [1,]};
            {dg      , 2, 2, [1,]};
            {dg      , 2, 3, [1,]};
            {dg      , 2, 4, [1,]};
            {dg      , 3, 0, [1,]};
            {dg      , 3, 1, [1,]};
            {dg      , 3, 2, [1,]};
            {dg      , 3, 3, [1,]};
            {dg      , 3, 4, [1,]};
            {dg      , 3, 5, [1,]};
            {p0,       1, 0, [1,]};
            {p0,       2, 0, [1,]};
            {p0,       3, 0, [1,]};
        };
    end

    methods (Test)

        function test_nodality(t, factory)
            [fun, dim, degree] = factory{:};
            e = fun(dim, degree);
            action = e.evaluate_dual_basis(@e.tabulate_basis);
            action = reshape(action, [e.fe_space_dim, e.fe_space_dim]);
            t.verifyEqual(action, eye(e.fe_space_dim), 'AbsTol', 1e-12);
        end

        function test_derivatives(t)  %#ok<MANU>
            % FIXME: How to test it?
        end

        function test_facet_dofs(t)  %#ok<MANU>
            % FIXME: How to test it?
        end

        function test_entity_dofs(t)  %#ok<MANU>
            % FIXME: How to test it?
        end

        function test_value_shape(t, factory)
            [fun, dim, degree, expected_value_shape] = factory{:};
            e = fun(dim, degree);
            t.verifyEqual(e.value_shape, expected_value_shape);
        end

    end

end


function fun = lagrange()
    fun = @fe.create_lagrange_element;
end


function fun = dg()
    fun = @fe.create_discontinuous_lagrange_element;
end


function fun = nedelec()
    fun = @fe.create_nedelec_element;
end


function fun = rt()
    fun = @fe.create_raviart_thomas_element;
end


function fun = p0()
    fun = @(dim, degree) fe.create_p0_element(dim);
end
