classdef TestTransferMatrix < matlab.unittest.TestCase

    properties (TestParameter)
        factory = {
            {2, nedelec, 1, dg, 1};
            {2, nedelec, 2, dg, 2};
            {3, nedelec, 1, dg, 1};
            {3, nedelec, 2, dg, 2};
        };
    end

    methods (Test)

        function test_transfer_matrix(t, factory)
            % Test function create_transfer_matrix()
            %
            % This rather a test showing how to create transfer matrix
            % as it is not testing a public function from the package

            [dim, ofun, odegree, dfun, ddegree] = factory{:};

            o = ofun(dim, odegree);
            d = dfun(dim, ddegree);

            M = create_transfer_matrix(o, d);

            for i = 1:o.fe_space_dim
                % i-the basis function in origin element
                fi = @(x) util.paren(o.tabulate_basis(x), i, 1:dim);

                for j = 1:d.fe_space_dim
                    for k = 1:dim
                        % k-th component of j-th dual basis member in destination element
                        Fjk = @(v) util.paren(d.evaluate_dual_basis(v), k, j);

                        % Check that a component in M is matching
                        t.verifyEqual(M(i, k, j), Fjk(fi), 'AbsTol', 1e-15);
                    end
                end
            end

        end

    end

end


function M = create_transfer_matrix(origin_element, destination_element)

    o = origin_element;
    d = destination_element;

    dim = o.simplex.dim;
    assert(d.simplex.dim == dim);

    % Value size
    vs = dim;  % NB: Now assuming vector-valued elements

    basis_flat = @(x) reshape(o.tabulate_basis(x), o.fe_space_dim*vs, 1);
    M = d.evaluate_dual_basis(basis_flat);
    M = reshape(M, o.fe_space_dim, vs, d.fe_space_dim);

end


function fun = dg()
    fun = @fe.create_discontinuous_lagrange_element;
end


function fun = nedelec()
    fun = @fe.create_nedelec_element;
end
