function cols = select_linearly_independent_measurements_impl(Mtx, Mrx)
    % Select linearly independent measurements
    %
    % This is a lower-level function. Normally, user would use
    %
    %   [Mtx, Mrx] = select_linearly_independent_measurements(Mtx, Mrx);
    %
    % For matrices Mtx, Mrx, which define measurements, select
    % their columns which make the measurements linearly independent.

    M = to_mat(Mtx, Mrx);
    if issparse(M)
        cols = select_subset_svds(M);
    else
        cols = select_subset_svd(M);
    end
    cols = sort(cols);
end


function M = to_mat(M1, M2)
    % Computes the product
    %
    %   M_{mk} = M1_{ik} * M2_{jk},
    %
    % where m is given "lexicographically" by (i,j),
    % and there is no summation over k.

    assert(size(M1, 2) == size(M2, 2));
    [ni, ~ ] = size(M1);
    [nj, nk] = size(M2);
    M = zeros(ni*nj, nk, 'like', M1(1)*M2(1));
    for k = 1:nk
        M(:, k) = kron(M1(:, k), M2(:, k));
    end
end


function cols = select_subset_svd(A, tol)
    % Select columns of A such that selected submatrix has full rank
    %
    % See Algorithm 12.2.1 in
    %
    %   @book {MR1417720,
    %       AUTHOR = {Golub, Gene H. and Van Loan, Charles F.},
    %        TITLE = {Matrix computations},
    %       SERIES = {Johns Hopkins Studies in the Mathematical Sciences},
    %      EDITION = {Third},
    %    PUBLISHER = {Johns Hopkins University Press, Baltimore, MD},
    %         YEAR = {1996},
    %        PAGES = {xxx+698},
    %         ISBN = {0-8018-5413-X; 0-8018-5414-8},
    %      MRCLASS = {65-02 (65Fxx)},
    %     MRNUMBER = {1417720},
    %          URL = {https://web.mit.edu/ehliu/Public/sclark/Golub%20G.H.,%20Van%20Loan%20C.F.-%20Matrix%20Computations.pdf},
    %   }

    % Compute SVD
    [~, S, V] = svd(full(A), 'econ');

    % Compute rank
    s = diag(S);
    clear S;
    if nargin==1
        tol = max(size(A)) * eps(max(s));
    end
    rnk = sum(s > tol);
    clear s;

    % Compute QR column permutation
    [~, ~, p] = qr(V(:, 1:rnk)', 0);
    clear V;

    % Return linearly independent columns
    cols = p(1:rnk);
end


function cols = select_subset_svds(A)
    % Select columns of A such that selected submatrix has full rank
    %
    % See Algorithm 12.2.1 in
    %
    %   @book {MR1417720,
    %       AUTHOR = {Golub, Gene H. and Van Loan, Charles F.},
    %        TITLE = {Matrix computations},
    %       SERIES = {Johns Hopkins Studies in the Mathematical Sciences},
    %      EDITION = {Third},
    %    PUBLISHER = {Johns Hopkins University Press, Baltimore, MD},
    %         YEAR = {1996},
    %        PAGES = {xxx+698},
    %         ISBN = {0-8018-5413-X; 0-8018-5414-8},
    %      MRCLASS = {65-02 (65Fxx)},
    %     MRNUMBER = {1417720},
    %          URL = {https://web.mit.edu/ehliu/Public/sclark/Golub%20G.H.,%20Van%20Loan%20C.F.-%20Matrix%20Computations.pdf},
    %   }

    % Compute SVD
    [~, ~, V] = svds(A, min(size(A)), 'smallestnz');

    % Compute rank
    rnk = size(V, 2);

    % Compute QR column permutation
    [~, ~, p] = qr(V(:, 1:rnk)', 0);
    clear V;

    % Return linearly independent columns
    cols = p(1:rnk);
end
