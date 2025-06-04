function lat = barycentric_lattice(dim, q)
    % Create simplex lattice in Z^dim; simplex is
    % a plane intersecting a cube in barycentric coords.

    cube = cartesian(0:q, dim+1);
    simplex = cube(sum(cube, 2)==q, :);
    lat = simplex;
end


function C = cartesian(seq, n)
    % Create Cartesian power seq^n, i.e., cartesian('abc', 2)
    % gives 'aa' 'ab' 'ac' 'ba' 'bb' 'bc' 'ca' 'cb' 'cc'

    args = cell(1, n);
    args(:) = {seq};
    [F{1:n}] = ndgrid(args{:});
    for i=n:-1:1
        G(:,i) = F{i}(:);
    end
    C = unique(G , 'rows');
end
