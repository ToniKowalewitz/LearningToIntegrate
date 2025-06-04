function mesh = generate_half_disc_mesh(n)
    % Build half-disc mesh
    %
    % SYNTAX
    %   mesh = generate_half_disc_mesh(n)
    %
    % INPUT/OUTPUT PARAMETERS
    %
    %   n    ... Number of divisions
    %   mesh ... Instance of Mesh class

    num_vertices = (3*n + 2)*(n + 1)/2;
    num_cells = 3*n*n;

    vertex_coords = zeros(2, num_vertices);
    cells = zeros(3, num_cells, 'uint32');

    % Vertex at origin
    vertex_coords(:, 1) = [0; 0];

    % Remaining vertices
    v = 2;
    for i = 1:n
        for j = 0:3*i
            r = i/n;
            th = pi*j/(3*i);
            x = r*cos(th);
            y = r*sin(th);
            vertex_coords(:, v) = [x; y];
            v = v + 1;
        end
    end
    assert(v == num_vertices + 1);

    % Build cells
    c = 1;
    for i = 1:n
        for k = 0:2
            for j = 0:(i*2 - 2)
                if mod(j, 2) == 0
                    i1 = 1 + (3*i - 1)*i/2     + k*i     + floor(j/2)    ;
                    i2 = 1 + (3*i - 1)*i/2     + k*i     + floor(j/2) + 1;
                    i0 = 1 + (3*i - 4)*(i-1)/2 + k*(i-1) + floor(j/2)    ;
                else
                    i0 = 1 + (3*i - 4)*(i-1)/2 + k*(i-1) + floor(j/2)    ;
                    i1 = 1 + (3*i - 4)*(i-1)/2 + k*(i-1) + floor(j/2) + 1;
                    i2 = 1 + (3*i - 1)*i/2     + k*i     + floor(j/2) + 1;
                end
                cells(:, c) = [i0; i1; i2];
                c = c + 1;
            end
        end
    end
    assert(c == num_cells + 1);

    mesh = meshing.Mesh(2, vertex_coords, cells);
end
