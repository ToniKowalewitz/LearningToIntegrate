function mesh = generate_unit_cube_mesh(divisions)
    % Build unit cube mesh
    %
    % SYNTAX
    %   mesh = generate_unit_cube_mesh(divisions)
    %
    % INPUT/OUTPUT PARAMETERS
    %
    %   divisions ... Number of divisions along axes, i.e., [nx, ny, ...]
    %   mesh      ... Instance of Mesh class
    %
    % REMARKS
    %
    %   generate_unit_cube_mesh([3, 5]) results in unit square mesh of
    %   triangles with 3 divisiones along x and 5 along y axis. Analogically,
    %   generate_unit_cube_mesh([3, 4, 5]) gives mesh of tetrahedrons.

    if numel(divisions) == 2
      mesh = generate_unit_2cube_mesh(divisions(1), divisions(2));
    elseif numel(divisions) == 3
      mesh = generate_unit_3cube_mesh(divisions(1), divisions(2), divisions(3));
    else
      error('Dimensions %d not implemented', numel(divisions));
    end
end


function mesh = generate_unit_2cube_mesh(nx, ny)
    % Build vertex coordinates
    x = linspace(0, 1, nx + 1);
    y = linspace(0, 1, ny + 1);
    [x, y] = ndgrid(x, y);
    num_vertices = (nx + 1)*(ny + 1);
    vertex_coords = zeros(2, num_vertices);
    vertex_coords(1, :) = reshape(x, num_vertices, 1);
    vertex_coords(2, :) = reshape(y, num_vertices, 1);

    % Build cell-to-vertex topology
    cells = zeros(3, 2*nx*ny, 'uint32');
    for iy = 0:ny-1
        for ix = 0:nx-1
            v0 = 1 + iy*(nx + 1) + ix;
            v1 = v0 + 1;
            v2 = v0 + (nx + 1);
            v3 = v1 + (nx + 1);

            c0 = 1 + 2*(iy*nx + ix);
            cells(:, c0+0) = [v0, v1, v2];
            cells(:, c0+1) = [v1, v2, v3];
        end
    end

    mesh = meshing.Mesh(2, vertex_coords, cells);
end


function mesh = generate_unit_3cube_mesh(nx, ny, nz)
    % Build vertex coordinates
    x = linspace(0, 1, nx + 1);
    y = linspace(0, 1, ny + 1);
    z = linspace(0, 1, nz + 1);
    [x, y, z] = ndgrid(x, y, z);
    num_vertices = (nx + 1)*(ny + 1)*(nz + 1);
    vertex_coords = zeros(3, num_vertices);
    vertex_coords(1, :) = reshape(x, num_vertices, 1);
    vertex_coords(2, :) = reshape(y, num_vertices, 1);
    vertex_coords(3, :) = reshape(z, num_vertices, 1);

    % Build cell-to-vertex topology
    cells = zeros(4, 6*nx*ny*nz, 'uint32');
    for iz = 0:nz-1
        for iy = 0:ny-1
            for ix = 0:nx-1
                v0 = 1 + iz*(nx + 1)*(ny + 1) + iy*(nx + 1) + ix;
                v1 = v0 + 1;
                v2 = v0 + (nx + 1);
                v3 = v1 + (nx + 1);
                v4 = v0 + (nx + 1)*(ny + 1);
                v5 = v1 + (nx + 1)*(ny + 1);
                v6 = v2 + (nx + 1)*(ny + 1);
                v7 = v3 + (nx + 1)*(ny + 1);

                c0 = 1 + 6*(iz*nx*ny + iy*nx + ix);
                cells(:, c0+0) = [v0, v1, v3, v7];
                cells(:, c0+1) = [v0, v1, v5, v7];
                cells(:, c0+2) = [v0, v4, v5, v7];
                cells(:, c0+3) = [v0, v2, v3, v7];
                cells(:, c0+4) = [v0, v4, v6, v7];
                cells(:, c0+5) = [v0, v2, v6, v7];
            end
        end
    end

    mesh = meshing.Mesh(3, vertex_coords, cells);
end
