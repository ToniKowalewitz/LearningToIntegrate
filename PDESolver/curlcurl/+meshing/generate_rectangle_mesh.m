function mesh = generate_rectangle_mesh(xmin, xmax, ymin, ymax, varargin)
    % Build unit cube mesh using Matlab's PDE toolbox
    %
    % SYNTAX
    %   mesh = generate_rectangle_mesh(xmin, xmax, ymin, ymax, varargin)
    %
    % INPUT/OUTPUT PARAMETERS
    %
    %   xmin     ... Minimal x-coordinate
    %   xmax     ... Maximal x-coordinate
    %   ymin     ... Minimal y-coordinate
    %   ymax     ... Maximal y-coordinate
    %   varargin ... Any additional arguments to Matlab's 'generateMesh()'
    %   mesh     ... Instance of Mesh class
    %
    % EXAMPLES
    %
    %   mesh = meshing.generate_rectangle_mesh(0, 2, 0, 4);
    %   figure();
    %   meshing.plot_mesh(mesh);
    %
    %   mesh = meshing.generate_rectangle_mesh(0, 2, 0, 4, 'Hmax', 0.1);
    %   figure();
    %   meshing.plot_mesh(mesh);

    model = createpde(1);
    rectangle = [3; 4; xmin; xmax; xmax; xmin; ymax; ymax; ymin; ymin];
    g = decsg(rectangle);
    geometryFromEdges(model, g);
    generateMesh(model, 'GeometricOrder', 'linear', varargin{:});
    mesh = meshing.Mesh(2, model.Mesh.Nodes, uint32(model.Mesh.Elements));
end
