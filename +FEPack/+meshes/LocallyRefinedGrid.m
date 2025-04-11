%> @file LocallyRefinedGrid.m
%> @brief Contains the meshes.LocallyRefinedGrid class.
% =========================================================================== %
%> @brief Class for rectangular grids with possible local refinement
%>
%> A mesh here consists of nodes and elements (rectangle connecting
%> the nodes) of a rectangle.
% =========================================================================== %
classdef LocallyRefinedGrid < FEPack.FEPackObject
  % FEPack.meshes.LocallyRefinedGrid < FEPack.meshes.Mesh

  properties (SetAccess = protected)
    
    %> @brief Dimension of the geometry
    dimension = 0;

    %> @brief Total number of vertices in the mesh
    numPoints = 0;

    %> @brief Matrix that contains the coordinates of the points. The k-th
    %> point is represented by the k-th line
    points = [];

    %> @brief Total number of elements in the mesh
    numElements = 0;

    %> @brief Matrix that contains the indices of the elements. The k-th
    %> point is represented by the k-th element. Each line contains the
    %> indices of its vertices
    elements = [];

  end

  methods

    % ============= %
    % Create a mesh %
    % ============= %
    function mesh = LocallyRefinedGrid(dimension, axisLimits, numNodes, name)
      % LocallyRefinedGrid constructor for uniform grid
      %
      % INPUTS:  * dimension (integer) is the dimension, only 1 or 2;
      %          * axisLimits (dimension x 2 matrix) contains the bounds in each direction;
      %          * numNodes (dimension-sized vector) contains the numbers of nodes in each direction;
      %
      % OUTPUTS: * mesh (LocallyRefinedGrid), the mesh.

      % Default arguments
      if (nargin < 4), randomName(mesh); end

      mesh.dimension = dimension;
      
      if (dimension == 1)
        mesh.numPoints = numNodes(1);
        mesh.points = linspace(axisLimits(1, 1), axisLimits(1, 2), numNodes(1));
        mesh.numElements = numNodes(1)-1;
        mesh.elements = [(1:numNodes(1)-1)', (2:numNodes(1))'];
      elseif (dimension == 2)
        Nx = numNodes(1);
        Ny = numNodes(2);
        mesh.numPoints = Nx * Ny;

        X = linspace(axisLimits(1, 1), axisLimits(1, 2), Nx);
        Y = linspace(axisLimits(2, 1), axisLimits(2, 2), Ny);
        mesh.points = zeros(mesh.numPoints, 2);
        for idI = 1:Nx
          for idJ = 1:Ny
            idP = sub2ind([Nx, Ny], idI, idJ);
            mesh.points(idP, :) = [X(idI), Y(idJ)];
          end
        end

        mesh.numElements = (Nx - 1) * (Ny - 1);
        mesh.elements = zeros(mesh.numElements, 4);
        for idCol = 1:Nx-1
          for idRow = 1:Ny-1
            idP1 = sub2ind([Nx, Ny], idCol  , idRow);
            idP2 = sub2ind([Nx, Ny], idCol+1, idRow);
            idP3 = sub2ind([Nx, Ny], idCol+1, idRow+1);
            idP4 = sub2ind([Nx, Ny], idCol  , idRow+1);

            idE = sub2ind([Nx-1, Ny-1], idCol, idRow);

            mesh.elements(idE, :) = [idP1, idP2, idP3, idP4];
          end
        end
      else
        error(['dimension must be either 1 or 2. The value ', int2str(dimension), ' is not allowed.']);
      end
      
      if (nargin >= 4)
        mesh.name = name;
      end
    end

    function visualize(mesh, show_numbering)

      if (nargin < 2)
        show_numbering = false;
      end

      figure;
      set(groot,'defaultAxesTickLabelInterpreter','latex');
      set(groot,'defaulttextinterpreter','latex');
      set(groot,'defaultLegendInterpreter','latex');

      if (mesh.dimension == 1)

        plot([mesh.points(1), mesh.points(end)], [0 0], 'k');
        hold on;
        plot(mesh.points, 0, 'ro', 'MarkerSize', 8);
        
        if (show_numbering)
          for idP = 1:mesh.numPoints
            text(mesh.points(idP), 0.03, int2str(idP), 'Color', 'red', 'FontSize', 14);
          end
          
          for idE = 1:mesh.numElements
            P1 = mesh.points(mesh.elements(idE, 1));
            P2 = mesh.points(mesh.elements(idE, 2));
            text(0.5*(P1 + P2), -0.03, int2str(idE), 'Color', 'black', 'FontSize', 14);
          end
        end
        ylim([-0.2, 0.2]);

      elseif (mesh.dimension == 2)
        for idE = 1:mesh.numElements
          P1 = mesh.points(mesh.elements(idE, 1), :);
          P2 = mesh.points(mesh.elements(idE, 2), :);
          P3 = mesh.points(mesh.elements(idE, 3), :);
          P4 = mesh.points(mesh.elements(idE, 4), :);

          plot([P1(1) P2(1) P3(1) P4(1)], [P1(2) P2(2) P3(2) P4(2)], 'k');
          hold on;
        end

        plot(mesh.points(:, 1), mesh.points(:, 2), 'ro', 'MarkerSize', 8);

        if (show_numbering)
          for idP = 1:mesh.numPoints
            text(mesh.points(idP, 1), mesh.points(idP,2)+0.02, int2str(idP), 'Color', 'red', 'FontSize', 14);
          end

          for idE = 1:mesh.numElements
            P1 = mesh.points(mesh.elements(idE, 1), :);
            P3 = mesh.points(mesh.elements(idE, 3), :);

            text(0.5*(P1(1) + P3(1)), 0.5*(P1(2) + P3(2)), int2str(idE), 'Color', 'black', 'FontSize', 14);
          end
        end

        % if 
      else
        error(['dimension must be either 1 or 2. The value ', int2str(mesh.dimension), ' is not allowed.']);
      end
      set(gca, 'DataAspectRatio', [1 1 1]);

    end

    % function childmesh = toDomain(mesh)
    %
    % end
  end

end
