%> @file MeshRectangle.m
%> @brief Contains the meshes.MeshRectangle class.
% =========================================================================== %
%> @brief Class for periodic meshes of rectangles
%>
%> A mesh here consists of nodes and elements (triangles and edges connecting
%> the nodes) of a rectangle. We only consider periodic meshes
% =========================================================================== %
classdef MeshRectangle < FEPack.meshes.Mesh
  % FEPack.meshes.MeshRectangle < FEPack.meshes.Mesh

  properties

  end

  methods

    % ============= %
    % Create a mesh %
    % ============= %
    function mesh = MeshRectangle(is_structured, BBx, BBy, numNodesX,...
                                  numNodesY, side_names, name)
      % MeshRectangle constructor for rectangle mesh
      %
      % INPUTS:  * is_structured (boolean) indicates if the mesh is
      %            structured or not;
      %          * BBx (2-sized vector) contains the x-coordinates
      %            of the bounding box;
      %          * BBy (2-sized vector) contains the y-coordinates
      %            of the bounding box;
      %          * numNodesX (integer) is the number of nodes on the x-edges;
      %          * numNodesY (integer) is the number of nodes on the y-edges;
      %          * side_names (4x1 string) contains the side names (optional)
      %
      % OUTPUTS: * mesh (MeshRectangle), the mesh.

      % Default arguments
      if (nargin < 7), randomName(mesh); end
      if (nargin < 6), side_names = {'ymin'; 'xmax'; 'ymax'; 'xmin'}; end
      if (nargin < 5), numNodesY = 8; end
      if (nargin < 4), numNodesX = 8; end
      if (nargin < 3), BBy = [0.0, 1.0]; end
      if (nargin < 2), BBx = [0.0, 1.0]; end
      if (nargin < 1), is_structured = false; end

      % Generate the .geo file
      fileID = fopen('FEPackmesh.geo', 'w');

      geomsg = ['Include "', FEPack.FEPackObject.pathCpp, '/+FEPack/+tools/FEPackGmsh_macros.geo";\n\n',...
                'h0 = 0.1;\n',...
                'is_structured = ', int2str(is_structured), ';\n',...
                'x1 = ', num2str(BBx(1), '%0.8f'), '; y1 = ', num2str(BBy(1), '%0.8f'), '; z1 = 0.0;\n',...
                'x2 = ', num2str(BBx(2), '%0.8f'), '; y2 = ', num2str(BBy(1), '%0.8f'), '; z2 = 0.0;\n',...
                'x3 = ', num2str(BBx(2), '%0.8f'), '; y3 = ', num2str(BBy(2), '%0.8f'), '; z3 = 0.0;\n',...
                'x4 = ', num2str(BBx(1), '%0.8f'), '; y4 = ', num2str(BBy(2), '%0.8f'), '; z4 = 0.0;\n',...
                'numNodesX = ', int2str(numNodesX), '; numNodesY = ', int2str(numNodesY), ';\n\n',...
                'domain_name = "rect";\n',...
                'side_name1 = "', side_names{1}, '";\nside_name2 = "', side_names{2}, '";\n',...
                'side_name3 = "', side_names{3}, '";\nside_name4 = "', side_names{4}, '";\n\n',...
                'Call FEPack_Rectangle;\n\n',...
                'Physical Point(1) = domain_0[];\n',...
                'Physical Line("', side_names{1}, '") = domain_1[];\n',...
                'Physical Line("', side_names{2}, '") = domain_2[];\n',...
                'Physical Line("', side_names{3}, '") = domain_3[];\n',...
                'Physical Line("', side_names{4}, '") = domain_4[];\n',...
                'Physical Surface("rect")= domain_5[];\n\n',...
                'Mesh.Format = 50;\n',...
                'Mesh.ElementOrder = 1;\n',...
                'Mesh.MshFileVersion = 2.2;'
               ];

      fprintf(fileID, geomsg);
      fclose(fileID);

      % Use Gmsh
      system([FEPack.FEPackObject.pathBash, '/gmsh-4.10.5-Linux64/bin/gmsh FEPackmesh.geo -2']);
      FEPackmesh;

      % Construct the mesh
      mesh.dimension = 2;
      mesh.numEdgeNodes = [numNodesX; numNodesY; 0];
      mesh.numPoints = msh.nbNod;
      mesh.points = msh.POS;
      mesh.numSegments = size(msh.LINES, 1);
      mesh.segments = msh.LINES(:, 1:2);
      mesh.numTriangles = size(msh.TRIANGLES, 1);
      mesh.triangles = msh.TRIANGLES(:, 1:3);
      mesh.refSegments = msh.LINES(:, 3);
      mesh.refTriangles = sparse(mesh.numTriangles, 1);

      if (nargin >= 8)
        mesh.name = name;
      end

      % Construct maps between edge nodes and subdomains
      % The domains are ordered as : xmax - xmin - ymax - ymin
      evec = [2; 3; 4; 5];
      Icoo = [1; 2; 1; 2];
      numD = [4; 1; 3; 2];

      for idom = 1:4
        % Maps between edge nodes
        pts = unique(mesh.segments(mesh.refSegments == evec(idom), :));
        [~, cle] = sort(mesh.points(pts, Icoo(idom)));
        mesh.maps{numD(idom)} = pts(cle);

        % Subdomains
        mesh.domains{numD(idom)} = FEPack.meshes.FEDomain(mesh, side_names{idom}, 1, evec(idom), pts(cle));
      end
      mesh.domains{5} = FEPack.meshes.FEDomain(mesh, 'volumic', 2, 0);

    end

    % function childmesh = toDomain(mesh)
    %
    % end



  end

end
