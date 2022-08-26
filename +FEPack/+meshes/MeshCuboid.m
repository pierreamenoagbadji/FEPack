%> @file MeshCuboid.m
%> @brief Contains the meshes.MeshCuboid class.
% =========================================================================== %
%> @brief Class for periodic meshes of cuboids
%>
%> A mesh here consists of nodes and elements (triangles and edges connecting
%> the nodes) of a cuboid. We only consider periodic meshes
% =========================================================================== %
classdef MeshCuboid < FEPack.meshes.Mesh
  % FEPack.meshes.MeshCuboid < FEPack.meshes.Mesh

  properties

  end

  methods

    % ============= %
    % Create a mesh %
    % ============= %
    function mesh = MeshCuboid(is_structured, BBx, BBy, BBz, numNodesX,...
                               numNodesY, numNodesZ, side_names, name)
      % MeshCuboid constructor for mesh of a cuboid
      %
      % INPUTS:  * is_structured (boolean) indicates if the mesh is
      %            structured or not;
      %          * BBx (2-sized vector) contains the x-coordinates
      %            of the bounding box;
      %          * BBy (2-sized vector) contains the y-coordinates
      %            of the bounding box;
      %          * BBz (2-sized vector) contains the z-coordinates
      %            of the bounding box;
      %          * numNodesX (integer) is the number of nodes on the x-edges;
      %          * numNodesY (integer) is the number of nodes on the y-edges;
      %          * numNodesZ (integer) is the number of nodes on the y-edges;
      %          * side_names (8x1 string) contains the side names (optional)
      %
      % OUTPUTS: * mesh (MeshCuboid), the mesh.

      % Default arguments
      if (nargin < 9), randomName(mesh); end
      if (nargin < 8), side_names = {'xmax'; 'xmin'; 'ymax'; 'ymin'; 'zmax'; 'zmin'}; end
      if (nargin < 7), numNodesZ = 4; end
      if (nargin < 6), numNodesY = 4; end
      if (nargin < 5), numNodesX = 4; end
      if (nargin < 4), BBz = [0.0, 1.0]; end
      if (nargin < 3), BBy = [0.0, 1.0]; end
      if (nargin < 2), BBx = [0.0, 1.0]; end
      if (nargin < 1), is_structured = false; end

      % Generate the .geo file
      fileID = fopen('FEPackmesh.geo', 'w');

      geomsg = ['Include "', FEPack.FEPackObject.pathCpp, '/+FEPack/+tools/FEPackGmsh_macros.geo";\n\n',...
                'h0 = 0.1;\n',...
                'is_structured = ', int2str(is_structured), ';\n',...
                'x1 = ', num2str(BBx(1), '%0.8f'), '; y1 = ', num2str(BBy(1), '%0.8f'), '; z1 = ', num2str(BBz(1), '%0.8f'), ';\n',...
                'x2 = ', num2str(BBx(2), '%0.8f'), '; y2 = ', num2str(BBy(1), '%0.8f'), '; z2 = ', num2str(BBz(1), '%0.8f'), ';\n',...
                'x3 = ', num2str(BBx(2), '%0.8f'), '; y3 = ', num2str(BBy(2), '%0.8f'), '; z3 = ', num2str(BBz(1), '%0.8f'), ';\n',...
                'x4 = ', num2str(BBx(1), '%0.8f'), '; y4 = ', num2str(BBy(2), '%0.8f'), '; z4 = ', num2str(BBz(1), '%0.8f'), ';\n',...
                'x5 = ', num2str(BBx(1), '%0.8f'), '; y5 = ', num2str(BBy(1), '%0.8f'), '; z5 = ', num2str(BBz(2), '%0.8f'), ';\n',...
                'x6 = ', num2str(BBx(2), '%0.8f'), '; y6 = ', num2str(BBy(1), '%0.8f'), '; z6 = ', num2str(BBz(2), '%0.8f'), ';\n',...
                'x7 = ', num2str(BBx(2), '%0.8f'), '; y7 = ', num2str(BBy(2), '%0.8f'), '; z7 = ', num2str(BBz(2), '%0.8f'), ';\n',...
                'x8 = ', num2str(BBx(1), '%0.8f'), '; y8 = ', num2str(BBy(2), '%0.8f'), '; z8 = ', num2str(BBz(2), '%0.8f'), ';\n',...
                'numNodesX = ', int2str(numNodesX), '; numNodesY = ', int2str(numNodesY), '; numNodesZ = ', int2str(numNodesZ), ';\n\n',...
                'domain_name = "rect";\n',...
                'side_name1 = "', side_names{1}, '";\nside_name2 = "', side_names{2}, '";\n',...
                'side_name3 = "', side_names{3}, '";\nside_name4 = "', side_names{4}, '";\n',...
                'side_name5 = "', side_names{5}, '";\nside_name6 = "', side_names{6}, '";\n\n',...
                'Call FEPack_Cuboid;\n\n',...
                'Physical Surface("', side_names{1}, '") = domain_1[];\n',...
                'Physical Surface("', side_names{2}, '") = domain_2[];\n',...
                'Physical Surface("', side_names{3}, '") = domain_3[];\n',...
                'Physical Surface("', side_names{4}, '") = domain_4[];\n',...
                'Physical Surface("', side_names{5}, '") = domain_5[];\n',...
                'Physical Surface("', side_names{6}, '") = domain_6[];\n',...
                'Physical Volume("rect")= domain_7[];\n\n',...
                'Mesh.Format = 50;\n',...
                'Mesh.ElementOrder = 1;\n',...
                'Mesh.MshFileVersion = 2.2;'
               ];

      fprintf(fileID, geomsg);
      fclose(fileID);

      % Use Gmsh
      system([FEPack.FEPackObject.pathBash, '/gmsh-4.10.5-Linux64/bin/gmsh FEPackmesh.geo -3']);
      FEPackmesh;

      % Construct the mesh
      mesh.dimension = 3;
      mesh.numEdgeNodes = [numNodesX; numNodesY; numNodesZ];
      mesh.numPoints = msh.nbNod;
      mesh.points = msh.POS;
      mesh.numTriangles = size(msh.TRIANGLES, 1);
      mesh.triangles = msh.TRIANGLES(:, 1:3);
      mesh.numTetrahedra = size(msh.TETS, 1);
      mesh.tetrahedra = msh.TETS(:, 1:4);
      mesh.refTriangles = msh.TRIANGLES(:, 4);
      mesh.refTetrahedra = sparse(mesh.numTetrahedra, 1);

      if (nargin >= 8)
        mesh.name = name;
      end

      % Construct maps between edge nodes and subdomains
      Icoo = [2, 3; 2, 3; 3, 1; 3, 1; 1, 2; 1, 2];
      for idom = 1:6
        % Maps between edge nodes
        pts = unique(mesh.triangles(mesh.refTriangles == idom, :));
        [~, cle] = sortrows(mesh.points(pts, :), Icoo(idom, :));
        mesh.maps{idom} = pts(cle);

        % Subdomains
        mesh.domains{idom} = FEPack.meshes.FEDomain(mesh, side_names{idom}, 2, idom, pts(cle));
      end
      mesh.domains{7} = FEPack.meshes.FEDomain(mesh, 'volumic', 3, 0);

    end

    % function childmesh = toDomain(mesh)
    %
    % end



  end

end
