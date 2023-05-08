%> @file MeshSegment.m
%> @brief Contains the meshes.MeshSegment class.
% =========================================================================== %
%> @brief Class for one-dimensional meshes
%>
%> A mesh here consists of nodes and elements (segments connecting the nodes)
%> of a bounded interval.
% =========================================================================== %
classdef MeshSegment < FEPack.meshes.Mesh
  % FEPack.meshes.MeshSegment < FEPack.meshes.Mesh

  properties

  end

  methods

    % ============= %
    % Create a mesh %
    % ============= %
    function mesh = MeshSegment(varargin)
      % MeshSegment constructor for segment mesh
      % The first argument (string) indicates the type of mesh:
      %     - for a uniform mesh, the syntax is
      %         mesh = meshes.MeshSegment('uniform', a, b, N, FEorder, side_names, name)
      %       a, b are the bounds of the segment and N the number of nodes
      %
      %     - for a mesh generated from points, the syntax is
      %         mesh = meshes.MeshSegment('vertices', points, FEorder, side_names, name)
      %       vertices is a list of the nodes
      %
      %     FEorder (integer) is the order of FE;
      %     name (optional) is the name of the mesh
      if (nargin < 1)

        % Default argmuments
        numPoints = 8;
        points = linspace(0, 1, numPoints).';
        randomName(mesh);
        side_names = {'xmin'; 'xmax'};

      elseif strcmpi(varargin{1}, 'uniform')

        % Generate a mesh of equispaced points
        numPoints = varargin{4};
        points = sort(linspace(varargin{2}, varargin{3}, numPoints).');

        % Give a random name if nothing is specified
        if (nargin < 7)
          randomName(mesh);
        end

        if (nargin < 6)
          side_names = {'xmin'; 'xmax'};
        else
          side_names = varargin{6};
        end

        % Default FE order
        if (nargin < 5)
          FEorder = 1;
        end

      elseif strcmpi(varargin{1}, 'vertices')

        % Generate a mesh from a list of vertices
        points = unique(varargin{2});
        numPoints = length(points);

        % Give a random name if nothing is specified
        if (nargin < 5)
          randomName(mesh);
        end

        if (nargin < 4)
          side_names = {'xmin'; 'xmax'};
        else
          side_names = varargin{4};
        end

        % Default FE order
        if (nargin < 3)
          FEorder = 1;
        end

      else

        error(['Type de maillage ''', varargin{1}, ''' non reconnu. ',...
               'Les options possibles sont ''uniform'' et ''vertices''.']);

      end

      % Properties of the mesh
      mesh.dimension = 1;
      mesh.numEdgeNodes = [1; 0; 0];
      mesh.numPoints = numPoints;
      mesh.points = [points(:), zeros(numPoints, 2)];
      mesh.numSegments = numPoints - 1;
      mesh.segments = [(1:numPoints-1).', (2:numPoints).'];
      mesh.refPoints = [1; zeros(numPoints-2, 1); 2];
      mesh.refSegments = zeros(mesh.numSegments, 1);
      mesh.attachFEorder(FEorder);
      mesh.maps{1} = numPoints;
      mesh.maps{2} = 1;
      numD = [2; 1]; % The domains are ordered as : xmax - xmin
      for idom = 1:2
        mesh.domains{idom} = FEPack.meshes.FEDomain(mesh, side_names{numD(idom)}, 0, idom, mesh.maps{idom});
      end
      mesh.domains{3} = FEPack.meshes.FEDomain(mesh, 'volumic', 1, 0);
      mesh.mapdomains = [mesh.domains{1}.reference, mesh.domains{2}.reference];

    end

  end

end
