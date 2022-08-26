%> @file FEDomain.m
%> @brief Contains the meshes.FEDomain class.
% =========================================================================== %
%> @brief Abstract class for all domains
%>
% =========================================================================== %
classdef FEDomain < FEPack.FEPackObject
  % FEPack.meshes.Mesh < FEPack.FEPackObject

  properties (SetAccess = protected)

    %> @brief Dimension of the domain
    dimension = 0;

    %> @brief Mesh
    mesh = [];

    %> @brief Reference of the domain
    reference = 0;

    %> @brief Number of elements within the domain
    numElts = 0;

    %> @brief indices of elements within the domain
    idelements = [];

    %> @brief elements
    elements = [];

    %> @brief Number of points within the domain
    numPoints = 0;

    %> @brief Indices of the points composing the domain
    points = [];

    % %> @brief Indices of the segments composing the domain
    % segments = [];
    %
    % %> @brief Indices of the triangles composing the domain
    % triangles = [];

  end

  methods

    function FEdom = FEDomain(mesh, domain_name, dimension, reference, points)

      FEdom.name = domain_name;
      FEdom.dimension = dimension;
      FEdom.mesh = mesh;
      FEdom.reference = reference;

      switch (dimension)
      case 0
        FEdom.idelements = find(mesh.refPoints == reference);
        FEdom.elements = FEdom.idelements;
      case 1
        FEdom.idelements = find(mesh.refSegments == reference);
        FEdom.elements = mesh.segments(FEdom.idelements, :);
      case 2
        FEdom.idelements = find(mesh.refTriangles == reference);
        FEdom.elements = mesh.triangles(FEdom.idelements, :);
      case 3
        FEdom.idelements = find(mesh.refTetrahedra == reference);
        FEdom.elements = mesh.tetrahedra(FEdom.idelements, :);
      end
      FEdom.numElts = length(FEdom.elements);

      if (nargin < 5)
        points = unique(FEdom.elements);
      end
      FEdom.points = points;
      FEdom.numPoints = length(FEdom.points);
      % FEdom.points = find(mesh.refPoints == reference);
      % FEdom.segments = find(mesh.refSegments == reference);
      % FEdom.triangles = find(mesh.refTriangles == reference);

    end

  end

end
