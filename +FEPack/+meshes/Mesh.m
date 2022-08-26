%> @file Mesh.m
%> @brief Contains the meshes.Mesh class.
% ======================================================================
%> @brief Abstract class for all meshes
%>
%> Each mesh consists of nodes (points in dimension 1, 2, or 3) which
%> are connected by elements.
%>
%> The abstract Mesh class provides the interface for this
%> functionality by storing the dimension of the nodes, as well as their
%> coordinates.
%>
%> Specific informations regarding the elements can be found in the
%> subclasses.
% ======================================================================
classdef Mesh < FEPack.FEPackObject
    % FEPack.meshes.Mesh < FEPack.FEPackObject

    properties (SetAccess = protected)

        %> @brief Total number of nodes in the mesh
        numNodes = 0;

        %> @brief Dimension of the geometry
        dimen = 0;

        %> @brief Matrix that contains the coordinates of the nodes. The k-th
        %> node is represented by the k-th line
        nodes = [];

    end

    methods (Abstract)

    end

end
