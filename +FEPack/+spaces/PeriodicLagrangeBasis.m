%> @file PeriodicLagrangeBasis.m
%> @brief Contains the pdes.PeriodicLagrangeBasis class.
% =========================================================================== %
%> @brief Abstract class for all domains
%>
% =========================================================================== %
classdef PeriodicLagrangeBasis < FEPack.spaces.SpectralBasis
  % FEPack.spaces.PeriodicLagrangeBasis < FEPack.spaces.SpectralBasis

  properties (SetAccess = protected)

  end

  methods

    function sp = PeriodicLagrangeBasis(domain)
      % PERIODICLAGRANGEBASIS

      if ~(isa(domain.mesh, 'FEPack.meshes.MeshSegment') ||...
           isa(domain.mesh, 'FEPack.meshes.MeshRectangle') ||...
           isa(domain.mesh, 'FEPack.meshes.MeshCuboid'))
           error('Je ne peux définir des fonctions de base de Lagrange périodisées que sur un segment, un rectangle, ou un cube.');
      end

      % Compute periodic Lagrange basis
      % ///////////////////////////////
      u = FEPack.pdes.PDEObject;
      ecs = FEPack.pdes.EssentialConditions;

      for idI = 1:domain.mesh.dimension
        % Periodic conditions on each face
        ecs = ecs & assignEcs((u|domain.mesh.domains{2*idI-1}) - (u|domain.mesh.domains{2*idI}), 0.0);
      end

      ecs.applyEcs;    % Compute the projection matrix
      Q = ecs.P(:, domain.IdPoints); % Restrict the projection matrix to the domain

      % Find and eliminate zero rows in the restricted projection matrix
      [Innz, ~, ~] = find(Q);
      zerosRows = (1:size(Q, 1)); zerosRows(Innz) = [];
      Q(zerosRows, :) = [];

      % Construct the basis functions
      [idI, idJ, ~] = find(Q);
      phis = sparse(idI, idJ, 1).';

      % Set parameters
      % //////////////
      sp@FEPack.spaces.SpectralBasis(domain, phis, size(phis, 2));
      sp.is_interpolated = 1;
      sp.computeBasisMatrices(0);

    end

  end

end
