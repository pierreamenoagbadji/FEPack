%> @file EssentialBCLagrangeBasis.m
%> @brief Contains the pdes.EssentialBCLagrangeBasis class.
% =========================================================================== %
%> @brief Abstract class for all domains
%>
% =========================================================================== %
classdef EssentialBCLagrangeBasis < FEPack.spaces.SpectralBasis
  % FEPack.spaces.EssentialBCLagrangeBasis < FEPack.spaces.SpectralBasis

  properties (SetAccess = protected)

  end

  methods

    function sp = EssentialBCLagrangeBasis(domain, ecs_fun)
      % ESSENTIALBCLAGRANGEBASIS
      %
      % Lagrange basis associated to a particular essential condition
      % 
      % INPUTS: domain, the domain on which the basis is defined;
      %         ecs_fun, function handle of the form ecs_fun = @(u, dom) ...
      %                  returns the essential boundary conditions.
      %                  For example,
      % Dirichlet conditions: @(u, dom) (u|dom.mesh.domains{idI} == 0) for any idI
      % Periodic conditions:  @(u, dom) (u|dom.mesh.domains{2*idI-1}) - (u|dom.mesh.domains{2*idI}) == 0 for any idI

      if ~(isa(domain.mesh, 'FEPack.meshes.MeshSegment') ||...
           isa(domain.mesh, 'FEPack.meshes.MeshRectangle') ||...
           isa(domain.mesh, 'FEPack.meshes.MeshCuboid'))
           error('Je ne peux définir des fonctions de base de Lagrange périodisées que sur un segment, un rectangle, ou un cube.');
      end

      % Compute periodic Lagrange basis
      % ///////////////////////////////
      u = FEPack.pdes.PDEObject;
      ecs = ecs_fun(u, domain);
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
