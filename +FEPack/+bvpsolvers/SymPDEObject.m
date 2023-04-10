%> @file SymPDEObject.m
%> @brief Contains the pdes.SymPDEObject class.
% =========================================================================== %
%> @brief class used for symbolic boundary conditions
% =========================================================================== %
classdef SymPDEObject < FEPack.FEPackObject
  % FEPack.pdes.SymPDEObject < FEPack.FEPackObject

  properties (SetAccess = protected)

    %> @brief Indicates if normal dierivative or not
    is_normal_derivative = 0;

  end

  methods

    % Normal derivative
    function dudn = dn(u)
      dudn = copy(u);
      dudn.is_normal_derivative = 1;
    end

    % For boundary conditions
    function ecs = onDomain(u, domain)
      if isa(domain, 'FEPack.meshes.FEDomain')

        % Symbolic boundary conditions
        ecs = FEPack.bvpsolvers.SymBoundaryCondition;
        ecs.domains = {domain};
        [~, domId] = find(domain.mesh.mapdomains == domain.reference);
        alpha = [0 0];
        alpha(domId) = 1;
        
        if (u.is_normal_derivative)
          ecs.gamma0 = {[0 0]};
          ecs.gamma1 = {alpha};
        else
          ecs.gamma0 = {alpha};
          ecs.gamma1 = {[0 0]};
        end

        % IdPoints = domain.IdPoints;
        % m = size(IdPoints, 1);

        % ecs.C = sparse(1:m, IdPoints, 1, m, domain.mesh.numPoints);

      else

        error(['L''operateur | n''est pas compatible avec le terme ', ...
               'de droite choisi (qui est de type', class(domain), ').']);

      end % if
    end

    function ecs = or(u, domain)
      ecs = onDomain(u, domain);
    end

    % % For user-defined symbolic boundary conditions
    % function bcs = onxmin(u)
    %   bcs = FEPack.pdes.SymBoundaryCondition;
    %   if (u.is_normal_derivative)
    %     bcs.dudnBCxmin = {1};
    %   else
    %     bcs.uBCxmin = {1};
    %   end
    % end

  end % methods

end
