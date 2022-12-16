%> @file PDEObject.m
%> @brief Contains the pdes.PDEObject class.
% =========================================================================== %
%> @brief class for PDE object
% =========================================================================== %
classdef PDEObject < FEPack.pdes.LinOperator % FEPack.FEPackObject
  % FEPack.pdes.PDEObject < FEPack.FEPackObject

  properties (SetAccess = protected)

    %> @brief Indicates if normal dierivative or not
    is_normal_derivative = 0;

  end

  methods

    function u = PDEObject   % Initialise PDEObject as 0-order term
      u.alpha = {[1 0 0 0]};
    end

    function v = dual(u)
      v = copy(u);
      v.is_dual = 1;
    end

    function dudn = dn(u)
      dudn = copy(u);
      dudn.is_normal_derivative = 1;
    end

    % For operators
    function op = id(u)   % 0-order term
      op = FEPack.pdes.LinOperator;
      op.is_dual = u.is_dual;
      op.alpha = {[1 0 0 0]};
    end

    function op = dx(u)   % x-derivative
      op = FEPack.pdes.LinOperator;
      op.is_dual = u.is_dual;
      op.alpha = {[0 1 0 0]};
    end

    function op = dy(u)   % y-derivative
      op = FEPack.pdes.LinOperator;
      op.is_dual = u.is_dual;
      op.alpha = {[0 0 1 0]};
    end

    function op = dz(u)   % z-derivative
      op = FEPack.pdes.LinOperator;
      op.is_dual = u.is_dual;
      op.alpha = {[0 0 0 1]};
    end

    function op = grad(u)   % gradient
      op = FEPack.pdes.LinOperator;
      op.is_dual = u.is_dual;
      op.alpha = {[0 1 0 0; 0 0 1 0; 0 0 0 1]};
    end

    function op = gradDir(u, vec)   % directional gradient
      op = vec * grad(u);
    end

    % For boundary conditions
    function ecs = onDomain(u, domain)
      if isa(domain, 'FEPack.meshes.FEDomain')

        % For essential conditions
        ecs = FEPack.pdes.EssentialConditions;
        IdPoints = domain.IdPoints;
        m = size(IdPoints, 1);

        ecs.C = sparse(1:m, IdPoints, 1, m, domain.mesh.numPoints);

      elseif (ischar(domain) ||...
              strcmpi(domain, 'xmin') || strcmpi(domain, 'xmax') ||...
              strcmpi(domain, 'ymin') || strcmpi(domain, 'ymax') ||...
              strcmpi(domain, 'zmin') || strcmpi(domain, 'zmax'))

        % For symbolic boundary conditions
        ecs = FEPack.pdes.SymBC;

        if (u.is_normal_derivative)
          ecs.traceNeumann.scalar.(domain) = 1.0;
        else
          ecs.traceDirichlet.scalar.(domain) = 1.0;
          ecs.traceDirichlet.operator.(domain) = u.fun;
        end

      else

        error('L''operateur | n''est pas compatible avec le terme de droite choisi.');

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
