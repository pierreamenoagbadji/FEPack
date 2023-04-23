%> @file PDEObject.m
%> @brief Contains the pdes.PDEObject class.
% =========================================================================== %
%> @brief class for PDE object
% =========================================================================== %
classdef PDEObject < FEPack.pdes.LinOperator % FEPack.FEPackObject
  % FEPack.pdes.PDEObject < FEPack.FEPackObject

  properties (SetAccess = protected)

  end

  methods

    function u = PDEObject   % Initialise PDEObject as 0-order term
      u.alpha = {[1 0 0 0]};
    end

    function v = dual(u)
      v = copy(u);
      v.is_dual = 1;
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

    function op = grad(u, dim)   % gradient
      if (dim < 1 || dim > 3)
        error('La dimension doit être comprise entre 1 et 3; %d n''est donc pas autorisé.', dim);
      end
      op = FEPack.pdes.LinOperator;
      op.is_dual = u.is_dual;
      op.alpha = {[zeros(dim, 1), eye(dim), zeros(dim, 3-dim)]};
    end

    function op = grad1(u)
      op = grad(u, 1);
    end

    function op = grad2(u)
      op = grad(u, 2);
    end

    function op = grad3(u)
      op = grad(u, 3);
    end

    % For boundary conditions
    function ecs = onDomain(~, domain)
      if isa(domain, 'FEPack.meshes.FEDomain')

        % For essential conditions
        ecs = FEPack.pdes.EssentialConditions;
        IdPoints = domain.IdPoints;
        m = size(IdPoints, 1);

        ecs.C = sparse(1:m, IdPoints, 1, m, domain.mesh.numPoints);

      else

        error(['L''operateur | n''est pas compatible avec le terme ', ...
               'de droite choisi (qui est de type', class(domain), ').']);

      end % if
    end

    function ecs = or(u, domain)
      ecs = onDomain(u, domain);
    end

  end % methods

end
