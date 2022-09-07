%> @file EssentialConditions.m
%> @brief Contains the pdes.EssentialConditions class.
% =========================================================================== %
%> @brief class for PDE object
% =========================================================================== %
classdef EssentialConditions < FEPack.FEPackObject
  % FEPack.pdes.EssentialConditions < FEPack.FEPackObject

  properties (SetAccess = public)

    %> @brief Constraints matrix
    C = [];

    %> @brief Right-hand side
    rhs = [];

    %> @brief Projection matrix
    P = [];

    %> @brief Reduced constraints right-hand side
    b = [];

  end


  methods

    function ecsRes = plus(ecsA, ecsB)
      ecsRes = copy(ecsA);
      ecsRes.C = ecsA.C + ecsB.C;
    end

    function ecsRes = mtimes(T, ecs)
      ecsRes = copy(ecs);
      ecsRes.C = T * ecs.C;
    end

    function ecsRes = uminus(ecs)
      ecsRes = (-1) * ecs;
    end

    function ecsRes = minus(ecsA, ecsB)
      ecsRes = ecsA + (-1) * ecsB;
    end

    function ecsRes = and(ecsA, ecsB)
      ecsRes = copy(ecsA);
      ecsRes.C = [ecsA.C; ecsB.C];
      ecsRes.rhs = [ecsA.rhs; ecsB.rhs];
    end

    function ecsRes = assignEcs(ecs, rhs)
      ecsRes = copy(ecs);

      if (size(rhs, 1) == 1)
        rhs = ones(size(ecs.C, 1), 1) * rhs;
      end
      ecsRes.rhs = rhs;
    end

    % Essential conditions
    function applyEcs(ecs, almostzero)

      if (nargin < 2)
        almostzero = 1e-12;
      end

      if (size(ecs.C, 1) ~= size(ecs.rhs, 1))
        error(['Le nombre de lignes de la matrice des contraintes (',...
               int2str(size(ecs.C, 1)), ') doit être égal à la taille du second membre (',...
               int2str(size(ecs.rhs, 1)), ').']);
      end

      % Reduce the system of constraints to minimal form
      % ////////////////////////////////////////////////
      % QR decomposition with permutation
      [Q, R, permut] = qr(ecs.C, 'vector');
      G = Q' * ecs.rhs;   % Q is a unitary matrix, so Q' = inv(Q)

      % Find the zeros elements in the diagonal of R. The corresponding lines
      % represent redundant constraints. Thus they have to be removed
      if (size(R, 1) == 1)
        redConst = (abs(R(1, 1)) < almostzero);
      else
        redConst = (abs(diag(R)) < almostzero);
      end

      if (max(redConst))
        % There are redundant constraints
        warning('on');
        warning('Les contraintes redondantes ont été supprimées.');

        % Constraints that are not compatible with the right-hand side
        numIncompConst = length(find(max(abs(G(redConst))) > almostzero));
        if (numIncompConst > 0)
          % Corresponds to constraints of type 0 = b
          warning([int2str(numIncompConst), ' contraintes sont incompatibles ',...
                  'avec le second membre. Elles ont été supprimées.']);
        end

        % Remove constraints
        R(redConst, :) = [];
        G(redConst, :) = [];
      end

      % Construct the projection matrices
      % /////////////////////////////////
      p = size(R, 1);
      N = size(R, 2);
      numRHS = size(G, 2);

      % Invert the upper-diagonal invertible part of R
      X = R(:, 1:p) \ [R(:, p + 1:end), G];
      G = X(:, end-numRHS+1:end);
      X = X(:, 1:end-numRHS);

      % The sets of eliminated and reduced indices
      IdsE = permut(1:p);
      IdsR = permut(1+p:end);

      % Projection matrix and right-hand side contribution
      ecs.P = sparse(1:(N - p), IdsR, 1, N - p, N);
      ecs.P(:, IdsE) = - X.';

      ecs.b = sparse(N, numRHS);
      ecs.b(IdsE, :) = G;

    end

  end

end
