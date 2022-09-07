%> @file LinOperator.m
%> @brief Contains the pdes.LinOperator class.
% =========================================================================== %
%> @brief class for PDE object
% =========================================================================== %
classdef LinOperator < FEPack.FEPackObject
  % FEPack.pdes.LinOperator < FEPack.FEPackObject

  properties (SetAccess = public)

    %> @brief Indicates if the operator is on the unknown or not
    is_dual = 0;

    %> @brief n-by-4 matrix alpha such that
    % op(U)i = alpha(i,1)*u + alpha(i,2)*dx(u) + alpha(i,2)*dy(u) + alpha(i,2)*dz(u)
    alpha = [];

    %> @brief Multiplicative coefficient (for unknown)
    fun = @(x) ones(size(x, 1), 1);

  end

  methods (Static)



  end

  methods

    function opRes = plus(opA, opB)
      if (opA.is_dual ~= opB.is_dual)
        error(['Problème de compatibilité : les opérateurs à ajouter ',...
               'doivent être appliquées à des variables de même type (',...
               'soit toutes primales, soit toutes duales).']);
      end
      opRes = copy(opA);
      opRes.alpha = opA.alpha + opB.alpha;
      opRes.fun = @(x) opA.fun(x) + opB.fun(x);
    end

    function res = mtimes(lhs, rhs)
      if isa(lhs, 'function_handle')

        % Product of function and operator
        res = copy(rhs);
        res.fun = @(P) lhs(P) .* rhs.fun(P);

      elseif isa(lhs, 'double')

        % Product of an operator and a scalar
        res = copy(rhs);
        numL = size(lhs, 2);
        lhs = [lhs, zeros(size(lhs, 1), size(rhs.alpha, 1)-numL)];
        res.alpha = lhs * rhs.alpha;

      elseif isa(lhs, 'FEPack.pdes.LinOperator')

        % Product of operator on unknown and operator on test function.
        if ( lhs.is_dual)
          error('L''opérateur de gauche doit être appliqué à une inconnue.');
        end
        if (~rhs.is_dual)
          error('L''opérateur de droite doit être appliqué à une fonction test.');
        end

        res = FEPack.pdes.Form;
        res.alpha_u = lhs.alpha;
        res.alpha_v = rhs.alpha;
        res.fun = lhs.fun;

      else

        error('Vous tentez de multiplier une instance LinOperator par une classe incompatible.');

      end
    end

    function opRes = uminus(op)
      opRes = (-1) * op;
    end

    function opRes = minus(opA, opB)
      opRes = opA + (-1) * opB;
    end

  end

end
