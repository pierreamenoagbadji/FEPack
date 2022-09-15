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
    alpha = {[]};

    %> @brief Multiplicative coefficient (for unknown)
    fun = {[]};

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
      opRes.alpha = [opA.alpha; opB.alpha];
      opRes.fun = [opA.fun; opB.fun];
    end

    function res = mtimes(lhs, rhs)

      if isa(lhs, 'function_handle')

        % Product of function and operator
        res = copy(rhs);

        for idT = 1:length(rhs.fun)
          if isempty(rhs.fun{idT})
            res.fun{idT} = @(P) lhs(P);
          else
            res.fun{idT} = @(P) lhs(P) .* rhs.fun{idT}(P);
          end
        end

      elseif isa(lhs, 'double')

        % Product of an operator and a scalar
        res = copy(rhs);

        for idT = 1:length(rhs.alpha)
          if (length(lhs) == 1)
            % If lhs is a scalar
            lhsT = lhs;
          else
            % If lhs is a vector or a matrix
            numL = size(lhs, 2);
            lhsT = [lhs, zeros(size(lhs, 1), size(rhs.alpha{idT}, 1)-numL)];
          end
          res.alpha{idT} = lhsT * rhs.alpha{idT};
        end

      elseif isa(lhs, 'FEPack.pdes.LinOperator')

        % Product of operator on unknown and operator on test function.
        if ( lhs.is_dual)
          error('L''opérateur de gauche doit être appliqué à une inconnue.');
        end
        if (~rhs.is_dual)
          error('L''opérateur de droite doit être appliqué à une fonction test.');
        end

        % Function multiplied to an operator on test function is not taken into account
        if max(~cellfun(@isempty, rhs.fun))
          warning('on');
          warning(['Pour une forme bilinéaire, les fonctions multipliées aux ',...
                  'opérateurs sur fonction test ne sont pas prises en compte.']);
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
