%> @file SymPDEObject.m
%> @brief Contains the pdes.SymPDEObject class.
% =========================================================================== %
%> @brief class used for symbolic boundary conditions
% =========================================================================== %
classdef SymPDEObject < FEPack.pdes.PDEObject
  % FEPack.pdes.SymPDEObject < FEPack.pdes.PDEObject

  properties (SetAccess = protected)
    
    %> @brief Indicates if normal dierivative or not
    is_normal_derivative = 0;

    %> @brief Matrix-represented operator
    op = [];

  end

  methods

    function u = SymPDEObject   % Initialise SymPDEObject as 0-order term
      u.alpha = [0 0];
      u.fun = [];
    end

    % Normal derivative
    function dudn = dn(u)
      dudn = copy(u);
      dudn.is_normal_derivative = 1;
    end

    function normal_derivative_not_allowed(u)
      if (u.is_normal_derivative)
        error('Cette opération n''est pas autorisée sur la dérivée normale.');
      end
    end

    % Multiplication 
    function ures = mtimes(lhs, u)
      ures = copy(u);
      
      if isa(lhs, 'function_handle')

        % Multiplication by function
        normal_derivative_not_allowed(u);

        if isempty(u.fun)
          ures.fun = @(P) lhs(P);
        else
          ures.fun = @(P) lhs(P) * ures.fun(P);
        end

      elseif (~isa(lhs, 'double'))

        % Multiplication should be by a function or a matrix (incl. scalar)
        error(['La multiplication ne peut se faire par une instance ', class(lhs),...
               '; seules les types ''function_handle'' et ''double'' sont autorisés.']);
      
      elseif (length(lhs) == 1)
        
        ures.alpha = lhs * u.alpha;
         
      else

        % Operator represented by matrix
        normal_derivative_not_allowed(u);
        
        if isempty(u.op)
          ures.op = lhs;
        else
          ures.op = lhs * ures.op;
        end
        
      end

    end

    % For boundary conditions
    function symBC = onDomain(u, domain)
      if isa(domain, 'FEPack.meshes.FEDomain')

        % Symbolic boundary conditions
        symBC = FEPack.bvpsolvers.SymBoundaryCondition;
        symBC.domains = {domain};
        [~, domId] = find(domain.mesh.mapdomains == domain.reference);
        u.alpha = [0 0];
        u.alpha(domId) = 1;
        
        if (u.is_normal_derivative)
          symBC.gamma0 = {[0 0]};
          symBC.gamma1 = {u.alpha};
        else
          symBC.gamma0 = {u.alpha};
          symBC.gamma1 = {[0 0]};
        end

        symBC.fun = {u.fun};
        symBC.op = {u.op};

      else

        error(['L''operateur | n''est pas compatible avec le terme ', ...
               'de droite choisi (qui est de type', class(domain), ').']);

      end % if
    end

    function symBC = or(u, domain)
      symBC = onDomain(u, domain);
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
