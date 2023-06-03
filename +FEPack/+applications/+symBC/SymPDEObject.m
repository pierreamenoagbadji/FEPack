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

    %> @note The attribute fun is a 'function_handle' (for multiplicative
    %> coefficient) or a matrix (representing a linear operator)
  end

  methods

    function u = SymPDEObject   % Initialise SymPDEObject as 0-order term
      u.alpha = 1;
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
      
      if isa(lhs, 'double') && (length(lhs) == 1)

        % Trivial multiplication by a scalar
        ures.alpha = lhs * u.alpha;

      else

        % Multiplication should be by a function or a scalar
        error(['La multiplication ne peut se faire par une instance ', class(lhs),...
               '; seuls les scalaires sont autorisés.']);
      
      end 

    end

    % For boundary conditions
    function symBC = onDomain(u, domain)
      if isa(domain, 'FEPack.meshes.FEDomain')

        % Symbolic boundary conditions
        symBC = FEPack.applications.symBC.SymBoundaryCondition;
        symBC.domain = {domain};
        
        if (u.is_normal_derivative)
          symBC.gamma0 = {0}; % {[0 0]};
          symBC.gamma1 = {u.alpha}; % {u.alpha};
        else
          symBC.gamma0 = {u.alpha}; % {u.alpha};
          symBC.gamma1 = {0}; % {[0 0]};
        end

      else

        error(['L''operateur | n''est pas compatible avec le terme ', ...
               'de droite choisi (qui est de type ', class(domain), ').']);

      end % if
    end

    function symBC = or(u, domain)
      symBC = onDomain(u, domain);
    end

    % % Boundary condition involving operator represented by matrix
    % function symBC = T_U(u, domain, Tmat, spectralB, representation)
    %   % INPUTS: * u, FEPack.applications.symBC.SymPDEObject.
    %   %         * domain, FEPack.meshes.FEDomain object, the domain on which
    %   %           the condition is defined.
    %   %         * Tmat, a matrix that represents the operator applied to the
    %   %           unknown.
    %   %         * spectralB, SpectralBasis object.
    %   %         * representation, a string between 'weak evaluation' and
    %   %           'projection', which specifies the definition of T.
    %   %
    %   % OUTPUTS: * symBC, FEPack.applications.symBC.SymBoundaryCondition.
    %   normal_derivative_not_allowed(u);
    %   symBC = FEPack.applications.symBC.SymBoundaryCondition;
    %   symBC.domain = {domain};
    %   symBC.gamma0 = {u.alpha};
    %   symBC.gamma1 = {0};
      
    %   if (~isempty(u.fun))

    %     error('l''inconnue semble être déjà multipliée par quelque chose');

    %   else
        
    %     % Set the operator
    %     symBC.fun = {{Tmat, spectralB, representation}};

    %   end

    % end

  end % methods

end
