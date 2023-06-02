%> @file SymBoundaryCondition.m
%> @brief Contains the pdes.SymBoundaryCondition class.
% =========================================================================== %
%> @brief class for user-defined symbolic essential conditions
% =========================================================================== %
classdef SymBoundaryCondition < FEPack.FEPackObject
  % FEPack.applications.symBC.SymBoundaryCondition < FEPack.FEPackObject

  properties (SetAccess = public)
    
    %> @brief Domain on which the condition is defined
    domain = {};

    %> @brief Right-hand sides
    rhs = [];

    %> @brief Coefficients of trace on boundaries
    gamma0 = {};

    %> @brief Coefficients of normal trace on boundaries
    gamma1 = {};
    
    %> @brief For Robin boundary conditions of the forms
    %> 1. alpha*dn(u) + f(x)*u = b   (fun = f is a function handle)
    %> 2. alpha*dn(u) +   T(u) = b   (fun = T is a matrix) 
    fun = cell(1, 2);

  end


  methods

    function check_domains_compatibility(symBCa, symBCb)
      if ((symBCa.domain{1} == symBCb.domain{1}) == 0)
        error('Les domaines doivent être les mêmes ou se faire face.');
      end
    end

    function rhs_verification(symBC, mode)

      if ((mode == 0) && (~isempty(symBC.rhs)))
        error('Cette opération n''est valable que sur des conditions sans second membre.');
      end
      
      if ((mode == 1) && isempty(symBC.rhs) && ~isempty(symBC.domain))
        error('Cette opération n''est valable que sur des conditions avec second membre.');
      end

    end

    function symBCres = plus(symBCa, symBCb)
      % Preliminary verifications
      check_domains_compatibility(symBCa, symBCb);
      rhs_verification(symBCa, 0);
      rhs_verification(symBCb, 0);

      symBCres = copy(symBCa);

      symBCres.domain = [symBCres.domain; symBCb.domain];
      symBCres.gamma0 = [symBCres.gamma0; symBCb.gamma0];
      symBCres.gamma1 = [symBCres.gamma1; symBCb.gamma1];
      symBCres.fun = [symBCres.fun; symBCb.fun];
      % symBCres.op = [symBCres.op; symBCb.op];


      % symBCres.gamma0 = [symBCa.gamma0; symBCb.gamma0];
      % symBCres.gamma1 = [symBCa.gamma1; symBCb.gamma1];


      % for idI = 1:2

      %   if isempty(symBCres.fun{1, idI})
      %     % BCa has no function
      %     symBCres.fun{1, idI} = symBCb.fun{1, idI};
      %   elseif ~isempty(symBCb.fun{1, idI})
      %     % BCa and BCb both has functions
      %     symBCres.fun{1, idI} = @(P) symBCres.fun{1, idI} + symBCb.fun{1, idI};
      %   end % remains the case where Bca has a function and BCb has no function: nothing changes.

      %   if isempty(symBCres.op{1, idI})
      %     symBCres.op{1, idI} = symBCb.op{1, idI};
      %   elseif ~isempty(symBCb.op{1, idI})
      %     symBCres.op{1, idI} = @(P) symBCres.op{1, idI} + symBCb.op{1, idI};
      %   end

      % end
    end

    function symBCres = mtimes(T, symBC)
      rhs_verification(symBC, 0);
      symBCres = copy(symBC);

      if ~(isa(T, 'double') && (length(T) == 1))
        error('La multiplication ne peut se faire que par un scalaire');
      end

      for idT = 1:length(symBC.domain)
        symBCres.gamma0{idT} = T * symBCres.gamma0{idT}; 
        symBCres.gamma1{idT} = T * symBCres.gamma1{idT};
      end 
      % symBCres.gamma0{1} = T * symBCres.gamma0{1}; 
      % symBCres.gamma1{1} = T * symBCres.gamma1{1}; 
    end

    function symBCres = uminus(symBC)
      symBCres = (-1) * symBC;
    end

    function symBCres = minus(symBCa, symBCb)
      symBCres = symBCa + (-1) * symBCb;
    end

    function symBCres = assignRHS(symBC, rhs)

      if (isa(rhs, 'double') || isa(rhs, 'function_handle') ||...
          isa(rhs, 'FEPack.spaces.SpectralBasis'))
        
        % Set right-hand side
        symBCres = copy(symBC);
        symBCres.rhs = rhs;

        % % Replace empty elements of symBCres.fun and symBCres.op by 1
        % for idI = 1:2
        %   if isempty(symBCres.fun{1, idI})
        %      symBCres.fun{1, idI} = 1;
        %   end

        %   if isempty(symBCres.op{1, idI})
        %      symBCres.op{1, idI} = 1;
        %   end
        % end
        
      else
      
        error('Seuls les scalaires ou les seconds membres sont acceptés.')
      
      end
    end

    function symBCres = eq(symBC, rhs)
    
      symBCres = assignRHS(symBC, rhs);
    
    end

    % AND is for concatenating symbolic boundary conditions
    function symBCres = and(symBCa, symBCb)
      
      if isa(symBCa, 'FEPack.applications.symBC.SymBoundaryCondition')
        if isempty(symBCa.domain), symBCa = {}; else, symBCa = {symBCa}; end
      end

      if isa(symBCb, 'FEPack.applications.symBC.SymBoundaryCondition')
        if isempty(symBCb.domain), symBCb = {}; else, symBCb = {symBCb}; end
      end

      % Make sure symBCa and symBCb are both assigned 
      % a rhs (if they are not empty)
      for idB = 1:length(symBCa)
        rhs_verification(symBCa{idB}, 1)
      end

      for idB = 1:length(symBCb)
        rhs_verification(symBCb{idB}, 1)
      end
      
      % Concatenate conditions
      symBCres = [symBCa; symBCb];
      
    end

  end

end
