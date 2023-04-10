%> @file SymBoundaryCondition.m
%> @brief Contains the pdes.SymBoundaryCondition class.
% =========================================================================== %
%> @brief class for user-defined symbolic essential conditions
% =========================================================================== %
classdef SymBoundaryCondition < FEPack.FEPackObject
  % FEPack.pdes.SymBoundaryCondition < FEPack.FEPackObject

  properties (SetAccess = public)
    
    %> @brief Domains on which the condition is defined
    domains = {};

    %> @brief Right-hand sides
    rhs = {};

    %> Coefficients of trace on boundaries
    gamma0 = {};

    %> Coefficients of normal trace on boundaries
    gamma1 = {};
    
    %> Function f only for Robin boundary conditions of the form
    %> alpha*dn(u) + f(x)*u = b
    fun = cell(1, 2);

    %> Operator T only for Robin boundary conditions of the form
    %> alpha*dn(u) + T(u) = b
    op = cell(1, 2);

  end


  methods

    function check_domains_compatibility(symBCa, symBCb)
      if ((symBCa.domains{1} == symBCb.domains{1}) == 0)
        error('Les domaines doivent être égaux ou se faire face.');
      end
    end

    function unit_conditions_no_rhs(symBC)
      if ((length(symBC.domains) ~= 1) ||...
         ~(isempty(symBC.rhs)) ||...
          (length(symBC.gamma0) ~= 1) ||...
          (length(symBC.gamma1) ~= 1))
        error('Cette opération n''est valable que sur des conditions unitaires.');
      end
    end

    function symBCres = plus(symBCa, symBCb)
      % Preliminary verifications
      check_domains_compatibility(symBCa, symBCb);
      unit_conditions_no_rhs(symBCa);
      unit_conditions_no_rhs(symBCb);

      symBCres = copy(symBCa);

      symBCres.gamma0 = {symBCa.gamma0{1} + symBCb.gamma0{1}};
      symBCres.gamma1 = {symBCa.gamma1{1} + symBCb.gamma1{1}};
    end

    function symBCres = mtimes(T, symBC)
      unit_conditions_no_rhs(symBC);
      symBCres = copy(symBC);

      symBCres.gamma0{1} = T * symBCres.gamma0{1}; 
      symBCres.gamma1{1} = T * symBCres.gamma1{1}; 
    end

    function symBCres = uminus(symBC)
      symBCres = (-1) * symBC;
    end

    function symBCres = minus(symBCa, symBCb)
      symBCres = symBCa + (-1) * symBCb;
    end

    function symBCres = assignEcs(symBC, rhs)
      unit_conditions_no_rhs(symBC);
      symBCres = copy(symBC);
      symBCres.rhs = {rhs};
    end

    function symBCres = eq(symBC, rhs)
      symBCres = assignEcs(symBC, rhs);
    end

    % AND is for concatenating symbolic boundary conditions
    function symBCres = and(symBCa, symBCb)
      symBCres = copy(symBCa);
      symBCres.domains = [symBCres.domains; symBCb.domains];
      symBCres.rhs = [symBCres.rhs; symBCb.rhs];
      symBCres.gamma0 = [symBCres.gamma0; symBCb.gamma0];
      symBCres.gamma1 = [symBCres.gamma1; symBCb.gamma1];
      symBCres.fun = [symBCres.fun; symBCb.fun];
      symBCres.op = [symBCres.op; symBCb.op];
    end

  end

end
