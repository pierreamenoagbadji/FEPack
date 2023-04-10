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
    
  end


  methods

    function check_domains_compatibility(symEcs, symEcsB)
      if ((symEcs.domains{1} == symEcsB.domains{1}) == 0)
        error('Les domaines doivent être égaux ou se faire face.');
      end
    end

    function symEcsRes = plus(symEcsA, symEcsB)
      % Preliminary verifications
      check_domains_compatibility(symEcsA, symEcsB);

      symEcsRes = copy(symEcsA);

      symEcsRes.gamma0 = {symEcsA.gamma0{1} + symEcsB.gamma0{1}};
      symEcsRes.gamma1 = {symEcsA.gamma1{1} + symEcsB.gamma1{1}};
    end

    function symEcsRes = mtimes(T, symEcs)
      symEcsRes = copy(symEcs);
      symEcsRes.gamma0 = {T * symEcs.gamma0{1}};
      symEcsRes.gamma1 = {T * symEcs.gamma1{1}};
    end

    function symEcsRes = uminus(symEcs)
      symEcsRes = (-1) * symEcs;
    end

    function symEcsRes = minus(symEcsA, symEcsB)
      symEcsRes = symEcsA + (-1) * symEcsB;
    end

    function symEcsRes = assignEcs(symEcs, rhs)
      symEcsRes = copy(symEcs);
      symEcsRes.rhs = rhs;
    end

    function symEcsRes = eq(symEcs, rhs)
      symEcsRes = assignEcs(symEcs, rhs);
    end

    % AND is for concatenating symbolic boundary conditions
    function symEcsRes = and(symEcsA, symEcsB)
      symEcsRes = copy(symEcsA);
      symEcsRes.domains = [symEcsRes.domains; symEcsB.domains];
      symEcsRes.rhs = [symEcsRes.rhs; symEcsB.rhs];
      symEcsRes.gamma0 = [symEcsRes.gamma0; symEcsB.gamma0];
      symEcsRes.gamma1 = [symEcsRes.gamma1; symEcsB.gamma1];
    end

  end

end
