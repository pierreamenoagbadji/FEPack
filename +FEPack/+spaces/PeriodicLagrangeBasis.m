%> @file PeriodicLagrangeBasis.m
%> @brief Contains the pdes.PeriodicLagrangeBasis class.
% =========================================================================== %
%> @brief Abstract class for all domains
%>
% =========================================================================== %
classdef PeriodicLagrangeBasis < FEPack.spaces.EssentialBCLagrangeBasis
  % FEPack.spaces.PeriodicLagrangeBasis < FEPack.spaces.EssentialBCLagrangeBasis

  properties (SetAccess = protected)

  end

  methods

    function sp = PeriodicLagrangeBasis(domain)
      % PERIODICLAGRANGEBASIS

      sp@FEPack.spaces.EssentialBCLagrangeBasis(domain, @(u, dom) FEPack.spaces.PeriodicLagrangeBasis.periodicBC(u, dom));
      
    end

  end
  
  methods (Static)

    function ecs = periodicBC(u, dom)
      ecs = FEPack.pdes.EssentialConditions;

      for idI = 1:dom.mesh.dimension
        % Periodic conditions on each face
        ecs = ecs & assignEcs((u|dom.mesh.domains{2*idI-1}) - (u|dom.mesh.domains{2*idI}), 0.0);
      end
    end
    
  end

end
