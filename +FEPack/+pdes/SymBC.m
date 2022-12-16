%> @file SymBC.m
%> @brief Contains the pdes.SymBC class.
% =========================================================================== %
%> @brief class for user-defined symbolic boundary conditions
% =========================================================================== %
classdef SymBC < FEPack.FEPackObject
  % FEPack.pdes.SymBC < FEPack.FEPackObject

  properties (SetAccess = public)

    %> @brief Coefficient of normal trace on cell boundaries
    traceNeumann = struct('scalar', struct(...
                                    'xmin', 0, 'ymin', 0, 'zmin', 0, ...
                                    'xmax', 0, 'ymax', 0, 'zmax', 0));

    %> @brief Coefficient of trace on cell boundaries
    traceDirichlet = struct('scalar',   struct(...
                                        'xmin', 0, 'ymin', 0, 'zmin', 0, ...
                                        'xmax', 0, 'ymax', 0, 'zmax', 0),...
                            'operator', struct(...
                                        'xmin', 0, 'ymin', 0, 'zmin', 0, ...
                                        'xmax', 0, 'ymax', 0, 'zmax', 0));

    %> @brief right-hand side
    rhs = 0;

  end


  methods

    function bcsRes = plus(bcsA, bcsB)

    end


  end

end
