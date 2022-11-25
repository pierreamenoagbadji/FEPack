%> @file PerHalfGuideBVP.m
%> @brief Contains the +solver.PerHalfGuideBVP class.
% =========================================================================== %
%> @brief class for solver of periodic half-guide problem
% =========================================================================== %
classdef PerHalfGuideBVP < FEPack.solver.BVPObject
  % FEPack.solver.PerHalfGuideBVP < FEPack.solver.BVPObject

  properties (SetAccess = protected)

    %> @brief Solution vector associated to positive side
    vecpos = [];

    %> @brief Solution vector associated to negative side
    vecneg = [];

    %> @brief Solution vector associated to interior domain
    vecint = [];

  end

  methods

  end

end
