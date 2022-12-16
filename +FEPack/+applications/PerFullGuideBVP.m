%> @file PerHalfGuideBVP.m
%> @brief Contains the +applications.PerHalfGuideBVP class.
% =========================================================================== %
%> @brief class for solver of periodic half-guide problem
% =========================================================================== %
classdef PerHalfGuideBVP < FEPack.applications.BVPObject
  % FEPack.applications.PerHalfGuideBVP < FEPack.applications.BVPObject

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
