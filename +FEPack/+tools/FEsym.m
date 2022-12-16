%> @file FEsym.m
%> @brief Contains the tools.FEsym class.
% =========================================================================== %
%> @brief class for symbolic and generic constants
% =========================================================================== %
classdef FEsym < FEPack.FEPackObject
  % FEPack.tools.FEsym < FEPack.FEPackObject

  properties (Constant)

    %> @brief Generic hypercube faces
    xmin_ = 'xmin';
    xmax_ = 'xmax';
    ymin_ = 'ymin';
    ymax_ = 'ymax';
    zmin_ = 'zmin';
    zmax_ = 'zmax';

  end

end
