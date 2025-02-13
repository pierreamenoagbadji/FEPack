%> @file HoneycombObject.m
%> @brief Contains the +applications.HoneycombObject class.
% =========================================================================== %
%> @brief class for examples of honeycomb potentials
% =========================================================================== %
classdef HoneycombObject < FEPack.FEPackObject
  % FEPack.applications.HoneycombObject < FEPack.FEPackObject

  properties (SetAccess = public)

    % Expression of the honeycomb potential
    fun = [];

    % Properties
    is_even = [];
    is_real = [];

    % Center of lattice
    center = [];

    % Lattice vectors
    vecPer1 = [];
    vecPer2 = [];

    % Dual lattice vectors
    dualVec1 = [];
    dualVec2 = [];

    % High symmetry quasi-momenta
    highSymK = [];
    brillouinVerts = [];

  end

  methods

    function obj = HoneycombObject(is_even, is_real, type, fun)

      % is_even (boolean): is the potential even?
      % is_real (boolean): is the potential real-valued?
      % type (string): options are 'atomic', 'optical', 'trigonometric', and 'custom'
      % fun (function handle): only need if the type is 'custom'; gives a user-specified expression of the potential

      obj.is_even = is_even;
      obj.is_real = is_real;

      % Center of lattice
      obj.center = [0; 0];

      % Lattice vectors
      obj.vecPer1 = [0.5*sqrt(3);  0.5];
      obj.vecPer2 = [0.5*sqrt(3); -0.5];

      % Dual lattice vectors
      obj.dualVec1 = (4*pi/sqrt(3)) * [0.5;  0.5*sqrt(3)];
      obj.dualVec2 = (4*pi/sqrt(3)) * [0.5; -0.5*sqrt(3)];

      % High symmetry quasi-momenta
      obj.highSymK = (obj.dualVec1 - obj.dualVec2) / 3; % High-symmetry quasi-momentum
      obj.brillouinVerts = [ obj.highSymK,...
                           -obj.highSymK + obj.dualVec1,...
                            obj.highSymK + obj.dualVec2,...
                           -obj.highSymK,...
                            obj.highSymK - obj.dualVec1,...
                           -obj.highSymK - obj.dualVec2];

      if strcmpi(type, 'atomic')

        % display(obj.myfun(5))
        obj.fun = @(x) obj.atomicPotential(x);

      elseif strcmpi(type, 'optical')

        obj.fun = @(x) obj.opticalPotential(x);

      elseif strcmpi(type, 'trigonometric')

        obj.fun = @(x) obj.trigonometricPolynomial(x);

      elseif strcmpi(type, 'custom')

        if (nargin < 4)
          error(['Option', type, ' was selected: a function handle has to be provided.']);
        end

        obj.fun = fun;

      else

        error(['Option ', type, ' unrecognized; only options are ''atomic'', ''optical'', ''trigonometric'', ''custom''']);

      end

    end

    function val = atomicPotential(obj, x, centers, amps, rads)

      % centers: 2 x Nc vector
      % amps:   Nc x 1  vector
      % rads:   Nc x 1  vector
      if (nargin < 5), rads = 0.4; end
      if (nargin < 4), amps = 10; end
      if (nargin < 3), centers = [0; 0]; end

      vecPer1 = obj.vecPer1(:);
      vecPer2 = obj.vecPer2(:);

      R = [vecPer1, vecPer2];
      T = R \ eye(2);
      normR = norm(R, 'fro');

      Nc = size(centers, 2);
      val = zeros(size(x, 1), 1);

      for idC = 1:Nc

        Xmod = (R * (mod(T * (x(:, 1:2).' - centers(:, idC)) + 0.5, 1) - 0.5)).';
        normXmod = sqrt(Xmod(:, 1).^2 + Xmod(:, 2).^2);

        val = val +...
              amps(idC) * FEPack.tools.cutoff(normXmod, -rads(idC)*normR, rads(idC)*normR);

      end

    end

    function val = opticalPotential(obj, x)
      
      if (obj.is_even)
        Xmod = cos(x(:, 1:2) *  obj.dualVec1) + ...
               cos(x(:, 1:2) *  obj.dualVec2) + ...
               cos(x(:, 1:2) * (obj.dualVec1  + obj.dualVec2));

        val = FEPack.tools.cutoff(Xmod, -1.5, 1.5);
      else
        Xmod = sin(x(:, 1:2) *  obj.dualVec1) + ...
               sin(x(:, 1:2) *  obj.dualVec2) + ...
               sin(x(:, 1:2) * (obj.dualVec1  + obj.dualVec2));

        val = tanh(Xmod);
      end

    end

    function val = trigonometricPolynomial(obj, x)
      
      if (obj.is_even)
        val = cos(x(:, 1:2) *  obj.dualVec1) +...
              cos(x(:, 1:2) *  obj.dualVec2) +...
              cos(x(:, 1:2) * (obj.dualVec1  + obj.dualVec2));
      else
        val = sin(x(:, 1:2) *  obj.dualVec1) +...
              sin(x(:, 1:2) *  obj.dualVec2) +...
              sin(x(:, 1:2) * (obj.dualVec1  + obj.dualVec2));
      end

    end

  end

end
