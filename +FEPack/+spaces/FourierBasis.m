%> @file FourierBasis.m
%> @brief Contains the pdes.FourierBasis class.
% =========================================================================== %
%> @brief Abstract class for all domains
%>
% =========================================================================== %
classdef FourierBasis < FEPack.spaces.SpectralBasis
  % FEPack.spaces.FourierBasis < FEPack.spaces.SpectralBasis

  properties (SetAccess = protected)

    %> @brief Fourier indices
    FourierIds = [];

  end

  methods (Static)

    function val = FourierBasisFun(P, n, FourierIdsX, FourierIdsY, FourierIdsZ, periods)

      dP = size(P, 2);
      P = [P, zeros(size(P, 1), 3-dP)];

      FourierIdsX = FourierIdsX(:).'; dx = length(FourierIdsX);
      FourierIdsY = FourierIdsY(:).'; dy = length(FourierIdsY);
      FourierIdsZ = FourierIdsZ(:).'; dz = length(FourierIdsZ);

      [Ix, Iy, Iz] = ind2sub([dx, dy, dz], n(:)');

      val = exp(2i*pi*(P(:, 1)*FourierIdsX(Ix)/periods(1) +...
                       P(:, 2)*FourierIdsY(Iy)/periods(2) +...
                       P(:, 3)*FourierIdsZ(Iz)/periods(3)));

    end

  end

  methods

    function sp = FourierBasis(domain, maxId, fine_evaluation)
      % FOURIERBASIS generates Fourier basis functions whose Fourier indices
      % along x (resp. y and z) range from -maxId(1) to maxId(1) (resp.
      % -maxId(2) to maxId(2) and -maxId(3) to maxId(3))
      %
      % INPUTS: * maxId, 3-sized vector that contains the extrema of the
      %           Fourier indices along x, y, and z.
      %           In you do not want Fourier functions along z for instance,
      %           simply set maxId(3) to 0.
      %         * fine_evaluation, integer that indicates how the matrices
      %           are to be constructed.
      %           - If fine_evaluation = 0, the matrices are computed by
      %             approximating each function with its interpolation.
      %           - If fine_evaluation = 1 and if phis is a function handle,
      %             the matrices are computed using a quadrature rule.

      if (nargin < 3)
        fine_evaluation = 0;
      end

      % There should be more than 1 and less or equal than domDim coefficients
      % for the maximum of indices
      if (isempty(maxId) || (length(maxId) > domain.mesh.dimension))
        error(['Le nombre de coefficients du shift ',...
              'doit être > 0 et <= à la dimension du maillage.']);
      end

      if (fine_evaluation)
        warning('on');
        warning('Attention : l''évaluation précise est plus coûteuse et souvent non-nécessaire');
      end

      % Fill the other coordinates of the shift with zeros
      LmaxId = length(maxId);
      maxId = [maxId(:); zeros(3 - LmaxId, 1)];

      sp.domain = domain;

      sp.FourierIds.X = (-maxId(1):maxId(1)); dx = length(sp.FourierIds.X);
      sp.FourierIds.Y = (-maxId(2):maxId(2)); dy = length(sp.FourierIds.Y);
      sp.FourierIds.Z = (-maxId(3):maxId(3)); dz = length(sp.FourierIds.Z);

      periods = max(domain.mesh.points(domain.IdPoints, :), [], 1) -...
                min(domain.mesh.points(domain.IdPoints, :), [], 1);
      periods(periods == 0) = 1;
      
      sp.numBasis = dx * dy * dz;
      sp.is_interpolated = 0;
      sp.phis = @(P, n) FEPack.spaces.FourierBasis.FourierBasisFun(P, n, sp.FourierIds.X, sp.FourierIds.Y, sp.FourierIds.Z, periods);

      sp.massmat = speye(sp.numBasis);

      % Compute projection matrix
      if (~fine_evaluation)

        MM = FEPack.pdes.Form.intg_U_V(sp.domain);
        mm = MM(sp.domain.IdPoints, sp.domain.IdPoints);
        phis_mat = sp.phis(sp.domain.mesh.points(sp.domain.IdPoints, :), 1:sp.numBasis);

        sp.projmat = phis_mat' * mm;

      else

        sp.projmat = zeros(sp.numBasis, sp.domain.numPoints);

        for idK = 1:sp.numBasis

          MM = FEPack.pdes.Form.intg_U_V(sp.domain, @(x) conj(sp.phis(x, idK)));
          mm = MM(sp.domain.IdPoints, sp.domain.IdPoints);
          sp.projmat(idK, :) = (mm * ones(sp.domain.numPoints, 1))';

        end

      end

      % Inverse of mass matrix
      sp.invmassmat = speye(sp.numBasis);

      % FE-to-spectral matrix
      sp.FE_to_spectral = sp.projmat;

      % % Attach spectral basis to domain
      % domain.attachSpectralBasis(sp);

    end

    %% AUCUN INTERET POUR FOURIER
    % function AA = intg_shiftU_V(varargin)
    %   % function AA = intg_shiftU_V(sp, tau, fun, quadRule)
    %   %
    %   % INTG_SHIFTU_V computes the matrix associated to the form
    %   %
    %   %     intg_domain f(x) * phi_j(x + tau) * overline{phi_i(x)} dx
    %   %
    %   % INPUTS: * sp, the FourierBasis object,
    %   %         * tau, n-sized vector where n equals the mesh dimension
    %   %         * fun (function handle) the coefficient in the integral
    %   %           WARNING: fun must take in argument a nx3 matrix, and return
    %   %                    a nx1 vector.
    %   %         * quadRule (QuadratureObject) is optional. For more information,
    %   %           see +tools/QuadratureObject.m.
    %   %           WARNING: if provided, quadRule.points are expected to
    %   %                    span [0, 1].

    %   sp  = varargin{1};
    %   dom = sp.domain;
    %   tau = varargin{2};

    %   % The size of the shift should be more than 1 and less or equal than the
    %   % mesh dimension
    %   if (isempty(tau) || (length(tau) > dom.mesh.dimension))
    %     error(['Le nombre de coefficients du shift ',...
    %           'doit être > 0 et <= à la dimension du maillage.']);
    %   end

    %   % Fill the other coordinates of the shift with zeros
    %   Ltau = length(tau);
    %   tau = [tau(:).', zeros(1, 3 - Ltau)];

    %   if (nargin < 3)

    %     % When no symbol is provided, it is chosen equal to 1 by default,
    %     % in which case the expression of the matrix is simplified.
    %     AA = diag(sp.phis(tau, (1:sp.numBasis)));

    %   else

    %     fun = varargin{3};

    %     if (nargin < 4)
    %       quadRule = FEPack.tools.QuadratureObject(1);
    %     end

    %     % The integrals are computed using a composite rule
    %     coo_names = 'XYZ';

    %     for Icoo = 1:3
    %       FourierId = sp.FourierIds.(coo_names(Icoo));
    %       N = 2 * max(FourierId) + 1;

    %       BB = [min(dom.mesh.points(:, Icoo)), max(dom.mesh.points(:, Icoo))];
    %       compRules.(coo_names(Icoo)) = FEPack.tools.QuadratureObject.composite_segment('uniform', BB(1), BB(2), N, quadRule);
    %     end

    %     dX = length(sp.FourierIds.X); NquadX = compRules.X.numQuad;
    %     dY = length(sp.FourierIds.Y); NquadY = compRules.Y.numQuad;
    %     dZ = length(sp.FourierIds.Z); NquadZ = compRules.Z.numQuad;

    %     % Define the points and weights
    %     points = [kron(compRules.X.points(:), kron(ones(NquadY, 1), ones(NquadZ, 1))),...
    %               kron(ones(NquadX, 1), kron(compRules.Y.points(:), ones(NquadZ, 1))),...
    %               kron(ones(NquadX, 1), kron(ones(NquadY, 1), compRules.Z.points(:)))];

    %     weights = kron(compRules.X.weights(:), kron(ones(NquadY, 1), ones(NquadZ, 1))).*...
    %               kron(ones(NquadX, 1), kron(compRules.Y.weights(:), ones(NquadZ, 1))).*...
    %               kron(ones(NquadX, 1), kron(ones(NquadY, 1), compRules.Z.weights(:)));

    %     % Compute the matrix coefficients
    %     FourierIdsX = kron(sp.FourierIds.X, kron(ones(1, dY), ones(1, dZ)));
    %     FourierIdsY = kron(ones(1, dX), kron(sp.FourierIds.Y, ones(1, dZ)));
    %     FourierIdsZ = kron(ones(1, dX), kron(ones(1, dY), sp.FourierIds.Z));

    %     phis = exp(2i * pi * (points(:, 1) * FourierIdsX +...
    %                           points(:, 2) * FourierIdsY +...
    %                           points(:, 3) * FourierIdsZ));

    %     Intg = weights .* fun(points) .* exp(2i * pi * (tau(1)*FourierIdsX +...
    %                                                     tau(2)*FourierIdsY +...
    %                                                     tau(3)*FourierIdsZ));

    %     AA = phis' * (Intg .* phis);

    %   end

    % end

  end

end
