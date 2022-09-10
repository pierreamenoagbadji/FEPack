%> @file SpectralBasis.m
%> @brief Contains the pdes.SpectralBasis class.
% =========================================================================== %
%> @brief Abstract class for all domains
%>
% =========================================================================== %
classdef SpectralBasis < FEPack.FEPackObject
  % FEPack.spaces.SpectralBasis < FEPack.FEPackObject

  properties (SetAccess = protected)

    %> @brief Domain on which the basis functions are defined
    domain = [];

    %> @brief Number of basis functions
    numBasis = 0;

    %> @brief Shows if the basis functions are defined by their values at nodes
    %> or by an expression
    is_interpolated = 0;

    %> @brief Expression of basis functions
    phis = [];

    %> @brief L2 projection matrix onto P1 Lagrange FE spaces
    %> projmat_{i, k} = <w_i, phi_k>_L2, where w_j are Lagrange basis functions
    %> and phi_k are the spectral basis functions
    projmat = [];

    %> @brief Mass matrix
    massmat = [];

    %> @brief FE-to-Spectral matrix. If the FE components (phi(x1),...,phi(xN))
    %> of a function are stored in \vec{\phi}, then FEtoSp * \vec{\phi} returns
    %> the components of its projection in the spectral base.
    %> FE_to_spectral = massmat \ projmat.'
    FE_to_spectral = [];

  end

  methods

    function sp = SpectralBasis(domain, phis, numBasis)
      % SPECTRALBASIS
      %
      % INPUTS: * domain, FEPack.meshes.FEDomain object.
      %         * phis, contains the basis functions.
      %           - If phis is a function_handle, it must be of the form
      %                phis = @(x, n) ....
      %             where x is a Nx-by-1 set of points and n, a 1-by-Nn set of
      %             indices of the function. phis(x, n) returns a Nx-by-Nn
      %             matrix, where Nx and Nn are respectively the lengthes of x
      %             and n.
      %           - If phis is a matrix, then it must be of size N-by-m, where
      %             N is the number of DOFs. In that case, the k-th column
      %             contains the values of phi_k at the nodes of the domain.

      if (nargin > 0)

        % Check the type of basis functions
        if (isa(phis, 'function_handle'))
          is_interpolated = 0;
        else
          if (size(phis, 1) ~= domain.numPoints)
            error(['Le nombre de composantes des fonctions de base spectrales (',...
                   int2str(size(phis, 1)), ') doit coincider avec le nombre de ',...
                   'degrés de liberté associés au domaine (', ...
                   int2str(domain.numPoints), ').']);
          end
          is_interpolated = 1;
        end

        % Set domain and basis functions
        sp.domain = domain;
        sp.phis = phis;
        sp.is_interpolated = is_interpolated;
        sp.numBasis = numBasis;

      end

    end

    function computeBasisMatrices(sp, fine_evaluation)
      % COMPUTEBASISMATRICES computes the mass matrix and the L2
      % projection matrix onto P1 Lagrange FE spaces associated to a basis.
      %
      % INPUTS: * sp, FEPack.pdes.SpectralBasis object
      %         * fine_evaluation, integer that indicates how the matrices
      %           are to be constructed.
      %           - If fine_evaluation = 0, the matrices are computed by
      %             approximating each function with its interpolation.
      %           - If fine_evaluation = 1 and if phis is a function handle,
      %             the matrices are computed using a quadrature rule.

      % Mode of evaluation
      if (nargin < 2)
        fine_evaluation = 0;
      end

      if (fine_evaluation)
        warning('on');
        warning('Attention : l''évaluation précise est plus coûteuse et souvent non-nécessaire');
      end

      % Set matrices
      if (~fine_evaluation || sp.is_interpolated)

        MM = FEPack.pdes.Form.intg_U_V(sp.domain);
        mm = MM(sp.domain.IdPoints, sp.domain.IdPoints);

        if (sp.is_interpolated)
          phis_mat = sp.phis;
        else
          phis_mat = sp.phis(sp.domain.mesh.points(sp.domain.IdPoints, :), 1:sp.numBasis);
        end

        sp.projmat = mm * phis_mat;
        sp.massmat = phis_mat' * sp.projmat;

      else

        sp.projmat = zeros(sp.domain.numPoints, sp.numBasis);
        sp.massmat = zeros(sp.numBasis);

        for idK = 1:sp.numBasis

          MM = FEPack.pdes.Form.intg_U_V(sp.domain, @(x) sp.phis(x, idK));
          mm = MM(sp.domain.IdPoints, sp.domain.IdPoints);
          sp.projmat(:, idK) = mm * ones(sp.domain.numPoints, 1);

          for idL = 1:sp.numBasis

            NN = FEPack.pdes.Form.intg_U_V(sp.domain, @(x) sp.phis(x, idL) .* conj(sp.phis(x, idK)));
            nn = NN(sp.domain.IdPoints, sp.domain.IdPoints);
            sp.massmat(idK, idL) = ones(1, sp.domain.numPoints) * nn * ones(sp.domain.numPoints, 1);

          end

        end

      end

      % FE-to-spectral matrix
      sp.FE_to_spectral = sp.massmat \ sp.projmat.';

    end

  end

end
