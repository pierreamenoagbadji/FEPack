%> @file PDEObject.m
%> @brief Contains the pdes.PDEObject class.
% =========================================================================== %
%> @brief class for PDE object
% =========================================================================== %
classdef PDEObject < FEPack.FEPackObject
  % FEPack.pdes.PDEObject < FEPack.FEPackObject

  properties (SetAccess = protected)

  end

  methods (Static)

    % Shape functions
    function phis = shapeFunctions(P, dimension, alpha)
      % SHAPEFUNCTIONS Expression of derivatives of Lagrange P1 shape functions.
      % phis = (P, d, alpha) where P is a point and alpha a (d+1)-vector,
      % returns the linear combination
      %
      %   alpha_1 * phi + \sum_{i = 1}^d alpha_{i+1} * (d_xi phi)
      %
      % where phi is a Lagrange P1 shape function.
      %
      % INPUTS: * P, the points in which the shape functions are evaluated.
      %           P is of size N-by-d
      %         * dimension, the dimension of the points (1, 2, or 3)
      %         * alpha, a (d+1)-vector that contains the coefficients of the
      %           linear combination of the derivatives of the shape functions.
      %
      % OUTPUTS: * phis, (d+1)-by-N matrix containing the shape functions

      d = dimension;

      if (length(alpha) ~= (d+1))
        error('alpha doit avoir le même nombre de colonnes que la dimension du domaine + 1.');
      end

      N = size(P, 1);
      phis = zeros(N, d+1);

      phis(:, 1) = alpha(1) * (1 - P * ones(d, 1)) - ones(N, d) * alpha(2:end)';
      phis(:, 2:end) = alpha(1) * P + kron(ones(N, 1), alpha(2:end));

    end

    % Elementary matrices
    function Aelem = mat_elem(P, domDim, volDim, alpha_u, alpha_v, fun, quadRule)
      % Compute elementary matrix whose components are given by
      %
      %  \int_{T} fun(x) * [alpha_1 * u + \sum_{i = 1}^d alpha_{i+1} * (d_xi u)]
      %                  * [ beta_1 * v + \sum_{i = 1}^d  beta_{i+1} * (d_xi v)]
      %
      % where T is a given element, and where u and v are the Lagrange P1 basis
      % functions associated to the nodes of the element T.
      %
      % INPUTS: * P, 3x(d+1) matrix containing the coordinates of the element's
      %           vertices, where d is the domain dimension.
      %         * domDim, the dimension of the domain.
      %           WARNING: not to be confused with the dimension of the
      %                    volumic geometry. For instance, for a segment,
      %                    domDim = 1.
      %         * volDim, the dimension of the volumic geometry
      %         * alpha_u, a (d+1)-vector that contains the coefficients of the
      %           linear combination of the unknown and its derivatives.
      %         * alpha_v, a (d+1)-vector that contains the coefficients of the
      %           linear combination of the test function and its derivatives.
      %         * fun (function handle) the coefficient in the integral
      %           WARNING: fun must take in argument a nx3 matrix, and return
      %                    a nx1 vector.
      %         * quadRule (QuadratureObject) is optional. For more information,
      %           see +tools/QuadratureObject.m
      %
      % OUTPUTS: Aelem, a (domDim+1)-by-(domDim+1) matrix.

      % Each point must have volDim coordinates at least and 3 coordinates at most
      if ((size(P, 1) < volDim) || (size(P, 1) > 3))
        error(['Les points doivent un nombre de coordonnées égal à 3 au ',...
               'plus, et à la dimension volumique (', int2str(volDim), ') au moins.']);
      end

      % The domain dimension should match the number of points in the element
      if (domDim + 1 ~= size(P, 2))
        error(['La dimension du domaine + 1 (', int2str(domDim + 1), ') et ',...
               'le nombre de points (', int2str(size(P, 2)), ') de l''élément ',...
               'doivent coincider.']);
      end

      % The domain dimension should not exceed the volumic dimension
      if (domDim > volDim)
        error(['La dimension du domaine (', int2str(domDim), ') ne peut ',...
               'pas être supérieure à la dimension volumique (', int2str(volDim), ')']);
      end

      % There should be more than 1 and less or equal than domDim + 1 coefficients
      % for the linear combination
      if (isempty(alpha_u) || (length(alpha_u) > domDim + 1))
        error(['Le nombre de coefficients de la combinaison linéaire de u ',...
              'doit être > 0 et <= à la dimension du domaine + 1.']);
      end
      if (isempty(alpha_v) || (length(alpha_v) > domDim + 1))
        error(['Le nombre de coefficients de la combinaison linéaire de v ',...
              'doit être > 0 et <= à la dimension du domaine + 1.']);
      end

      % Preliminary adjustments and default values
      Lu = length(alpha_u); alpha_u = [alpha_u(:); zeros(domDim+1-Lu, 1)];
      Lv = length(alpha_v); alpha_v = [alpha_v(:); zeros(domDim+1-Lv, 1)];
      dP = size(P); P = [P; zeros(3-dP(1), dP(2))];

      if (nargin < 6)
        fun = @(x) ones(size(x, 1), 1);
      end
      if (nargin < 7)
        quadRule = FEPack.tools.QuadratureObject(domDim);
      end
      Xquad = quadRule.points;
      Wquad = quadRule.weights;

      % Map to reference element
      mapToRel.A = -P(1:volDim, 1) + P(1:volDim, 2:end);
      mapToRel.B =  P(1:volDim, 1);
      switch (domDim)
      case 1
        mapToRel.J = sqrt((P(:, 2) - P(:, 1))' * (P(:, 2) - P(:, 1)));
      case 2
        Tu = cross(P(:, 2) - P(:, 1), P(:, 3) - P(:, 1));
        mapToRel.J = sqrt(Tu' * Tu);
      case 3
        mapToRel.J = abs(det(mapToRel.A));
      end

      % Weighted function
      weightedFun = Wquad .* fun((mapToRel.A * Xquad + mapToRel.B).').';

      % Shape functions
      coeffs.u = alpha_u;
      coeffs.v = alpha_v;
      fieldnames = ['u'; 'v'];

      for idI = 1:2
        beta = coeffs.(fieldnames(idI));

        if (domDim == volDim)
          % 0-order term and partial derivatives
          alpha = [beta(1), (mapToRel.A \ beta(2:end)).'];
        else
          % 0-order term and tangential derivatives
          alpha = [beta(1), beta(2:end).' ./ sqrt(diag(mapToRel.A' * mapToRel.A)).' ];
        end

        phis.(fieldnames(idI)) = FEPack.pdes.PDEObject.shapeFunctions(Xquad.', domDim, alpha);
      end

      % Compute the elementary matrix
      Aelem = mapToRel.J * (phis.v' * diag(weightedFun) * phis.u);

    end

    % Matrix assembly
    function Aglob = assembleFEmatrices(domain, Aloc)
      % ASSEMBLEFEMATRICES Finite elements matrix assembly
      % Aglob = ASSEMBLEFEMATRICES(domain, Aloc) where domain is a
      % FEPack.meshes.FEDomain object, and where Aloc is a function handle
      % computes a sparse matrix that results from the assembly of local
      % matrices computed on the elements of domain, via Aloc.
      %
      % INPUTS: * domain, FEPack.meshes.FEDomain object, the domain for which
      %           the assembly process is performed.
      %         * Aloc a function handle. Given a nx3 matrix corresponding
      %           containg the coordinates of the n points of an element,
      %           Aloc should return the n-by-n associated elementary matrix.
      %
      % OUTPUTS: * Aglob, a N-by-N matrix, where N is the number of DOFs.

      N = domain.mesh.numPoints;   % Number of degrees of freedom
      domDim = domain.dimension;
      dd = (domDim + 1) * (domDim + 1);

      II = zeros(domain.numElts*dd, 1);
      JJ = zeros(domain.numElts*dd, 1);
      VV = zeros(domain.numElts*dd, 1);

      index_II = kron(ones(domDim + 1, 1), (1:domDim + 1).');
      index_JJ = kron((1:domDim + 1).', ones(domDim + 1, 1));

      for ielts = 1:domain.numElts

        % Nodes composing to the element
        P = domain.mesh.points(domain.elements(ielts, :), :).';

        % Elementary matrix associated to the element
        Aelem = Aloc(P);

        % Save the elementary matrix
        index = (dd*(ielts-1)+1):(dd*(ielts-1)+dd);
        II(index) = domain.elements(ielts, index_II);
        JJ(index) = domain.elements(ielts, index_JJ);
        VV(index) = Aelem(:);

      end

      Aglob = sparse(II, JJ, VV, N, N);

    end

    % Global matrices
    function Aglob = global_matrix(varargin)

      % function Aglob = GLOBAL_MATRIX(domain, alpha_u, alpha_v, fun, quadRule)
      % computes the FE matrix whose components are given by
      %
      %  \int_{domain} fun(x) * [alpha_1 * u + \sum_{i = 1}^d alpha_{i+1} * (d_xi u)]
      %                       * [ beta_1 * v + \sum_{i = 1}^d  beta_{i+1} * (d_xi v)]
      %
      % where domain is a given domain, and where u and v are the Lagrange P1
      % basis functions.
      %
      % INPUTS: * domain, FEPack.meshes.FEDomain object, the domain on which
      %           the integrals are evaluated.
      %         * alpha_u, a (d+1)-vector that contains the coefficients of the
      %           linear combination of the unknown and its derivatives.
      %         * alpha_v, a (d+1)-vector that contains the coefficients of the
      %           linear combination of the test function and its derivatives.
      %         * fun (function handle) the coefficient in the integral
      %           WARNING: fun must take in argument a nx3 matrix, and return
      %                    a nx1 vector.
      %         * quadRule (QuadratureObject) is optional. For more information,
      %           see +tools/QuadratureObject.m
      %
      % OUTPUTS: * Aglob, a N-by-N matrix, where N is the number of DOFs.

      dom = varargin{1};

      if (((dom.dimension < dom.mesh.dimension) || (dom.dimension == 0)) && ...
           (~isempty(find(varargin{2}(2:end) == 1, 1)) || ~isempty(find(varargin{3}(2:end) == 1, 1))))
        error(['Calculer l''intégrale d''une dérivée tangentielle sur un ',...
               'domaine surfacique n''est pas autorisé. Pour lever cette ',...
               'erreur, revoir le calcul des matrices élémentaires.']);
      end

      if (dom.dimension ~= 0)

        if (nargin < 4)
          Aloc = @(P) FEPack.pdes.PDEObject.mat_elem(P, dom.dimension, dom.mesh.dimension, varargin{2}, varargin{3});
        elseif (nargin < 5)
          Aloc = @(P) FEPack.pdes.PDEObject.mat_elem(P, dom.dimension, dom.mesh.dimension, varargin{2}, varargin{3}, varargin{4});
        else
          Aloc = @(P) FEPack.pdes.PDEObject.mat_elem(P, dom.dimension, dom.mesh.dimension, varargin{2}, varargin{3}, varargin{4}, varargin{5});
        end

      else

        % Mass matrix on a point
        if (nargin < 4)
          Aloc = @(P) 1.0;
        else
          Aloc = @(P) varargin{4}(P);
        end

      end

      Aglob = FEPack.pdes.PDEObject.assembleFEmatrices(varargin{1}, Aloc);

    end

    % More user-friendly aliases
    function AA = intg_U_V(varargin)
      % function AA = intg_U_V(domain, fun, quadRule)
      AA = FEPack.pdes.PDEObject.global_matrix(varargin{1}, 1, 1, varargin{2:end});
    end

    function AA = intg_U_DxV(varargin)
      % function AA = intg_U_DxV(domain, fun, quadRule)
      AA = FEPack.pdes.PDEObject.global_matrix(varargin{1}, [1 0], [0 1], varargin{2:end});
    end

    function AA = intg_U_DyV(varargin)
      % function AA = intg_U_DyV(domain, fun, quadRule)
      AA = FEPack.pdes.PDEObject.global_matrix(varargin{1}, [1 0 0], [0 0 1], varargin{2:end});
    end

    function AA = intg_U_DzV(varargin)
      % function AA = intg_U_DzV(domain, fun, quadRule)
      AA = FEPack.pdes.PDEObject.global_matrix(varargin{1}, [1 0 0 0], [0 0 0 1], varargin{2:end});
    end

    function AA = intg_DxU_V(varargin)
      % function AA = intg_DxU_V(domain, fun, quadRule)
      AA = FEPack.pdes.PDEObject.global_matrix(varargin{1}, [0 1], [1 0], varargin{2:end});
    end

    function AA = intg_DxU_DxV(varargin)
      % function AA = intg_DxU_DxV(domain, fun, quadRule)
      AA = FEPack.pdes.PDEObject.global_matrix(varargin{1}, [0 1], [0 1], varargin{2:end});
    end

    function AA = intg_DxU_DyV(varargin)
      % function AA = intg_DxU_DyV(domain, fun, quadRule)
      AA = FEPack.pdes.PDEObject.global_matrix(varargin{1}, [0 1 0], [0 0 1], varargin{2:end});
    end

    function AA = intg_DxU_DzV(varargin)
      % function AA = intg_DxU_DzV(domain, fun, quadRule)
      AA = FEPack.pdes.PDEObject.global_matrix(varargin{1}, [0 1 0 0], [0 0 0 1], varargin{2:end});
    end

    function AA = intg_DyU_V(varargin)
      % function AA = intg_DyU_V(domain, fun, quadRule)
      AA = FEPack.pdes.PDEObject.global_matrix(varargin{1}, [0 0 1], [1 0 0], varargin{2:end});
    end

    function AA = intg_DyU_DxV(varargin)
      % function AA = intg_DyU_DxV(domain, fun, quadRule)
      AA = FEPack.pdes.PDEObject.global_matrix(varargin{1}, [0 0 1], [0 1 0], varargin{2:end});
    end

    function AA = intg_DyU_DyV(varargin)
      % function AA = intg_DyU_DyV(domain, fun, quadRule)
      AA = FEPack.pdes.PDEObject.global_matrix(varargin{1}, [0 0 1], [0 0 1], varargin{2:end});
    end

    function AA = intg_DyU_DzV(varargin)
      % function AA = intg_DyU_DzV(domain, fun, quadRule)
      AA = FEPack.pdes.PDEObject.global_matrix(varargin{1}, [0 0 1 0], [0 0 0 1], varargin{2:end});
    end

    function AA = intg_DzU_V(varargin)
      % function AA = intg_DzU_V(domain, fun, quadRule)
      AA = FEPack.pdes.PDEObject.global_matrix(varargin{1}, [0 0 0 1], [1 0 0 0], varargin{2:end});
    end

    function AA = intg_DzU_DxV(varargin)
      % function AA = intg_DzU_DxV(domain, fun, quadRule)
      AA = FEPack.pdes.PDEObject.global_matrix(varargin{1}, [0 0 0 1], [0 1 0 0], varargin{2:end});
    end

    function AA = intg_DzU_DyV(varargin)
      % function AA = intg_DzU_DyV(domain, fun, quadRule)
      AA = FEPack.pdes.PDEObject.global_matrix(varargin{1}, [0 0 0 1], [0 0 1 0], varargin{2:end});
    end

    function AA = intg_DzU_DzV(varargin)
      % function AA = intg_DzU_DzV(domain, fun, quadRule)
      AA = FEPack.pdes.PDEObject.global_matrix(varargin{1}, [0 0 0 1], [0 0 0 1], varargin{2:end});
    end

    function AA = intg_gradU_gradV(varargin)
      % function AA = intg_gradU_gradV(domain, matfun, quadRule)
      % WARNING: Here, matfun is a function handle that must take in argument
      %          a Nx3 vector, and return a (N*d)xd block matrix that looks
      %          for instance like
      %                      (A    C)
      %                      (B    D)
      %          where A is a Nx1 vector which is the (1, 1) component of matfun
      %          applied to the Nx3 input. d is the dimension (which in the
      %          model above equals 2). An example of such function would be
      %
      %                 matfun = @(P) [cos(P(:, 1)), sin(P(:, 2));...
      %                                exp(P(:, 1)), zeros(size(P, 1), 1)]
      dom = varargin{1};
      domDim = dom.dimension;
      if (domDim == 0), domDim = 1; end

      if (nargin < 2)
        % Identity matrix by default
        varargin{2} = @(P) kron(eye(domDim), ones(size(P, 1), 1));
      end

      AA = sparse(0);
      entrees = [varargin(1), 0, 0, varargin(2:end)];

      for idI = 1:domDim
        for idJ = 1:domDim
          % Extract the desired component of the matrix function
          EI = @(P) sparse(1:size(P,1), ((idI-1)*size(P,1)+1):(idI*size(P,1)), 1, size(P,1), domDim*size(P,1));
          EJ = @(P) sparse(idJ, 1, 1, domDim, 1);

          % Modify the inputs
          entrees{2} = zeros(1, domDim+1); entrees{2}(idI+1) = 1;
          entrees{3} = zeros(1, domDim+1); entrees{3}(idJ+1) = 1;
          entrees{4} = @(P) EI(P) * varargin{2}(P) * EJ(P);

          % Update the matrix
          AA = AA + FEPack.pdes.PDEObject.global_matrix(entrees{:});
        end
      end
    end

    function AA = intg_TU_V(domain, Tmat, spectralB, representation)
      % function AA = INTG_TU_V(domain, Tmat, space, representation)
      % Computes the FE matrix associated to the integral
      %
      %       \intg_domain (Tu) * v
      %
      % where T is an operator applied to u. T is represented by a matrix by
      % means of a spectral basis (phi_k)_k. The components of the matrix are
      % the Tmat_{k,l} which are either given by
      %
      %       (1) weak evaluation: Tmat_{k,l} = <Tphi_l, phi_k>; or
      %
      %       (2) projection and representation:
      %             Proj(Tphi_l) = Tmat_{1,l} phi_1 + .... + Tmat_{N, l} phi_N
      %          where Proj(Tphi_l) is the projection of phi_l in the (phi_k)_k.
      %
      % INPUTS: * domain, FEPack.meshes.FEDomain object, the domain on which
      %           the integrals are evaluated.
      %         * Tmat, a matrix that represents the operator applied to the
      %           unknown.
      %         * spectralB, SpectralBasis object.
      %         * representation, a string between 'weak evaluation' and
      %           'projection', which specifies the definition of T.
      %
      % OUTPUTS: * AA, a N-by-N matrix, where N is the number of DOFs.

      if isa(Tmat, 'function_handle')

        % The operator is a multiplication by a function
        AA = intg_U_V(domain, Tmat);

      elseif (length(Tmat) == 1)

        % Trivial case of multiplication by a scalar
        AA = intg_U_V(domain);

      else

        if strcmpi(representation, 'projection')

          % Deduce the weakly evaluated matrix from the projected one
          Tmat = spectralB.massmatInv * Tmat;

        elseif ~strcmpi(representation, 'weak evaluation')

          % Only 'projection' and 'weak evaluation' are allowed
          error(['La variable evaluation ne peut valoir que ''weak',...
                 'evaluation'' ou ''projection''.']);

        end

        % If not done already, compute the matrices associated to
        % the spectral basis
        if isempty(spectralB.massmat)
          spectralB.computeBasisMatrices(0);
        end

        % Deduce the matrix
        Proj = spectralB.massmatInv * spectralB.projmat.';
        AA = Proj' * Tmat * Proj;

      end

    end

  end

  methods

    function v = dual(u, domain)

      v = u;
      v.is_dual = 1;

    end

    % For essential conditions
    function ecs = onDomain(~, domain)

      ecs = FEPack.pdes.EssentialConditions;
      IdPoints = domain.IdPoints;
      m = size(IdPoints, 1);

      ecs.C = sparse(1:m, IdPoints, 1, m, domain.mesh.numPoints);

    end

    function ecs = or(u, domain)

      ecs = onDomain(u, domain);

    end



  end

end
