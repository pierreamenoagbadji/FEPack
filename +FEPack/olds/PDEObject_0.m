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
      %           P is of size Nxd
      %         * dimension, the dimension of the points (1, 2, or 3)
      %         * vec a (d+1)-vector that contains the coefficients of the
      %           linear combination of the derivatives of the shape functions.
      %
      % OUTPUTS: * phis, (d+1)xN matrix containing the shape functions

      if isempty(find((0:dimension) == Iord, 1))
        error(['L''indice de la dérivée partielle doit être un entier ',...
               'inférieur à la dimension (ici ', int2str(dimension), ').']);
      end

      N = size(P, 1);
      phis = zeros(dimension + 1, N);


      if (Iord == 0)

        phis(1, :) = 1;
        for idI = 1:dimension
          phis(idI+1, :) =  P(:, idI).';
          phis(1, :) = phis(1, :) - P(:, idI).';
        end

      else

        phis(1, :) = -ones(1, N);
        phis(Iord + 1, :) = ones(1, N);

      end

    end

    % Elementary matrices
    function Aelem = mat_elem(domDimension, P, Iord_u, Iord_v, fun, quadRule)

      % Compute elementary matrix whose components are given by
      %
      %       \int_{T} fun(x) * (d^m_xi) u * conj((d^n_xj) v)
      %
      % INPUTS: * domDimension, the dimension of the domain.
      %           WARNING: not to be confused with the dimension of the
      %                    volumic geometry. For instance, for a segment,
      %                    domDimension = 1.
      %         * P, (d+1)x3 matrix containing the coordinates of the element's
      %           vertices, where d is the domain dimension.
      %         * Iord_u gives the derivation component for the unknown:
      %           - Iord_u = 0: there is no derivation.
      %           - Iord_u > 0: derivative of the unknown with respect to
      %             the Iord_u-th variable.
      %         * Iord_v is similar to Iord_u, but for test function
      %         * fun (function handle) the coefficient in the integral
      %           WARNING: fun must take in argument a nx3 matrix, and return
      %                    a nx1 vector.
      %         * quadRule (QuadratureObject) is optional. For more information,
      %           see +tools/QuadratureObject.m

      % The dimension should match the number of points in the element
      d = domDimension;
      if (d + 1 ~= size(P, 1))
        error(['La dimension + 1 (', int2str(d + 1), ') et le ',...
               'nombre de points (', int2str(size(P, 1)), ') de l''élément ',...
               'doivent coincider.']);
      end

      if (nargin < 5)
        fun = @(x) ones(size(x, 1), 1);
      end

      % Quadrature rule
      if (nargin < 6)
        quadRule = FEPack.tools.QuadratureObject(domDimension);
      end
      Xquad = quadRule.points;
      Wquad = quadRule.weights;

      % Map to reference element
      Q = P.';
      mapRe.A = -Q(:, 1) + Q(:, 2:end);
      mapRe.B =  Q(:, 1);
      switch (d)
      case 1
        mapRe.J = sqrt((Q(:, 2) - Q(:, 1))' * (Q(:, 2) - Q(:, 1)));
      case 2
        Tu = cross(Q(:, 2) - Q(:, 1), Q(:, 3) - Q(:, 1));
        mapRe.J = sqrt(Tu' * Tu);
      case 3
        mapRe.J = abs(det(mapRe.A));
      end

      % Weighted function
      funWgt = Wquad .* fun((mapRe.A * Xquad + mapRe.B).').';

      % Shape functions
      phis_u = FEPack.pdes.PDEObject.shapeFunctions(Xquad.', d, Iord_u);
      phis_v = FEPack.pdes.PDEObject.shapeFunctions(Xquad.', d, Iord_v);

      if (Iord_u == 0)
        phis_u = FEPack.pdes.PDEObject.shapeFunctions(Xquad.', d, Iord_u);
      else
        invJt = (mapRe.A.') \ eye(size(mapRe.A));
        for Ju = 1:dimension
          phis_u = phis_u + invJt(Iord_u, Ju) * FEPack.pdes.PDEObject.shapeFunctions(Xquad.', d, Iord_u);

        phis_u =
        pu = FEPack.pdes.PDEObject.shapeFunctions(Xquad.', d, Iord_u);


      % Compute the elementary matrix
      Aelem = mapRe.J * (phis_v * diag(funWgt) * phis_u');

    end

    % Matrix assembly
    function Aglob = assembleFEmatrices(domain, Aloc)

      N = domain.mesh.numPoints;   % Number of degrees of freedom
      dim = domain.dimension;
      dd = (dim + 1) * (dim + 1);

      II = zeros(domain.numElts*dd, 1);
      JJ = zeros(domain.numElts*dd, 1);
      VV = zeros(domain.numElts*dd, 1);

      index_II = kron(ones(dim + 1, 1), (1:dim + 1).');
      index_JJ = kron((1:dim + 1).', ones(dim + 1, 1));

      for ielts = 1:domain.numElts

        % Nodes composing to the element
        P = domain.mesh.points(domain.elements(ielts, :), :);

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
    function AA = global_matrix(varargin)

      % function AA = GLOBAL_MATRIX(domain, Iord_u, Iord_v, fun, quadRule)
      dom = varargin{1};

      if (dom.dimension ~= 0)

        if ((dom.dimension < dom.mesh.dimension) && (varargin{2} || varargin{3}))
          warning('on');
          warning(['Vous tentez d''évaluer l''intégrale d''une dérivée ',...
                   'sur une portion du bord. Puisque les fonctions de base ne ',...
                   'sont pas C1 sur de tels domaines, attendez-vous à des ',...
                  'résultats faux.']);
        end

        if (nargin < 4)
          Aloc = @(P) FEPack.pdes.PDEObject.mat_elem(dom.dimension, P, varargin{2}, varargin{3});
        elseif (nargin < 5)
          Aloc = @(P) FEPack.pdes.PDEObject.mat_elem(dom.dimension, P, varargin{2}, varargin{3}, varargin{4});
        else
          Aloc = @(P) FEPack.pdes.PDEObject.mat_elem(dom.dimension, P, varargin{2}, varargin{3}, varargin{4}, varargin{5});
        end

      else

        % Mass matrix on a point
        if (varargin{2} || varargin{3})
          error(['Les intégrales de dérivées sur ',...
                 'des domaines de dimension 0 ne sont pas autorisés.']);
        elseif (nargin < 4)
          Aloc = @(P) 1.0;
        else
          Aloc = @(P) varargin{4}(P);
        end

      end

      AA = FEPack.pdes.PDEObject.assembleFEmatrices(varargin{1}, Aloc);

    end

    % More user-friendly aliases
    function AA = intg_U_V(varargin)
      % function AA = intg_U_V(domain, fun, quadRule)
      AA = FEPack.pdes.PDEObject.global_matrix(varargin{1}, 0, 0, varargin{2:end});
    end

    function AA = intg_U_DxV(varargin)
      % function AA = intg_U_DxV(domain, fun, quadRule)
      AA = FEPack.pdes.PDEObject.global_matrix(varargin{1}, 0, 1, varargin{2:end});
    end

    function AA = intg_U_DyV(varargin)
      % function AA = intg_U_DyV(domain, fun, quadRule)
      AA = FEPack.pdes.PDEObject.global_matrix(varargin{1}, 0, 2, varargin{2:end});
    end

    function AA = intg_U_DzV(varargin)
      % function AA = intg_U_DzV(domain, fun, quadRule)
      AA = FEPack.pdes.PDEObject.global_matrix(varargin{1}, 0, 3, varargin{2:end});
    end

    function AA = intg_DxU_V(varargin)
      % function AA = intg_DxU_V(domain, fun, quadRule)
      AA = FEPack.pdes.PDEObject.global_matrix(varargin{1}, 1, 0, varargin{2:end});
    end

    function AA = intg_DxU_DxV(varargin)
      % function AA = intg_DxU_DxV(domain, fun, quadRule)
      AA = FEPack.pdes.PDEObject.global_matrix(varargin{1}, 1, 1, varargin{2:end});
    end

    function AA = intg_DxU_DyV(varargin)
      % function AA = intg_DxU_DyV(domain, fun, quadRule)
      AA = FEPack.pdes.PDEObject.global_matrix(varargin{1}, 1, 2, varargin{2:end});
    end

    function AA = intg_DxU_DzV(varargin)
      % function AA = intg_DxU_DzV(domain, fun, quadRule)
      AA = FEPack.pdes.PDEObject.global_matrix(varargin{1}, 1, 3, varargin{2:end});
    end

    function AA = intg_DyU_V(varargin)
      % function AA = intg_DyU_V(domain, fun, quadRule)
      AA = FEPack.pdes.PDEObject.global_matrix(varargin{1}, 2, 0, varargin{2:end});
    end

    function AA = intg_DyU_DxV(varargin)
      % function AA = intg_DyU_DxV(domain, fun, quadRule)
      AA = FEPack.pdes.PDEObject.global_matrix(varargin{1}, 2, 1, varargin{2:end});
    end

    function AA = intg_DyU_DyV(varargin)
      % function AA = intg_DyU_DyV(domain, fun, quadRule)
      AA = FEPack.pdes.PDEObject.global_matrix(varargin{1}, 2, 2, varargin{2:end});
    end

    function AA = intg_DyU_DzV(varargin)
      % function AA = intg_DyU_DzV(domain, fun, quadRule)
      AA = FEPack.pdes.PDEObject.global_matrix(varargin{1}, 2, 3, varargin{2:end});
    end

    function AA = intg_DzU_V(varargin)
      % function AA = intg_DzU_V(domain, fun, quadRule)
      AA = FEPack.pdes.PDEObject.global_matrix(varargin{1}, 3, 0, varargin{2:end});
    end

    function AA = intg_DzU_DxV(varargin)
      % function AA = intg_DzU_DxV(domain, fun, quadRule)
      AA = FEPack.pdes.PDEObject.global_matrix(varargin{1}, 3, 1, varargin{2:end});
    end

    function AA = intg_DzU_DyV(varargin)
      % function AA = intg_DzU_DyV(domain, fun, quadRule)
      AA = FEPack.pdes.PDEObject.global_matrix(varargin{1}, 3, 2, varargin{2:end});
    end

    function AA = intg_DzU_DzV(varargin)
      % function AA = intg_DzU_DzV(domain, fun, quadRule)
      AA = FEPack.pdes.PDEObject.global_matrix(varargin{1}, 3, 3, varargin{2:end});
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
      d = dom.dimension;
      if (d == 0), d = 1; end

      if (nargin < 2)
        % Identity matrix by default
        varargin{2} = @(P) kron(eye(d), ones(size(P, 1), 1));
      end

      AA = sparse(0);
      entrees = [varargin(1), 0, 0, varargin(2:end)];

      for idI = 1:d
        for idJ = 1:d
          % Extract the desired component of the matrix function
          EI = @(P) sparse(1:size(P,1), ((idI-1)*size(P,1)+1):(idI*size(P,1)), 1, size(P,1), d*size(P,1));
          EJ = @(P) sparse(idJ, 1, 1, d, 1);

          % Modify the inputs
          entrees{2} = idI;
          entrees{3} = idJ;
          entrees{4} = @(P) EI(P) * varargin{2}(P) * EJ(P);

          % Update the matrix
          AA = AA + FEPack.pdes.PDEObject.global_matrix(entrees{:});
        end
      end
    end

  end

  methods

    function ecs = onDomain(~, domain)

      ecs = FEPack.pdes.EssentialConditions;
      points = domain.points;
      m = size(points, 1);

      ecs.C = sparse(1:m, points, 1, m, domain.mesh.numPoints);

    end

    function ecs = or(FEObj, domain)

      ecs = onDomain(FEObj, domain);

    end



  end

end
