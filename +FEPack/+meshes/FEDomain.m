%> @file FEDomain.m
%> @brief Contains the meshes.FEDomain class.
% =========================================================================== %
%> @brief Abstract class for all domains
%>
% =========================================================================== %
classdef FEDomain < FEPack.FEPackObject
  % FEPack.meshes.FEDomain < FEPack.FEPackObject

  properties (SetAccess = public)
    
    %> @brief (Optional) Spectral basis attached to the domain
    spectralBasis = [];

  end

  properties (SetAccess = protected)

    %> @brief Dimension of the domain
    dimension = 0;

    %> @brief Mesh
    mesh = [];

    %> @brief Reference of the domain
    reference = 0;

    %> @brief Number of elements within the domain
    numElts = 0;

    %> @brief indices of elements within the domain
    idelements = [];

    %> @brief elements
    elements = [];

    %> @brief Number of points within the domain
    numPoints = 0;

    %> @brief Indices of the points composing the domain
    IdPoints = [];

    % %> @brief Indices of the segments composing the domain
    % segments = [];
    %
    % %> @brief Indices of the triangles composing the domain
    % triangles = [];
    
    %> @brief Number of degreees of freedom associated to
    %> one element composing the domain
    numDOFloc = [];

  end

  methods

    function FEdom = FEDomain(mesh, domain_name, dimension, reference, IdPoints)

      FEdom.name = domain_name;
      FEdom.dimension = dimension;
      FEdom.mesh = mesh;
      FEdom.reference = reference;

      switch (dimension)
      case 0
        
        FEdom.idelements = find(mesh.refPoints == reference);
        FEdom.elements = FEdom.idelements;
        FEdom.numDOFloc = 1;
        
      case 1
        
        FEdom.idelements = find(mesh.refSegments == reference);
        FEdom.elements = mesh.segments(FEdom.idelements, :);
        FEdom.numDOFloc = mesh.FEorder + 1;
        
      case 2
        
        FEdom.idelements = find(mesh.refTriangles == reference);
        FEdom.elements = mesh.triangles(FEdom.idelements, :);
        FEdom.numDOFloc = (mesh.FEorder + 1)*(mesh.FEorder + 2)/2;
        
      case 3
        
        FEdom.idelements = find(mesh.refTetrahedra == reference);
        FEdom.elements = mesh.tetrahedra(FEdom.idelements, :);
        FEdom.numDOFloc = (mesh.FEorder + 1)*(mesh.FEorder + 2)*(mesh.FEorder + 3)/6;
        
      end
      FEdom.numElts = size(FEdom.elements, 1);

      if (nargin < 5)
        IdPoints = unique(FEdom.elements);
      end
      FEdom.IdPoints = IdPoints;
      FEdom.numPoints = length(FEdom.IdPoints);
      % FEdom.points = find(mesh.refPoints == reference);
      % FEdom.segments = find(mesh.refSegments == reference);
      % FEdom.triangles = find(mesh.refTriangles == reference);

    end

    % function attachSpectralBasis(FEdom, spectralB)

    %   % Preliminary verifications
    %   if ((spectralB.domain == FEdom) ~= 1)
    %     error('La base n''est pas définie sur le domaine auquel on souhaite l''attacher.');
    %   end

    %   FEdom.spectralBasis = spectralB;

    % end

    function val = eq(domA, domB)
      % Checks if domA and domB are equal or 
      % are facing each other (last one relevant for cuboids only)
      if (domA.mesh ~= domB.mesh)
        val = 0;
      elseif (domA.reference == domB.reference)
        val = 1;
      else
        [Ia, ~] = find(domA.mesh.mapdomains == domA.reference);
        [Ib, ~] = find(domB.mesh.mapdomains == domB.reference);

        if (Ia == Ib)
          val = 2;
        else
          val = 0;
        end
      end
    end

    function structLoc = locateInDomain(domain, P, almostzero)
      % function [val, structLoc] = LOCATEINDOMAIN(domain, P, almostzero)
      % computes the barycentric coordinates of points P in a given domain.
      %
      % INPUTS: * domain, FEPack.meshes.FEDomain object.
      %         * P, a N-by-3 matrix containing the coordinates of the points.
      %         * almostzero (optional) used to check if a point belongs to an
      %           element or not.
      %
      % OUTPUTS: * structLoc, a structure with fields
      %             -> elements, an N-sized vector. elements(i) represents the
      %                index of the element containing P(i, :).
      %                elements(i) = NaN if P(i, :) does not belong to the
      %                domain.
      %             -> barycoos, an N-by-d matrix that contains the
      %                barycentric coordinates of each point.
      %                barycoos(i,j) = NaN if P(i, :) does not belong to the
      %                domain.

      if (nargin < 4)
        almostzero = 1e-12;
      end

      N = size(P, 1);
      d = domain.dimension + 1;

      structLoc.elements = NaN(N, 1);
      structLoc.barycoos = NaN(N, d);

      idElt = 1;
      X = zeros(d, N);
      isLocated = zeros(N, 1);
      mm = -Inf*ones(1, N);
      fprintf('Localisation dans maillage.\n');

      while (true) % Loop through the elements
        % S = [x1, ...., x_d; y1, ...., y_d; z1, ...., z_d]
        S = domain.mesh.points(domain.elements(idElt, :).', :); % d-by-3

        % The barycentric coordinates are solution of the vector equation
        %       alpha_1*OS_1 + ... + alpha_d*OS_d = OP,
        % which can be written as a linear system.
        A = [S.'; ones(1, d)];   % 4-by-d
        B = [P.'; ones(1, N)];   % 4-by-N

        % In some cases (when domain.dimension < domain.mesh.dimension for
        % instance), A might admit some redundant/zero rows. Find them and
        % delete them
        [Q, R, perm] = qr(A, 'vector');   % A(:, perm) = Q * R
        G = Q' * B;   % Q is a unitary matrix, so Q' = inv(Q)

        % The vector perm is chosen so that ABS(DIAG(R)) is decreasing.
        % In other words, we only need the d-first rows of R as the 4-d next
        % ones are necessarily null.
        R = R(1:d, :);  % d-by-d matrix
        G = G(1:d, :);  % d-by-N matrix

        % Solve the linear system and deduce the barycentric coefficients
        X(perm, :) = R \ G;      % d-by-N matrix

        % Find if the barycentric coordinates are all between 0 and 1 for
        % some points
        critval = min(X.*(1 - X), [], 1);
        isInElement = critval > -almostzero;   % 1-by-N
        structLoc.elements(isInElement) = idElt;
        structLoc.barycoos(isInElement, :) = X(:, isInElement).';
        isLocated(isInElement) = 1;
        mm = max(mm, critval); % Useful when some of the points could not be located.

        % Stop the algorithm if all the points have been located
        if (isLocated)
          fprintf('Tous les points ont été trouvés.\n');
          break;
        else
          idElt = idElt + 1;
        end

        % Stop if all the elements have been tested
        if (idElt > domain.numElts)
          warning(['Certains points ne semblent pas appartenir au domaine. ',...
                   'Les valeurs correspondantes sont NaN.']);
          nanIds = find(isnan(structLoc.elements));
          fprintf('Ids\tValeur critère\n');
          fprintf(['%d\t%0.5e <= ', num2str(-almostzero, '%0.5e'), '\n'], [nanIds.'; mm(nanIds)]);
           warning('Il faudrait peut-être réajuster le seuil du critère.');
          break;
        end
      end

    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % REVOIR: J'ai l'impression qu'on ne l'utilise pas pour le
    % moment. En plus on peut simplifier cette fonction grâce à
    % locateInDomain()
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % function [val, structInterp] = interpolate(domain, P, U, tol)

    %   if (nargin < 4)
    %     tol = 1e-15;
    %   end

    %   if (nargin < 3)
    %     U = sparse(domain.mesh.numPoints, 1);
    %   end

    %   Nx = size(P, 1);
    %   Nu = size(U, 2);
    %   d = domain.mesh.dimension;

    %   if (domain.dimension ~= d)
    %     error('Je ne peux traiter que les cas où les dimensions du domaine et du maillage sont égales.');
    %   end

    %   val = zeros(Nx, Nu);
    %   structInterp.coo = zeros(d+1, Nx);
    %   structInterp.idelements = zeros(1, Nx);

    %   for idP = 1:Nx

    %     % Find the element the current point belongs to as well as its
    %     % barycentric coordinates
    %     idE = 1;
    %     % mm = -Inf;

    %     while (true)
    %       % Nodes of the element
    %       S = domain.mesh.points(domain.elements(idE, :).', :).'; % 3-by-(d+1)

    %       % Compute the barycentric coordinates by solving a linear system
    %       A = [S(1:d, :); ones(1, d + 1)]; % (d+1)-by-(d+1)
    %       B = [P(idP, 1:d).'; 1];      % (d+1)-by-1

    %       baryCoo = A \ B;             % (d+1)-by-1
    %       % mm = max(mm, min(baryCoo.*(1 - baryCoo)));

    %       % Stop criterion
    %       if min(baryCoo.*(1 - baryCoo) > -tol)
    %         break;
    %       else
    %         idE = idE + 1;
    %       end

    %       % Error if point not found
    %       if (idE > domain.numElts)
    %         % mm
    %         error('Le point (%d, %d, %d) ne semble pas appartenir au domaine.', P(idP, 1), P(idP, 2), P(idP, 3));
    %       end
    %     end

    %     % Compute the interpolated value
    %     nodes = domain.elements(idE, :);  % 1-by-(d+1)
    %     val(idP, :) = baryCoo.' * U(nodes', :);
    %     structInterp.coo(:, idP) = baryCoo;
    %     structInterp.idelements(idP) = idE;
    %   end

    % end

  end

end
