%> @file ReferenceElement.m
%> @brief Contains the spaces.ReferenceElement class.
%> @brief Lagrange elements
% =========================================================================== %
%> @brief Abstract class for reference elements
%> The finite element considered are the Lagrage Pk finite element.
%> Therefore,
%> 1. The 1D reference element is a segment
%> 2. The 2D reference element is a triangle
%> 3. The 3D reference element is a tetrahedron.
%>
% =========================================================================== %
classdef ReferenceElement < FEPack.FEPackObject
  % FEPack.spaces.ReferenceElement < FEPack.FEPackObject

  properties (SetAccess = protected)

    %> @brief Dimension of reference element
    dimension = 0;

    %> @brief FE Order of reference element
    FEorder = 0;

    %> @brief Number of degrees of freedom associated to element
    numDOFloc = 0;

    %> Points associated to element
    points = [];

    %> Shape coefficients

  end

  methods
    
    function refElt = ReferenceElement(dimension, FEorder)
      % refElt = REFERENCEELEMENT(dimension, FEorder)
      % Constructor for ReferenceElement class
      %
      % INPUTS: * dimension, the dimension of the reference element (between 0 and 3)
      %         * FEorder, Finite element order
      %
      % OUTPUTS: * refElt, ReferenceElement object
      %
      % 
      if (dimension < 0 || dimension > 3)
        error('La dimension doit être comprise entre 0 et 3.');
      end

      refElt.dimension = dimension;
      refElt.FEorder = FEorder;
      
      switch (dimension)
      case 0

        refElt.numDOFloc = 1;
        refElt.points = 1.0;
      
      case 1
      
        refElt.numDOFloc = FEorder + 1;
        refElt.points = [0.0, 1.0, (1:FEorder-1)/FEorder];
      
      case 2
      
        refElt.numDOFloc = (FEorder + 1) * (FEorder + 2) / 2;
        locDOFIds = FEPack.spaces.ReferenceElement.TriangleLocalNumbering(FEorder);
        refElt.points = (locDOFIds - 1.0) / FEorder;

      case 3

        refElt.numDOFloc = (FEorder + 1) * (FEorder + 2) * (FEorder + 3) / 6;
        locDOFIds = FEPack.spaces.ReferenceElement.TetrahedronLocalNumbering(FEorder);
        refElt.points = (locDOFIds - 1.0) / FEorder;

      end

    end

    % function phis = shapeFunctions()
  end

  methods (Static)

    function nums = TriangleLocalNumbering(FEorder)
      % nums = TRIANGLELOCALNUMBERING(FEorder)
      % Provides a numbering for the reference triangle
      % associated with finite elements of order FEorder.
      %
      % To each node composing a triangle element, 
      % one can associate a pair of indexes (i, j) in a
      % manner that corresponds to its placement:
      %
      %    (1,2)       (1,3)               (1,4)
      %     | \         |    \               |   \
      % (1,1)--(2,1)  (1,2) (2,2)          (1,3)  (3,2)
      %   *k = 1*       |         \          |          \
      %               (1,1)--(2,1)--(3,1)  (1,2)  (2,2)  (3,2)
      %                     *k = 2*          |                 \
      %                                    (1,1)--(2,1)--(3,1)--(4,1)
      %                                              *k = 3*
      %
      % TRIANGLELOCALNUMBERING maps this geometric ordering
      % to the following standard numbering:
      %
      %    3      3          3
      %   | \     | \        |  \
      %   1--2    6  5       6   8
      % *k = 1*   |   \      |     \
      %           1--4--2    9  10  5
      %           *k = 2*    |       \
      %                      1--4--7--2
      %                       *k = 3*
      %
      % INPUTS: * FEorder, the order of the finite element triangle.
      %
      % OUTPUTS: * nums, a numDOFloc-by-2 matrix such that
      %             nums(p) = (i(p), j(p)), that is, the
      %             geometrical position of the p-th DOF.

      nums = zeros(0, 2);
      p = FEorder;
      idPmin = 1;

      while (p > 0)
        idPmax = idPmin + p;

        % Triangle vertices that belong to two distinct edges
        nums = [nums; [idPmin, idPmin]]; %#ok
        nums = [nums; [idPmax, idPmin]]; %#ok
        nums = [nums; [idPmin, idPmax]]; %#ok

        % Interior edge nodes (they belong to one edge only)
        edgeId = (idPmin+1:idPmax-1);
        edgeRv = edgeId(end:-1:1);
        onesId = idPmin * ones(1, length(edgeId));

        matI = [edgeId; edgeRv; onesId];
        matJ = [onesId; edgeId; edgeRv];

        nums = [nums; [matI(:), matJ(:)]]; %#ok
        
        % Update
        p = p - 3;
        idPmin = idPmin + 1;
      end

      if (p == 0)
        % P0
        nums = [nums; [idPmin, idPmin]];
      end

    end

    function nums = TetrahedronLocalNumbering(FEorder)
      % nums = TETRAHEDRONLOCALNUMBERING(FEorder) allows to switch
      % Provides a numbering for the reference tetrahedron
      % associated with finite elements of order FEorder.
      %
      % Each DOF can be associated with a triple of integers
      % (i, j, k) that corresponds to its location in the element.
      % This functions allows to switch from this triple to
      % a more standard numbering.
      %
      % Sorry, I'm not good for 3D drawings in code comments,
      % but here is a link that shows a typical numbering
      % (very similar to the one adopted here):
      %
      % https://anum-maths.univ-rennes1.fr/melina/danielmartin/melina/www/doc_html/ef/numlocte.html
      %
      % INPUTS: * FEorder, the order of the finite element triangle.
      %
      % OUTPUTS: * nums, a numDOFloc-by-3 matrix such that
      %             nums(p) = (i(p), j(p), k(p)), that is, the
      %             geometrical position of the p-th DOF.

      nums = zeros(0, 3);
      p = FEorder;
      idPmin = 1;

      while (p > 0)
        idPmax = idPmin + p;

        % Tetrahedron vertices: each one belongs to 3 faces (and 3 edges)
        nums = [nums; [idPmin, idPmin, idPmin]]; %#ok
        nums = [nums; [idPmax, idPmin, idPmin]]; %#ok
        nums = [nums; [idPmin, idPmax, idPmin]]; %#ok
        nums = [nums; [idPmin, idPmin, idPmax]]; %#ok

        % Interior edge nodes: each one belongs to 2 faces
        edgeId = (idPmin+1):(idPmax-1);
        edgeRv = edgeId(end:-1:1);
        onesId = idPmin * ones(1, length(edgeId));

        matI = [edgeId; edgeRv; onesId; onesId; edgeRv; onesId];
        matJ = [onesId; edgeId; edgeRv; onesId; onesId; edgeRv];
        matK = [onesId; onesId; onesId; edgeId; edgeId; edgeId];
        
        nums = [nums; [matI(:), matJ(:), matK(:)]]; %#ok
        
        % Interior face nodes: each one belongs to 1 face
        faceId = idPmin + FEPack.spaces.ReferenceElement.TriangleLocalNumbering(p - 3)';
        onesId = idPmin * ones(1, size(faceId, 2));
        obliqueFaceId = (FEorder+1) - (idPmin-1) - (faceId(1, :)-1) - (faceId(2, :)-1);
        
        matI = [faceId(1, :); faceId(1, :); obliqueFaceId; onesId];
        matJ = [faceId(2, :); onesId;       faceId(1, :);  faceId(1, :)];
        matK = [onesId;       faceId(2, :); faceId(2, :);  faceId(2, :)];

        nums = [nums; [matI(:), matJ(:), matK(:)]]; %#ok

        % Update
        p = p - 4;
        idPmin = idPmin + 1;
      end

      if (p == 0)
        % P0
        nums = [nums; [idPmin, idPmin, idPmin]];
      end

    end

  end
end
