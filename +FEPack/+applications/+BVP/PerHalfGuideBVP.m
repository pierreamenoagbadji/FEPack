% PerHalfGuideBVP.m
% Contains the +applications.PerHalfGuideBVP class.
% =========================================================================== %
% class for solver of periodic half-guide problem
% =========================================================================== %
classdef PerHalfGuideBVP < FEPack.applications.BVPObject
  % FEPack.applications.PerHalfGuideBVP < FEPack.applications.BVPObject

  properties (SetAccess = protected)

    % Mesh of the periodc cell
    mesh = [];

    % Volume bilinear form
    volBilinearForm = [];

    % Volume linear form
    volLinearForm = [];

    % Symbolic boundary conditions
    boundaryConditions = [];

    % Solution
    solvec = [];

  end

  methods

    % Initialization of the box
    function solbox = PerHalfGuideBVP(dimension, varargin)
      
      if (nargin > 0)
        solbox.initialize(dimension, varargin{:});
      end

    end

    function initialize(solbox, dimension, varargin)
      % INITIALIZE(solbox, dimension)
      % sets up the parameters for a boundary value problem
      % posed in a cell of given dimension.
      %
      % The dimension can be followed by parameter/value pairs 
      % to specify additional properties. Here are the available
      % parameter/value pairs:
      %
      % (*) 'semiInfiniteDirection' (integer between 1 and dimension): 
      %      the infinite dimension of the half-guide.
      % 
      % =========== %
      % Mesh inputs %
      % =========== %
      % (*) 'structured' (logical variable): tells if the mesh is structured or not.
      %
      % (*) 'numEdgeNodes' (integer): number of nodes per cell edge.
      %
      % (*) 'boundingBox' (dimension-by-2 matrix): 
      %     (1D): [xmin, xmax]; (2D): [xmin, xmax; ymin, ymax]; 
      %     (3D): [xmin, xmax; ymin, ymax; zmin, zmax].
      %
      % ========== %      
      % PDE inputs %
      % ========== %      
      % (*) 'volumeBilinearIntegrand' (function_handle @(u, v) F(u, v)),
      %      eg. F(u, v) = grad(u) * grad(v) + u * v.
      %      
      % (*) 'VolumeLinearForm' (double or function_handle @(v) G(v)),
      %      eg. G(v) = v.
      %     
      % (*) 'BoundaryConditions': of the form
      %       [1D]: @(u, xmax, xmin) ... (eg. @(u, xmax, xmin) ((u|xmin) == 0) & ((dn(u)|xmax) == 1))
      %       [2D]: @(u, xmax, xmin, ymax, ymin) ...
      %       [3D]: @(u, xmax, xmin, ymax, ymin, zmax, zmin) ...
      
      % Preliminary check on the dimension
      validDimension = @(d) isnumeric(d) && isscalar(d) && (d > 0) && (d < 4);
      
      if ~validDimension(dimension)
          error('La dimension doit être un scalaire entier valant 1, 2, ou 3.');
      end

      % Default values
      default_semiInfiniteDirection = 1;
      default_numEdgeNodes = (2^(7-dimension)) * ones(1, dimension);
      default_boundingBox = [zeros(dimension, 1), ones(dimension, 1)];
      default_volBilinearIntg = @(u, v) grad(u)*grad(v) + u*v;
      default_volLinearIntg = @(v) v;


      % Validation
      valid_semiInfiniteDirection = @(d) isnumeric(d) && isscalar(d) && (d > 0) && (d < dimension);
      valid_numEdgeNodes = @(N) isnumeric(N) && isvector(N) && (length(N) == dimension);
      valid_boundingBox = @(BB) isnumeric(BB) && ismatrix(BB) && (size(BB, 1) == dimension) && (size(BB, 2) == 2);
      valid_volBilinearIntg = @(blf) isa(blf, 'function_handle');
      valid_volLinearIntg = @(lf) isa(lf, 'function_handle') || (isnumeric(lf) && isscalar(lf));
      

      % Input
      params = inputParser;
      params.addParameter('semiInfiniteDirection', default_semiInfiniteDirection, valid_semiInfiniteDirection);
      %
      params.addParameter('structured', false, @(x) islogical(x));
      params.addParameter('numEdgeNodes', default_numEdgeNodes, valid_numEdgeNodes);
      params.addParameter('boundingBox', default_boundingBox, valid_boundingBox);
      params.addParameter('volumeBilinearIntegrand', default_volBilinearIntg, valid_volBilinearIntg);
      params.addParameter('volumeLinearIntegrand', default_volLinearIntg, valid_volLinearIntg);







      defaultNumNodes = (2^(7-dimension)) * ones(1, dimension);
      defaultBB = [zeros(dimension, 1), ones(dimension, 1)];
      % By default the BVP
      % -Delta(u) + u = 1,
      % combined with homogeneous Neumann conditions on the boundary,
      % is solved
      defaultBilinearForm = @(dom, u, v) FEPack.pdes.Form.intg(dom, grad(u)*grad(v) + u*v);
      defaultLinearForm = @(dom, v) FEPack.pdes.Form.intg(dom, v);
      defaultBC = @(u, varargin) FEPack.applications.symBC.symNeumannBC(dimension, u, varargin{:}); 
      % defaultBasis = 'Lagrange';
      % expectedBases = {'Fourier', 'Lagrange'};
      % defaultDimBasis = defaultNumNodes/2;  % Only for Fourier basis

      % Argument validation
      validNumNodes = @(N) isnumeric(N) && isvector(N) && (length(N) == dimension);
      validBB = @(BB) isnumeric(BB) && ismatrix(BB) && (size(BB, 1) == dimension) && (size(BB, 2) == 2);
      validSymBC = @(symBC) isempty(symBC) || isa(symBC, 'FEPack.applications.symBC.SymBoundaryCondition');
      validBilinearForm = @(blf) isa(blf, 'function_handle');
      validLinearForm = @(lf) isa(lf, 'function_handle') || (isa(lf, 'double') && (length(lf) == 1));
      % validBasis = @(sb) any(validatestring(sb, expectedBases));

      % Add inputs
      params = inputParser;
      params.addParameter('structured', false, @(x) islogical(x));
      params.addParameter('numEdgeNodes', defaultNumNodes, validNumNodes);
      params.addParameter('BoundingBox', defaultBB, validBB);
      params.addParameter('VolumeBilinearForm', defaultBilinearForm, validBilinearForm);
      params.addParameter('VolumeLinearForm', defaultLinearForm, validLinearForm);
      params.addParameter('BoundaryConditions', defaultBC, validSymBC);
      % params.addParameter('SpectralBasis', defaultBasis, validBasis);
      % params.addParameter('dimBasis', defaultDimBasis, validNumNodes);
      
      % Parse the inputs
      parse(params, varargin{:});
      entrees = params.Results;
      
      % Construct the mesh
      switch (dimension)
      case 1
        solbox.mesh = FEPack.meshes.MeshSegment('uniform', entrees.BoundingBox(1), entrees.BoundingBox(2), entrees.numEdgeNodes);
      case 2
        solbox.mesh = FEPack.meshes.MeshRectangle(entrees.structured, entrees.BoundingBox(1, :), entrees.BoundingBox(2, :), entrees.numEdgeNodes(1), entrees.numEdgeNodes(2));
      case 3
        solbox.mesh = FEPack.meshes.MeshCuboid(entrees.structured, entrees.BoundingBox(1, :), entrees.BoundingBox(2, :), entrees.BoundingBox(3, :), entrees.numEdgeNodes(1), entrees.numEdgeNodes(2), entrees.numEdgeNodes(3));
      end

      % % Attach a spectral basis to the domains
      % for idI = 1:dimension
      %   if strcmpi(entrees.SpectralBasis, 'Lagrange')
          
      %     FEPack.spaces.PeriodicLagrangeBasis(solbox.mesh.domains{2*idI-1});
      %     FEPack.spaces.PeriodicLagrangeBasis(solbox.mesh.domains{2*idI});

      %   else

      %     FourierIds = diag(entrees.dimBasis);
      %     FEPack.spaces.FourierBasis(solbox.mesh.domains{2*idI-1}, FourierIds(idI, :));
      %     FEPack.spaces.FourierBasis(solbox.mesh.domains{2*idI},   FourierIds(idI, :));

      %   end
      % end 

      % Forms and boundary conditions
      solbox.volBilinearForm = entrees.VolumeBilinearForm;
      solbox.volLinearForm = entrees.VolumeLinearForm;
      solbox.boundaryConditions = entrees.BoundaryConditions;

      solbox.is_initialized = true;

    end

  end

end
