% PerHalfGuideBVP.m
% Contains the +applications.PerHalfGuideBVP class.
% =========================================================================== %
% class for solver of periodic half-guide problem
% =========================================================================== %
classdef PerHalfGuideBVP < FEPack.applications.BVP.BVPObject
  % FEPack.applications.PerHalfGuideBVP < FEPack.applications.BVP.BVPObject

  properties (SetAccess = protected)

    % Infinite direction of the guide
    semiInfiniteDirection = [];

    % Orientation of the guide
    orientation = [];

    % Mesh of the periodc cell
    mesh = [];

    % Volume bilinear integrand
    volumeBilinearIntegrand = [];

    % Volume linear integrand
    volumeLinearIntegrand = [];

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
      % (*) 'orientation' (1 or -1): tells if the guide goes to +infinity or -infinity.
      % 
      % =========== %
      % Mesh inputs %
      % =========== %
      % (*) 'mesh' (FEPack.meshes.Mesh object): pregenerated mesh
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
      % (*) 'volumeLinearIntegrand' (double or function_handle @(v) G(v)),
      %      eg. G(v) = v.
      % 
      % (*) 'boundaryConditions': of the form
      %       [1D]: @(u, xmax, xmin) ... (eg. @(u, xmax, xmin) ((u|xmin) == 0) & ((dn(u)|xmax) == 1))
      %       [2D]: @(u, xmax, xmin, ymax, ymin) ...
      %       [3D]: @(u, xmax, xmin, ymax, ymin, zmax, zmin) ...
      %
      % Preliminary check on the dimension
      validDimension = @(d) isnumeric(d) && isscalar(d) && (d > 0) && (d < 4);
      
      if ~validDimension(dimension)
          error('La dimension doit être un scalaire entier valant 1, 2, ou 3.');
      end

      % Default values
      default_semiInfiniteDirection = 1;
      default_orientation = 1;
      default_mesh = [];
      default_numEdgeNodes = (2^(7-dimension)) * ones(1, dimension);
      default_boundingBox = [zeros(dimension, 1), ones(dimension, 1)];
      default_volBilinearIntg = @(u, v) grad(u) * grad(v) + u * v;
      default_volLinearIntg = 0;
      default_boundaryCdts = @(varargin) FEPack.applications.BVP.PerHalfGuideBVP.DirichletPlusPeriodic(semiInfiniteDirection, dimension, varargin); 
      
      % Validation
      valid_semiInfiniteDirection = @(d) isnumeric(d) && isscalar(d) && (d > 0) && (d <= dimension);
      valid_orientation = @(o) max([-1 1] == o);
      valid_mesh = @(msh) (isa(msh, 'FEPack.meshes.Mesh') && msh.dimension == dimension) || isempty(msh);
      valid_numEdgeNodes = @(N) isnumeric(N) && isvector(N) && (length(N) == dimension);
      valid_boundingBox = @(BB) isnumeric(BB) && ismatrix(BB) && (size(BB, 1) == dimension) && (size(BB, 2) == 2);
      valid_volBilinearIntg = @(blf) isa(blf, 'function_handle');
      valid_volLinearIntg = @(lf) isa(lf, 'function_handle') || (isnumeric(lf) && isscalar(lf));
      valid_boundaryCdts = @(bc) isa(bc, 'function_handle');

      % Input
      params = inputParser;
      params.addParameter('semiInfiniteDirection', default_semiInfiniteDirection, valid_semiInfiniteDirection);
      params.addParameter('orientation', default_orientation, valid_orientation);                                              
      params.addParameter('mesh', default_mesh, valid_mesh);
      params.addParameter('structured', false, @(x) islogical(x));
      params.addParameter('numEdgeNodes', default_numEdgeNodes, valid_numEdgeNodes);
      params.addParameter('boundingBox', default_boundingBox, valid_boundingBox);
      params.addParameter('volumeBilinearIntegrand', default_volBilinearIntg, valid_volBilinearIntg);
      params.addParameter('volumeLinearIntegrand', default_volLinearIntg, valid_volLinearIntg);
      params.addParameter('boundaryConditions', default_boundaryCdts, valid_boundaryCdts);
      
      % Parse the inputs
      parse(params, varargin{:});
      entrees = params.Results;
      
      % Geometric arguments
      solbox.semiInfiniteDirection = entrees.semiInfiniteDirection;
      solbox.orientation = entrees.orientation;

      if ~isempty(entrees.mesh)

        % Pregenerated mesh
        solbox.mesh = entrees.mesh;

      else

        % Generate the mesh
        switch (dimension)
        case 1
          solbox.mesh = FEPack.meshes.MeshSegment('uniform', entrees.BoundingBox(1), entrees.BoundingBox(2), entrees.numEdgeNodes);
        case 2
          solbox.mesh = FEPack.meshes.MeshRectangle(entrees.structured, entrees.BoundingBox(1, :), entrees.BoundingBox(2, :), entrees.numEdgeNodes(1), entrees.numEdgeNodes(2));
        case 3
          solbox.mesh = FEPack.meshes.MeshCuboid(entrees.structured, entrees.BoundingBox(1, :), entrees.BoundingBox(2, :), entrees.BoundingBox(3, :), entrees.numEdgeNodes(1), entrees.numEdgeNodes(2), entrees.numEdgeNodes(3));
        end
      
      end
      
      % Forms and boundary conditions
      solbox.volumeBilinearIntegrand = entrees.volumeBilinearIntegrand;
      solbox.volumeLinearIntegrand = entrees.volumeLinearIntegrand;
      solbox.boundaryConditions = entrees.boundaryConditions;

      solbox.is_initialized = true;

    end

    function solve(solbox)
      % SOLVE(solbox)
      % solve the boundary value problem whose parameters
      % have been previously initialized
      
      % Make sure the object has been initialized
      if ~(solbox.is_initialized)
        error('Il faut d''abord utiliser FEPack.applications.BVP.PerHalfGuideBVP.initialize().');
      end

      % ========= %
      % FE matrix %
      % ========= %
      AA = FEPack.pdes.Form.intg(solbox.mesh.domain('volumic'), solbox.volumeBilinearIntegrand);
      
      % =================== %
      % Boundary conditions %
      % =================== %
      usym = FEPack.applications.symBC.SymPDEObject;
      symBC = solbox.boundaryConditions(usym, doms{1:2*dimension});

      % Find the boundary condition corresponding to the transverse interface
      
    end

  end

  methods (Static)

    function symbc = NeumannPlusPeriodic(semiInfiniteDirection, dimension, u, varargin)

      symbc = FEPack.applications.symBC.SymBoundaryCondition;

      % Neumann Conditions
      idNeu = semiInfiniteDirection;
      symbc = symbc & ((dn(u)|varargin{idNeu}) == 1);
      
      % Periodic conditions
      idPer = 1:dimension; idPer(idNeu) = [];
      for idI = 1:dimension-1
        symbc = symbc & (((    u|varargin{2*idPer(idI)-1}) - (    u|varargin{2*idPer(idI)})) == 0) &...
                        (((dn(u)|varargin{2*idPer(idI)-1}) - (dn(u)|varargin{2*idPer(idI)})) == 0);
      end

    end

    function symbc = DirichletPlusPeriodic(semiInfiniteDirection, dimension, u, varargin)

      symbc = FEPack.applications.symBC.SymBoundaryCondition;

      % Dirichlet Conditions
      idNeu = semiInfiniteDirection;
      symbc = symbc & ((u|varargin{idNeu}) == 1);
      
      % Periodic conditions
      idPer = 1:dimension; idPer(idNeu) = [];
      for idI = 1:dimension-1
        symbc = symbc & (((    u|varargin{2*idPer(idI)-1}) - (    u|varargin{2*idPer(idI)})) == 0) &...
                        (((dn(u)|varargin{2*idPer(idI)-1}) - (dn(u)|varargin{2*idPer(idI)})) == 0);
      end

    end

  end
end
