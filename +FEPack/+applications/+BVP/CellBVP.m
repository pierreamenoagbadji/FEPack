%> @file CellBVP.m
%> @brief Contains the +applications.BVP.CellBVP class.
% =========================================================================== %
%> @brief class for solver of boundary value problem in a cell
% =========================================================================== %
classdef CellBVP < FEPack.applications.BVP.BVPObject
  
  properties (SetAccess = protected)

    %> @brief Mesh of the domain
    mesh = [];

    %> @brief Bilinear form
    volBilinearForm = [];

    %> @brief Linear form
    volLinearForm = [];

    %> @brief Boundary conditions
    boundaryConditions = [];

    %> @brief Solution
    solvec = [];

  end

  % properties (SetAccess = private)

  %   %> @brief Has the object been initialized?
  %   is_initialized = false;

  % end

  methods

    % Initialization of the box
    function solbox = CellBVP(dimension, varargin)
      
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
      % (*) 'structured' (logical variable): tells if the mesh is structured or not.
      %
      % (*) 'numEdgeNodes' (integer): number of nodes per cell edge.
      %
      % (*) 'BoundingBox' (dimension-by-2 matrix): 
      %     (1D): [xmin, xmax]; (2D): [xmin, xmax; ymin, ymax]; 
      %     (3D): [xmin, xmax; ymin, ymax; zmin, zmax].
      %      
      % (*) 'VolumeBilinearForm' (function_handle): it is of the form
      %         @(dom, u, v) FEPack.pdes.Form.intg(dom, F(u, v)),
      %    where dom represents the volume cell and 
      %          F(u, v) is a bilinear function of u and v (eg. F(u, v) = u*v).
      %
      % (*) 'VolumeLinearForm' (double or function_handle): either constant or function of the form
      %         @(dom, v) FEPack.pdes.Form.intg(dom, G(v)),
      %    where dom represents the volume cell and 
      %         G(v) is a linear function of v (eg. G(v) = v).
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

    % Resolution
    function solve(solbox, verbose)

      if (nargin < 2)
        verbose = 1;
      end

      % Make sure the object has been initialized
      if ~(solbox.is_initialized)
        error('L''instance CellBVP doit avoir été initialisée via CellBVP.initialize().');
      end
      
      % Volume FE matrices without surface contributions
      N = solbox.mesh.numPoints;
      u = FEPack.pdes.PDEObject;
      v = dual(u);
      dimension = solbox.mesh.dimension;
      doms = solbox.mesh.domains;
      AA = solbox.volBilinearForm(doms{end}, u, v);
      
      if ((isa(solbox.volLinearForm, 'double'))  && (solbox.volLinearForm == 0))
        LL = sparse(N, 1);
      else
        LL = solbox.volLinearForm(doms{end}, v);
      end

      % Boundary conditions
      % solbox.boundaryConditions is a function_handle
      usym = FEPack.applications.symBC.SymPDEObject;
      symBC = solbox.boundaryConditions(usym, doms{1:2*dimension});
      
      % Make sure there are exactly as many conditions 
      % as the number of faces
      if (numel(symBC) ~= 2*dimension)
        
        error(['On regarde un problème aux limites dans un cube de dimension ', int2str(dimension), '; il faut donc ', int2str(2*dimension), ' conditions aux bords et pas ', int2str(numel(symBC)), '.']);
      
      end

      % ========================================== %
      % Find the boundary conditions corresponding %
      % to each pair of facing boundaries          %
      % ========================================== %
      IdBC = zeros(dimension, 2);
      domVar = {'x', 'y', 'z'};

      for idI = 1:dimension
        is_current_domain = false(1, numel(symBC));

        for idB = 1:numel(symBC)
          is_current_domain(idB) = (symBC{idB}.domain{1} == doms{2*idI-1});
        end
        
        numBCdom = length(find(is_current_domain));
        if (numBCdom ~= 2)

          % Make sure each pair of facing boundaries has exactly 2 conditions 
          error([num2str(numBCdom), ' condition(s) a (ont) été définie(s) sur les bords ', domVar{idI}, '=cte, au lieu de 2 conditions.']);

        end

        IdBC(idI, :) = find(is_current_domain);
      end

      % Extract the surface contributions or essential conditions
      % from each set of boundary conditions
      SS = sparse(N, N);
      LS = sparse(N, 1);
      ecs = FEPack.pdes.EssentialConditions;

      for idI = 1:dimension
        
        % % Current facing boundaries
        % Sigma_A = solbox.mesh.domains{2*idI-1};
        % Sigma_B = solbox.mesh.domains{2*idI};

        % For each pair of facing boundaries t = cst \in {ta, tb}, 
        % the condition is expressed in generic form
        %
        %   BCneu*dU + BCdir*U = BCrhs,
        %  
        % with  dU = [du/dn|ta; du/dn|tb] and 
        %        U = [u|ta; u|tb; fa(x)*u|ta, fb(x)*u|tb, Ta(u), Tb(u)],
        %
        % and where BCdir and BCneu are resp. a 2-by-6 and 
        % 2-by-2 matrices and BCrhs a 2-by-1 vector.
        BCdir = zeros(2, 6);
        BCneu = zeros(2, 2);
        BCrhs =  cell(2, 1);
        opDir  = cell(2, 2);
        funDir = cell(2, 2);
        
        for idK = 1:2 % over conditions
          for idL = 1:numel(symBC{IdBC(idI, idK)}.domain) % over terms
            
            G0 = symBC{IdBC(idI, idK)}.gamma0{idL};
            G1 = symBC{IdBC(idI, idK)}.gamma1{idL};
            idA = (doms{2*idI-1} == symBC{IdBC(idI, idK)}.domain{idL});

            % Neumann term
            BCneu(idK, idA) = BCneu(idK, idA) + G1;
            
            if isempty(symBC{IdBC(idI, idK)}.fun{idL})
              
              % Dirichlet term with no function or operator
              BCdir(idK, idA) = BCdir(idK, idA) + G0;
              
            elseif isa(symBC{IdBC(idI, idK)}.fun{idL}, 'function_handle')
              
              % Dirichlet term with a function
              BCdir(idK, idA+2) = 1;
              if isempty(funDir{idK, idA})
                
                funDir{idK, idA} = @(P) G0 * symBC{IdBC(idI, idK)}.fun{idL}(P);
              
              else
              
                funDir{idK, idA} = @(P) G0 * symBC{IdBC(idI, idK)}.fun{idL}(P) + funDir{idK, idA}(P);
              
              end
            
            else 

              % Dirichlet term with an operator
              BCdir(idK, idA+4) = 1;
              if isempty(opDir{idK, idA})
                
                opDir{idK, idA}{1} = G0 * symBC{IdBC(idI, idK)}.fun{idL}{1};
                opDir{idK, idA}{2} =      symBC{IdBC(idI, idK)}.fun{idL}{2};
                opDir{idK, idA}{3} =      symBC{IdBC(idI, idK)}.fun{idL}{3};

              else

                % Make sure the matrix to be added are compatible
                if ~strcmpi(class(opDir{idK, idA}{2}), class(symBC{IdBC(idI, idK)}.fun{idL}{2}))
                  
                  error(['La condition aux limites fait intervenir ',...
                         'deux opérateurs définis sur le même ',...
                         'domaine : ceux-ci doivent être évalués avec le même type de base.']);

                elseif ~strcmpi(opDir{idK, idA}{3}, symBC{IdBC(idI, idK)}.fun{idL}{3})
                  
                  error(['La condition aux limites fait intervenir ',...
                         'deux opérateurs définis sur le même ',...
                         'domaine : ceux-ci doivent être évalués de la même ',...
                         'manière (soit ''projection'' ou ''interpolation''.)']);

                else
                  
                  opDir{idK, idA}{1} = G0 * symBC{IdBC(idI, idK)}.fun{idL}{1} + opDir{idK, idA}{1};
                  opDir{idK, idA}{2} =      symBC{IdBC(idI, idK)}.fun{idL}{2};
                  opDir{idK, idA}{3} =      symBC{IdBC(idI, idK)}.fun{idL}{3};

                end % if
              end % if

            end % if
          
          end % idL

          BCrhs{idK} = symBC{IdBC(idI, idK)}.rhs;
          if isa(BCrhs{idK}, 'function_handle')

            % Compute the rhs if it is a function handle
            BCrhs{idK} = BCrhs{idK}(solbox.mesh.points);

          else
            
            % The rhs is a scalar
            BCrhs{idK} = BCrhs{idK} * ones(N, 1);

          end % if
        end % idK

        if (rank([BCneu, BCdir]) < 2)
          % ==================================== %
          % There are less than 2 valid          %
          % boundary conditions due to redudancy %
          % ==================================== %
          error(['Faces ', domVar{idI}, '=cte : les conditions aux limites ne sont pas admissibles']);

        end
        
        if (rank(BCneu) == 2)

          % ===================================================== %
          % Case where BCneu is invertible                        %
          % u'|ta and u'|tb therefore can be expressed            %
          % as affine combinations of u|ta, u|tb, etc...          %
          % This corresponds to generic Robin boundary conditions %
          % ===================================================== %
          if (verbose)

            fprintf(['Faces ', domVar{idI}, '=cte : Conditions de Robin reconnues\n']);

          end
          
          invBCneu = BCneu \ eye(2);
          modBCdir = invBCneu * BCdir;

          % % Surface contributions for u|ta, u|tb
          % SSaa = FEPack.pdes.Form.intg(Sigma_A, u * v);
          % SSbb = FEPack.pdes.Form.intg(Sigma_B, u * v);
          % AtoB = sparse(Sigma_A.IdPoints, Sigma_B.IdPoints, 1, N, N);
          
          % SS = SS + modBCdir(1, 1) * SSaa + modBCdir(1, 2) * SSaa * AtoB +...
          %           modBCdir(2, 2) * SSbb + modBCdir(2, 1) * SSbb * AtoB';
          % %
          % LS = LS + invBCneu(1, 1) * SSaa * BCrhs{1} + invBCneu(1, 2) * SSaa * AtoB  * BCrhs{2} +...
          %           invBCneu(2, 2) * SSbb * BCrhs{2} + invBCneu(2, 1) * SSbb * AtoB' * BCrhs{1};

          for idK = 1:2
            for idL = 1:2

              % Surface contributions for u|ta, u|tb
              % Skl := intg(SigmaK, (u|SigmaL)*v), where the
              % obvious identification between SigmaK and SigmaL
              % is made.
              Skl = FEPack.pdes.Form.intg(doms{2*idI-2+idK}, u * v) *...
                    sparse(doms{2*idI-2+idK}.IdPoints, doms{2*idI-2+idL}.IdPoints, 1, N, N);

              SS = SS + modBCdir(idK, idL) * Skl;
              LS = LS + invBCneu(idK, idL) * Skl * BCrhs{idL};
              
              % Surface contributions for fa*u|ta, fb*u|tb
              if ~isempty(funDir{idK, idL})

                Skl = FEPack.pdes.Form.intg(doms{2*idI-2+idK}, funDir{idK, idL} * u * v) *...
                      sparse(doms{2*idI-2+idK}.IdPoints, doms{2*idI-2+idL}.IdPoints, 1, N, N);

                SS = SS + modBCdir(idK, idL) * Skl;
                LS = LS + invBCneu(idK, idL) * Skl * BCrhs{idL};
              
              end

              % Surface contributions for Ta(u)|ta, Tb(u)|tb
              if ~isempty(opDir{idK, idL})

                Skl = FEPack.pdes.Form.intg_TU_V(doms{2*idI-2+idK}, opDir{idK, idL}{1:2}) *...
                      sparse(doms{2*idI-2+idK}.IdPoints, doms{2*idI-2+idL}.IdPoints, 1, N, N);

                SS = SS + modBCdir(idK, idL) * Skl;
                LS = LS + invBCneu(idK, idL) * Skl * BCrhs{idL};
              
              end

            end
          end
                
        elseif (rank(BCneu) == 0)
          
          % =================================================== %
          % Case where BCneu = 0                                %
          % The boundary condition involves only u|ta and u|tb. %
          % This corresponds to Dirichlet boundary conditions.  %
          % =================================================== %
          if (verbose)

            fprintf(['Faces ', domVar{idI}, '=cte : Conditions de Dirichlet reconnues\n']);

          end

          % Make sure there is no multiplicative function
          % or operator involved in the Dirichlet condition
          if (max(max(abs(BCdir(:, 3:end)))) ~= 0)
            error('Les conditions essentielles qui font intervenir une fonction de multiplication ou un opérateur ne sont pas acceptées.');
          end

          % Construct the essential boundary conditions
          invBC0 = BCdir(:, 1:2) \ eye(2);
          ecs = ecs & ((u|doms{2*idI-1}) == (invBC0(1, 1) * BCrhs{1} + invBC0(1, 2) * BCrhs{2}))...
                    & ((u|doms{2*idI  }) == (invBC0(2, 1) * BCrhs{1} + invBC0(2, 2) * BCrhs{2}));

        elseif (rank(BCneu) == 1)

          % ==================================================================== %
          % This a case encompasses variations of (quasi-)periodic conditions    %
          % mixed condition, or many exotic (Cauchy-like) boundary conditions    %
          % which I do not know how (or need) to solve.                          %
          %                                                                      %
          % The idea: since BCneu is of rank 1, it can be written as             %
          %             BCneu = Lv * Rv.',                                       %
          % where Lv and Rv are 2-by-1 vectors. Assume that L1 ≠ 0, and set      %
          %           BCdir = [A1; A2],  where A1, A2 are 1-by-6 vectors.        %
          % Then the boundary conditions can be written as                       %
          %                                                                      %
          %    Rv' * dU + (A1/Lv1) * U = phi1/Lv1              (a natural BC)    %
          %    [A2 - A1*(Lv2/Lv1)] * U = phi2 - phi1*(Lv2/Lv1) (an essential BC) % 
          %                                                                      %
          % A compatibility condition needs to be checked if we do not want      %
          % the normal traces of u to figure in the formulation.                 %
          % ==================================================================== %
          if (verbose)

            fprintf(['Faces ', domVar{idI}, '=cte : Conditions (quasi-)periodiques, mixtes, ou exotiques\n']);

          end

          % Write BCneu as a product of 2-by-1 vectors Lv and Rv
          % Find a non null element in BCneu. Here we choose the maximum
          [~, maxId] = max(abs(BCneu(:)));

          % Find the indices of the max (Im, Km) and the other indices (I0, K0)
          Km = 1 + floor((maxId(1)-1)/2);  % Rv(Km) will be non null
          Im = maxId(1) - 2*(Km - 1);      % Lv(Im) will be non null

          I0 = mod(Im, 2) + 1; % If Im = 1 then I0 = 2, and vice versa
          K0 = mod(Km, 2) + 1; % If Km = 1 then K0 = 2, and vice versa

          % Compute Lv and Rv such that BCneu = Lv * Rv.'
          Lv = zeros(2, 1);
          Lv(Im) = 1;
          Lv(I0) = BCneu(I0, Km) / BCneu(Im, Km);

          Rv = zeros(2, 1);
          Rv(Km) = BCneu(Im, Km);
          Rv(K0) = BCneu(Im, K0);

          % One can deduce the essential condition
          % and express it under the form A0 * u = b0
          % where u = [u|ta; u|tb; fa(x)*u|ta, fb(x)*u|tb, Ta(u), Tb(u)],
          % A0 is a 1-by-6 vector and b0 a scalar
          A0 = BCdir(I0, :) - BCdir(Im, :) * Lv(I0);
          b0 = BCrhs{I0}    - BCrhs{Im}    * Lv(I0);

          % Make sure there is no multiplicative function
          % or operator involved in the Dirichlet condition
          if (max(max(abs(A0(3:end)))) ~= 0)
            error('Les conditions essentielles qui font intervenir une fonction de multiplication ou un opérateur ne sont pas acceptées.');
          end

          % If we do not want the normal traces to appear
          % in the variational formulation, then the test
          % function also need to satisfy a given essential
          % condition.
          % Check if the essential conditions satisfied by
          % the solution and the test function are the same
          % or proportional
          if (rank([A0(1), A0(2); -conj(Rv(2)), -conj(Rv(1))]) ~= 1)

            % The case with non-equivalent essential conditions
            % is not handled
            error(sprintf(['Faces ', domVar{idI}, '=cte : ',...
                           'Conditions de type Cauchy reconnues.\n',...
                           'Si on veut faire disparaitre les traces ', ...
                           'normales de la formulation, les fonctions ', ...
                           'test doivent verifier une condition ', ...
                           'essentielle differente de celle de ', ...
                           'la solution.\nDe telles conditions ne sont ', ...
                           'pas prises en charge.']));  %#ok

          else

            % Essential conditions
            ecs = ecs & ((A0(1) * (u|doms{2*idI-1}) +...
                          A0(2) * (u|doms{2*idI  })) == b0);

            % Surface contributions
            for idL = 1:2
              % Surface contributions for u|ta, u|tb
              % Sml := intg(Sigma_{Km}, (u|SigmaL)*v), where the
              % obvious identification between SigmaK and SigmaL
              % is made.
              Smm = FEPack.pdes.Form.intg(doms{2*idI-2+Km}, u * v);
              KmToL  = sparse(doms{2*idI-2+Km}.IdPoints, doms{2*idI-2+idL}.IdPoints, 1, N, N);
              KmToIm = sparse(doms{2*idI-2+Km}.IdPoints, doms{2*idI-2+Im}.IdPoints,  1, N, N);

              SS = SS + (BCdir(Im, idL)/Rv(Km)) * Smm * KmToL;
              LS = LS + Smm * KmToIm * BCrhs{Im};

              % Surface contributions for fa*u|ta, fb*u|tb
              if ~isempty(funDir{Km, idL})

                Smm = FEPack.pdes.Form.intg(doms{2*idI-2+Km}, funDir{Km, idL} * u * v);
                SS = SS + (BCdir(Im, idL)/Rv(Km)) * Smm * KmToL;
                LS = LS + Smm * KmToIm * BCrhs{Im};
              
              end

              % Surface contributions for Ta(u)|ta, Tb(u)|tb
              if ~isempty(opDir{Km, idL})

                Smm = FEPack.pdes.Form.intg_TU_V(doms{2*idI-2+Km}, opDir{Km, idL}{1:2});
                SS = SS + (BCdir(Im, idL)/Rv(Km)) * Smm * KmToL;
                LS = LS + Smm * KmToIm * BCrhs{Im};
              
              end
            end

          end

        else

          error('Comment es-tu arrivé ici ???');

        end

      end % for idI

      % Apply essential conditions
      if isempty(ecs.C)
        % There is no essential condition
        ecs.P = speye(N);
        ecs.b = sparse(N, 1);
      else
        % The projection matrix has not been computed yet
        ecs.applyEcs;
      end

      % Surface contributions and elimination
      AA0 = ecs.P * (AA + SS) * ecs.P';
      LL0 = ecs.P * (LL + LS);

      % Compute the solution
      solbox.solvec = ecs.P' * (AA0 \ LL0) + ecs.b; 
      
      % Plot solution
      trisurf(solbox.mesh.triangles, solbox.mesh.points(:, 1), solbox.mesh.points(:, 2), real(solbox.solvec));
      hold on;
      view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
      set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);

    end
    % Representation of the output

  end

end
