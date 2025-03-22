function [dispmesh, val] = edge_states_rational(...
  BBkpar, numNodesK,...
  BBlambda, numNodesLambda,...
  op_pos, numNodesXpos, numNodesYpos, bc_basis_functions_pos,...
  op_neg, numNodesXneg, numNodesYneg, bc_basis_functions_neg,...
  op_int, BBint, numNodesXint, numNodesYint,...
  parallel_use, parallel_numcores, plot_coefficients)
  
  %% [mesh, val] = EDGE_STATES_RATIONAL()

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Meshes
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %% Mesh associated to parallel quasi-momentum / spectral parameter
  dispmesh = FEPack.meshes.MeshRectangle(1, BBkpar, BBlambda, numNodesK, numNodesLambda);

  %% ========= (+) half-guide
  % Mesh
  mesh_pos = FEPack.meshes.MeshRectangle(1, [0 1], [0  1], numNodesXpos, numNodesYpos);
  cell_pos = mesh_pos.domain('volumic');

  % Boundary conditions
  if strcmpi(bc_basis_functions_pos, 'Lagrange')
    init_BCstruct_pos.spB0 = FEPack.spaces.PeriodicLagrangeBasis(mesh_pos.domains{4});
    init_BCstruct_pos.spB1 = FEPack.spaces.PeriodicLagrangeBasis(mesh_pos.domains{3});
  else
    FourierIds = [0 0]; FourierIds(1) = numNodesXpos/4;
    init_BCstruct_pos.spB0 = FEPack.spaces.FourierBasis(mesh_pos.domains{4}, FourierIds);
    init_BCstruct_pos.spB1 = FEPack.spaces.FourierBasis(mesh_pos.domains{3}, FourierIds);
  end

  %% ========= (-) half-guide
  % Mesh
  mesh_neg = FEPack.meshes.MeshRectangle(1, [0 1], [0 -1], numNodesXneg, numNodesYneg);
  cell_neg = mesh_neg.domain('volumic');

  % Boundary conditions
  if strcmpi(bc_basis_functions_neg, 'Lagrange')
    init_BCstruct_neg.spB0 = FEPack.spaces.PeriodicLagrangeBasis(mesh_neg.domains{4});
    init_BCstruct_neg.spB1 = FEPack.spaces.PeriodicLagrangeBasis(mesh_neg.domains{3});
  else
    FourierIds = [0 0]; FourierIds(1) = numNodesXneg/4;
    init_BCstruct_neg.spB0 = FEPack.spaces.FourierBasis(mesh_neg.domains{4}, FourierIds);
    init_BCstruct_neg.spB1 = FEPack.spaces.FourierBasis(mesh_neg.domains{3}, FourierIds);
  end

  %% ========= Interior domain
  mesh_int = FEPack.meshes.MeshRectangle(1, [0 1], BBint, numNodesXint, numNodesYint);
  cell_int  = mesh_int.domain('volumic');
  edge0xInt = mesh_int.domain('xmin');
  edge1xInt = mesh_int.domain('xmax');
  edgeYminInt = mesh_int.domain('ymin');
  edgeYmaxInt = mesh_int.domain('ymax');

  %% ========= Plot coefficients
  if (nargin >= 19 && plot_coefficients)
    % Plot coefficient
    numX = 2;
    numY = 2;

    if isa(op_int.funP, 'function_handle')
      funP = op_int.funP;
    else
      funP = @(x) ones(size(x, 1), 1) * op_int.funP;
    end

    if isa(op_int.funQ, 'function_handle')
      funQ = op_int.funQ;
    else
      funQ = @(x) ones(size(x, 1), 1) * op_int.funQ;
    end
    
    figure;
    set(groot,'defaultAxesTickLabelInterpreter','latex');
    set(groot,'defaulttextinterpreter','latex');
    set(groot,'defaultLegendInterpreter','latex');
    for idX = -numX:numX-1
      for idY = -numY:numY-1
        X = mesh_pos.points(:, 1) + idX;
        Y = mesh_pos.points(:, 2) + idY;

        % P
        subplot(1, 2, 1);
        trisurf(mesh_pos.triangles, X, Y, funP([X, Y]));
        hold on;
        shading interp;
        view(2);
        set(gca, 'DataAspectRatio', [1 1 1], 'FontSize', 16);
        colorbar('TickLabelInterpreter', 'latex');
        xlabel('$P$');

        % Q
        subplot(1, 2, 2);
        trisurf(mesh_pos.triangles, X, Y, funQ([X, Y]));
        hold on;
        shading interp;
        view(2);
        set(gca, 'DataAspectRatio', [1 1 1], 'FontSize', 16);
        colorbar('TickLabelInterpreter', 'latex');
        xlabel('$Q$');
      end
    end

    pause;
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Finite Element matrices
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  u = FEPack.pdes.PDEObject; 
  v = dual(u);

  %% (+) half-guide
  op_pos.gradu_gradv = (op_pos.funP * (op_pos.Tmat' * grad2(u))) * (op_pos.Tmat' * grad2(v));
  op_pos.gradu_vectv = (op_pos.funP * (op_pos.Tmat' * grad2(u))) * (op_pos.vect * v);
  op_pos.vectu_gradv = (op_pos.funP * (op_pos.vect * u)) * (op_pos.Tmat' * grad2(v));
  op_pos.vectu_vectv = (op_pos.funP * (op_pos.vect * u)) * (op_pos.vect * v);
  op_pos.funQ_u_v    = (op_pos.funQ * u) * v;
  op_pos.funR_u_v    = (op_pos.funR * u) * v;

  fprintf('(+) half-guide FE matrices computation\n');
  mat_pos_gradu_gradv = FEPack.pdes.Form.intg(cell_pos, op_pos.gradu_gradv);
  mat_pos_gradu_vectv = FEPack.pdes.Form.intg(cell_pos, op_pos.gradu_vectv);
  mat_pos_vectu_gradv = FEPack.pdes.Form.intg(cell_pos, op_pos.vectu_gradv);
  mat_pos_vectu_vectv = FEPack.pdes.Form.intg(cell_pos, op_pos.vectu_vectv);
  mat_pos_funQ_u_v    = FEPack.pdes.Form.intg(cell_pos, op_pos.funQ_u_v);
  mat_pos_funR_u_v    = FEPack.pdes.Form.intg(cell_pos, op_pos.funR_u_v);
  fprintf('(+) half-guide FE matrices computation done.\n');

  %% (-) half-guide
  op_neg.gradu_gradv = (op_neg.funP * (op_neg.Tmat' * grad2(u))) * (op_neg.Tmat' * grad2(v));
  op_neg.gradu_vectv = (op_neg.funP * (op_neg.Tmat' * grad2(u))) * (op_neg.vect * v);
  op_neg.vectu_gradv = (op_neg.funP * (op_neg.vect * u)) * (op_neg.Tmat' * grad2(v));
  op_neg.vectu_vectv = (op_neg.funP * (op_neg.vect * u)) * (op_neg.vect * v);
  op_neg.funQ_u_v    = (op_neg.funQ * u) * v;
  op_neg.funR_u_v    = (op_neg.funR * u) * v;

  fprintf('(-) half-guide FE matrices computation\n');
  mat_neg_gradu_gradv = FEPack.pdes.Form.intg(cell_neg, op_neg.gradu_gradv);
  mat_neg_gradu_vectv = FEPack.pdes.Form.intg(cell_neg, op_neg.gradu_vectv);
  mat_neg_vectu_gradv = FEPack.pdes.Form.intg(cell_neg, op_neg.vectu_gradv);
  mat_neg_vectu_vectv = FEPack.pdes.Form.intg(cell_neg, op_neg.vectu_vectv);
  mat_neg_funQ_u_v    = FEPack.pdes.Form.intg(cell_neg, op_neg.funQ_u_v);
  mat_neg_funR_u_v    = FEPack.pdes.Form.intg(cell_neg, op_neg.funR_u_v);
  fprintf('(-) half-guide FE matrices computation done.\n');

  %% Interior domain
  op_int.gradu_gradv = (op_int.funP * (op_int.Tmat' * grad2(u))) * (op_int.Tmat' * grad2(v));
  op_int.gradu_vectv = (op_int.funP * (op_int.Tmat' * grad2(u))) * (op_int.vect * v);
  op_int.vectu_gradv = (op_int.funP * (op_int.vect * u)) * (op_int.Tmat' * grad2(v));
  op_int.vectu_vectv = (op_int.funP * (op_int.vect * u)) * (op_int.vect * v);
  op_int.funQ_u_v    = (op_int.funQ * u) * v;
  op_int.funR_u_v    = (op_int.funR * u) * v;

  fprintf('Interior domain FE matrices computation\n');
  mat_int_gradu_gradv = FEPack.pdes.Form.intg(cell_int, op_int.gradu_gradv);
  mat_int_gradu_vectv = FEPack.pdes.Form.intg(cell_int, op_int.gradu_vectv);
  mat_int_vectu_gradv = FEPack.pdes.Form.intg(cell_int, op_int.vectu_gradv);
  mat_int_vectu_vectv = FEPack.pdes.Form.intg(cell_int, op_int.vectu_vectv);
  mat_int_funQ_u_v    = FEPack.pdes.Form.intg(cell_int, op_int.funQ_u_v);
  mat_int_funR_u_v    = FEPack.pdes.Form.intg(cell_int, op_int.funR_u_v);
  fprintf('Interior domain FE matrices computation done.\n');

  ecs_int = (((u|edge0xInt) - (u|edge1xInt)) == 0.0);
  ecs_int.applyEcs;
  PPint = ecs_int.P;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Edge states
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Variables
  opts.computeSol = false;
  opts.solBasis = false;
  opts.verbose = 0;
  opts.stop_if_spectral_radius_is_one = true;
  err_prefix = '[propagationOperators.m][spectral_radius_is_one]';
  num_err = length(err_prefix);

  % Output initialization
  numD = dispmesh.numPoints;
  val = zeros(numD, 1);

  % Loop
  if (parallel_use)

    if (nargin < 18)
      error('parallel_use has been set to true; you have to specify the number of cores through parallel_numcores. Use the command feature(''numcores'') to see the number of available cores.');
    end

    % Parallel version
    pool = gcp('nocreate'); % Get current parallel pool
    if isempty(pool)
      parpool('local', parallel_numcores); % Start a parallel pool if none exists
    end
    fprintf('Using parfor (parallel execution)\n');

    kpars   = dispmesh.points(:, 1);
    lambdas = dispmesh.points(:, 2);

    parfor idI = 1:numD
      fprintf('%d sur %d\n', idI, numD);

      kpar   = kpars(idI);
      lambda = lambdas(idI); 

      try

        % ========= (+) half-guide
        % FE matrix
        AApos = mat_pos_gradu_gradv...
              + 1i * kpar     * mat_pos_vectu_gradv...
              - 1i * kpar     * mat_pos_gradu_vectv...
              + (kpar * kpar) * mat_pos_vectu_vectv...
              + mat_pos_funQ_u_v...
              - lambda * mat_pos_funR_u_v;

        % Boundary condition
        BCstruct_pos = init_BCstruct_pos;
        BCstruct_pos.BCu  = 1i*lambda;
        BCstruct_pos.BCdu = -1;
        %
        % BCstruct_pos.BCu  = 1;
        % BCstruct_pos.BCdu = 0;

        % Half-guide problem
        [~, ~, ~, ~, ~, BCstruct_pos, boundary_op_pos] = PeriodicHalfGuideBVP(mesh_pos, +1, 2, AApos, BCstruct_pos, 0, opts);

        % ========= (-) half-guide
        % FE matrix
        AAneg = mat_neg_gradu_gradv...
              + 1i * kpar     * mat_neg_vectu_gradv...
              - 1i * kpar     * mat_neg_gradu_vectv...
              + (kpar * kpar) * mat_neg_vectu_vectv...
              + mat_neg_funQ_u_v...
              - lambda * mat_neg_funR_u_v;

        % Boundary condition
        BCstruct_neg = init_BCstruct_neg;
        BCstruct_neg.BCu  = -1i*lambda;
        BCstruct_neg.BCdu = 1;
        %
        % BCstruct_neg.BCu  = 1;
        % BCstruct_neg.BCdu = 0;

        % Half-guide problem
        [~, ~, ~, ~, ~, BCstruct_neg, boundary_op_neg] = PeriodicHalfGuideBVP(mesh_neg, -1, 2, AAneg, BCstruct_neg, 0, opts);

        % ========= Interior problem
        AAint = mat_int_gradu_gradv +...
              +  1i * kpar    * mat_int_vectu_gradv...
              -  1i * kpar    * mat_int_gradu_vectv...
              + (kpar * kpar) * mat_int_vectu_vectv...
              + mat_int_funQ_u_v...
              - lambda * mat_int_funR_u_v;

        % Coefficients associated to boundary condition
        Nbpos = size(boundary_op_pos, 1);
        Nbneg = size(boundary_op_neg, 1);

        BCu_int_pos = (BCstruct_pos.BCu' * eye(Nbpos) + BCstruct_pos.BCdu * boundary_op_pos)...
                    \ (boundary_op_pos * BCstruct_pos.BCu - BCstruct_pos.BCdu' * eye(Nbpos));
        BCu_int_neg = (BCstruct_neg.BCu' * eye(Nbneg) + BCstruct_neg.BCdu * boundary_op_neg)...
                    \ (boundary_op_neg * BCstruct_neg.BCu - BCstruct_neg.BCdu' * eye(Nbneg));

        SSpos = FEPack.pdes.Form.intg_TU_V(edgeYmaxInt, BCu_int_pos, BCstruct_pos.spB0, 'projection');
        SSneg = FEPack.pdes.Form.intg_TU_V(edgeYminInt, BCu_int_neg, BCstruct_neg.spB0, 'projection');

        AAint0 = PPint * (AAint - SSpos - SSneg) * PPint.';
        
        % Compute infinity norm of resolvent
        val(idI) = max(max(abs(AAint0 \ eye(size(AAint0)))));

      catch err_obj
        
        % val(idI) = 0;
        if (length(err_obj.message) >= num_err) &&...
          (strcmp(err_obj.message(1:num_err), err_prefix))
          % Some eigenvalues of the Riccati operator are on the unit circle
          val(idI) = 0;
        else
          % Unknown error
          rethrow(err_obj);
        end
      
      end
    end

  else 

    kpars   = dispmesh.points(:, 1);
    lambdas = dispmesh.points(:, 2);

    for idI = 1:numD
      fprintf('%d sur %d\n', idI, numD);

      kpar   = kpars(idI);
      lambda = lambdas(idI); 

      try

        % ========= (+) half-guide
        % FE matrix
        AApos = mat_pos_gradu_gradv...
              + 1i * kpar     * mat_pos_vectu_gradv...
              - 1i * kpar     * mat_pos_gradu_vectv...
              + (kpar * kpar) * mat_pos_vectu_vectv...
              + mat_pos_funQ_u_v...
              - lambda * mat_pos_funR_u_v;

        % Boundary condition
        BCstruct_pos = init_BCstruct_pos;
        BCstruct_pos.BCu  = 1i*lambda;
        BCstruct_pos.BCdu = -1;
        %
        % BCstruct_pos.BCu  = 1;
        % BCstruct_pos.BCdu = 0;

        % Half-guide problem
        [~, ~, ~, ~, ~, BCstruct_pos, boundary_op_pos] = PeriodicHalfGuideBVP(mesh_pos, +1, 2, AApos, BCstruct_pos, 0, opts);

        % ========= (-) half-guide
        % FE matrix
        AAneg = mat_neg_gradu_gradv...
              + 1i * kpar     * mat_neg_vectu_gradv...
              - 1i * kpar     * mat_neg_gradu_vectv...
              + (kpar * kpar) * mat_neg_vectu_vectv...
              + mat_neg_funQ_u_v...
              - lambda * mat_neg_funR_u_v;

        % Boundary condition
        BCstruct_neg = init_BCstruct_neg;
        BCstruct_neg.BCu  = -1i*lambda;
        BCstruct_neg.BCdu = 1;
        %
        % BCstruct_neg.BCu  = 1;
        % BCstruct_neg.BCdu = 0;

        % Half-guide problem
        [~, ~, ~, ~, ~, BCstruct_neg, boundary_op_neg] = PeriodicHalfGuideBVP(mesh_neg, -1, 2, AAneg, BCstruct_neg, 0, opts);

        % ========= Interior problem
        AAint = mat_int_gradu_gradv +...
              +  1i * kpar    * mat_int_vectu_gradv...
              -  1i * kpar    * mat_int_gradu_vectv...
              + (kpar * kpar) * mat_int_vectu_vectv...
              + mat_int_funQ_u_v...
              - lambda * mat_int_funR_u_v;

        % Coefficients associated to boundary condition
        Nbpos = size(boundary_op_pos, 1);
        Nbneg = size(boundary_op_neg, 1);

        BCu_int_pos = (BCstruct_pos.BCu' * eye(Nbpos) + BCstruct_pos.BCdu * boundary_op_pos)...
                    \ (boundary_op_pos * BCstruct_pos.BCu - BCstruct_pos.BCdu' * eye(Nbpos));
        BCu_int_neg = (BCstruct_neg.BCu' * eye(Nbneg) + BCstruct_neg.BCdu * boundary_op_neg)...
                    \ (boundary_op_neg * BCstruct_neg.BCu - BCstruct_neg.BCdu' * eye(Nbneg));

        SSpos = FEPack.pdes.Form.intg_TU_V(edgeYmaxInt, BCu_int_pos, BCstruct_pos.spB0, 'projection');
        SSneg = FEPack.pdes.Form.intg_TU_V(edgeYminInt, BCu_int_neg, BCstruct_neg.spB0, 'projection');

        AAint0 = PPint * (AAint - SSpos - SSneg) * PPint.';
        
        % Compute infinity norm of resolvent
        val(idI) = max(max(abs(AAint0 \ eye(size(AAint0)))));

      catch err_obj
        
        % val(idI) = 0;
        if (length(err_obj.message) >= num_err) &&...
          (strcmp(err_obj.message(1:num_err), err_prefix))
          % Some eigenvalues of the Riccati operator are on the unit circle
          val(idI) = 0;
        else
          % Unknown error
          rethrow(err_obj);
        end
      
      end
    end 

  end

end