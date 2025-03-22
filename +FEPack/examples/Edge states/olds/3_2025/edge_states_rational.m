function [kpars, Egies, val] = edge_states_rational(parallel, HcObj, dirac, argspde, taskset)

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  %% Preliminaries
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  delta = argspde.delta;

  % Make sure the edge coefficients are coprime integers
  a1 = HcObj.edge.a1;
  b1 = HcObj.edge.b1;
  [G, x, y] = gcd(a1, b1);  % a1*x + b1*y = G
  a2 = -y; b2 = x;

  if (G ~= 1)
    error('a1 et b1 doivent être premiers entre eux.');
  end

  % Edge vectors
  edge_vec1      =  a1 * HcObj.vecPer1  + b1 * HcObj.vecPer2;
  edge_vec2      =  a2 * HcObj.vecPer1  + b2 * HcObj.vecPer2;
  edge_dual_vec1 =  b2 * HcObj.dualVec1 - a2 * HcObj.dualVec2;
  edge_dual_vec2 = -b1 * HcObj.dualVec1 + a1 * HcObj.dualVec2;


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  %% Range around Dirac point
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  [~, minEgy, maxEgy] = compute_dirac_range(...
                              dirac.numEigvals,...
                              HcObj.vecPer1, HcObj.vecPer2, HcObj.dualVec1, HcObj.dualVec2,...
                              HcObj.V, HcObj.W, delta,...
                              edge_vec1, edge_vec2, edge_dual_vec1, edge_dual_vec2,...
                              dirac.numNodesX, dirac.numNodesY, dirac.numNodesK, true,...
                              parallel.use, parallel.numcores, 1e-1, taskset.spectrum_pert_bulk.compute);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  %% Meshes
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %% (+) half-guide
  % ======= Mesh
  mesh_pos = FEPack.meshes.MeshRectangle(1, [0 1], [0  1], argspde.numNodesXpos, argspde.numNodesYpos);
  cell_pos = mesh_pos.domain('volumic');

  % ======= Boundary conditions
  if strcmpi(argspde.bc_basis_functions_pos, 'Lagrange')
    BCstruct_pos.spB0 = FEPack.spaces.PeriodicLagrangeBasis(mesh_pos.domains{4});
    BCstruct_pos.spB1 = FEPack.spaces.PeriodicLagrangeBasis(mesh_pos.domains{3});
  else
    FourierIds = [0 0]; FourierIds(1) = argspde.numNodesXpos/4;
    BCstruct_pos.spB0 = FEPack.spaces.FourierBasis(mesh_pos.domains{4}, FourierIds);
    BCstruct_pos.spB1 = FEPack.spaces.FourierBasis(mesh_pos.domains{3}, FourierIds);
  end

  %% (-) half-guide
  % ======= Mesh
  mesh_neg = FEPack.meshes.MeshRectangle(1, [0 1], [0 -1], argspde.numNodesXneg, argspde.numNodesYneg);
  cell_neg = mesh_neg.domain('volumic');
  
  % ======= Boundary conditions
  if strcmpi(argspde.bc_basis_functions_neg, 'Lagrange')
    BCstruct_neg.spB0 = FEPack.spaces.PeriodicLagrangeBasis(mesh_neg.domains{4});
    BCstruct_neg.spB1 = FEPack.spaces.PeriodicLagrangeBasis(mesh_neg.domains{3});
  else
    FourierIds = [0 0]; FourierIds(1) = argspde.numNodesXneg/4;
    BCstruct_neg.spB0 = FEPack.spaces.FourierBasis(mesh_neg.domains{4}, FourierIds);
    BCstruct_neg.spB1 = FEPack.spaces.FourierBasis(mesh_neg.domains{3}, FourierIds);
  end

  %% Interior domain
  yInt = ceil(1/(2*pi*argspde.delta));
  mesh_int = FEPack.meshes.MeshRectangle(1, [0 1], [-yInt yInt], argspde.numNodesXint, argspde.numNodesYint);
  cell_int  = mesh_int.domain('volumic');
  edge0xInt = mesh_int.domain('xmin');
  edge1xInt = mesh_int.domain('xmax');
  edgeYminInt = mesh_int.domain('ymin');
  edgeYmaxInt = mesh_int.domain('ymax');

  %% Mesh the parallel quasi-momentum range
  if (argspde.center_parallel_quasi_momentum_around_K)
    kpar_center = HcObj.highSymK.' * edge_vec1;
  else
    kpar_center = 0;
  end 
  numKpar = argspde.numKpar;
  kpars = kpar_center + linspace(argspde.minKpar, argspde.maxKpar, numKpar);

  %% Mesh the energy range
  % kappa_plus_inf  = HcObj.kappa(+inf);
  % kappa_minus_inf = HcObj.kappa(-inf);

  if (~isfield(argspde, 'minEgy') ||...
      (isfield(argspde, 'minEgy') && isempty(argspde.minEgy)) ||...
      ~isfield(argspde, 'maxEgy') ||...
      (isfield(argspde, 'maxEgy') && isempty(argspde.maxEgy)))

    % If no energy range has been specified, use the one predicted by theory
    argspde.minEgy = minEgy + 0.1*(maxEgy - minEgy);
    argspde.maxEgy = maxEgy - 0.1*(maxEgy - minEgy);

  end
  numEgy  = argspde.numEgy;
  Egies = linspace(argspde.minEgy, argspde.maxEgy, numEgy);

  opts.tol = 1e-2 * (maxEgy - minEgy);
  opts.computeSol = false;
  opts.solBasis = false;
  opts.verbose = 0;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  %% Coefficients
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  Rmat = [edge_vec1, edge_vec2];
  Tmat = [edge_dual_vec1'; edge_dual_vec2'] / (2*pi); % Inverse of Rmat
  Te1 = Tmat' * [1; 0];

  pot = @(x) HcObj.V(x) + argspde.delta * HcObj.kappa(argspde.delta * x(:, 1:2) * edge_dual_vec2) .* HcObj.W(x);
  trPot = @(x) pot((Rmat * x(:, 1:2)')');

  Vpot = @(x) HcObj.V((Rmat * x(:, 1:2)')');
  Wpot = @(x) HcObj.W((Rmat * x(:, 1:2)')');

  if (taskset.plot_coefficients)
    % Plot coefficient
    numX = 2;
    numY = 2;
    figure;
    set(groot,'defaultAxesTickLabelInterpreter','latex');
    set(groot,'defaulttextinterpreter','latex');
    set(groot,'defaultLegendInterpreter','latex');
    for idX = -numX:numX-1
      for idY = -numY:numY-1
        X = mesh_pos.points(:, 1) + idX;
        Y = mesh_pos.points(:, 2) + idY;

        % V
        subplot(1, 2, 1);
        trisurf(mesh_pos.triangles, X, Y, HcObj.V([X, Y]));
        hold on;
        shading interp;
        view(2);
        set(gca, 'DataAspectRatio', [1 1 1], 'FontSize', 16);
        colorbar('TickLabelInterpreter', 'latex');
        xlabel('$V$');

        % W
        subplot(1, 2, 2);
        trisurf(mesh_pos.triangles, X, Y, HcObj.W([X, Y]));
        hold on;
        shading interp;
        view(2);
        set(gca, 'DataAspectRatio', [1 1 1], 'FontSize', 16);
        colorbar('TickLabelInterpreter', 'latex');
        xlabel('$W$');
      end
    end

    figure;
    set(groot,'defaultAxesTickLabelInterpreter','latex');
    set(groot,'defaulttextinterpreter','latex');
    set(groot,'defaultLegendInterpreter','latex');
    for idX = -2*numX:2*numX-1
      for idY = -2*numY:2*numY-1
        X = mesh_pos.points(:, 1) + idX;
        Y = mesh_pos.points(:, 2) + idY;

        % Q
        trisurf(mesh_pos.triangles, X, Y, pot([X, Y]));
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
  %
  %% Finite element matrices
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  u = FEPack.pdes.PDEObject; 
  v = dual(u);

  % (+) half-guide
  fprintf('(+) half-guide FE matrices computation\n');
  mat_pos_gradu_gradv = FEPack.pdes.Form.intg(cell_pos, (Tmat' * grad2(u))   * (Tmat' * grad2(v)));
  mat_pos_gradu_vec1v = FEPack.pdes.Form.intg(cell_pos, (Tmat' * grad2(u))   * (Te1   * v));
  mat_pos_vec1u_gradv = FEPack.pdes.Form.intg(cell_pos, (Te1   * u) * (Tmat' * grad2(v)));
  mat_pos_vec1u_vec1v = FEPack.pdes.Form.intg(cell_pos, (Te1   * u) * (Te1   * v));
  mat_pos_Vpot_u_v    = FEPack.pdes.Form.intg(cell_pos, (Vpot  * u) * v);
  mat_pos_Wpot_u_v    = FEPack.pdes.Form.intg(cell_pos, (Wpot  * u) * v);
  mat_pos_u_v = FEPack.pdes.Form.intg(cell_pos, u * v);
  fprintf('(+) half-guide FE matrices computation done.\n');

  % (-) half-guide
  fprintf('(-) half-guide FE matrices computation\n');
  mat_neg_gradu_gradv = FEPack.pdes.Form.intg(cell_neg, (Tmat' * grad2(u))   * (Tmat' * grad2(v)));
  mat_neg_gradu_vec1v = FEPack.pdes.Form.intg(cell_neg, (Tmat' * grad2(u))   * (Te1   * v));
  mat_neg_vec1u_gradv = FEPack.pdes.Form.intg(cell_neg, (Te1   * u) * (Tmat' * grad2(v)));
  mat_neg_vec1u_vec1v = FEPack.pdes.Form.intg(cell_neg, (Te1   * u) * (Te1   * v));
  mat_neg_Vpot_u_v    = FEPack.pdes.Form.intg(cell_neg, (Vpot  * u) * v);
  mat_neg_Wpot_u_v    = FEPack.pdes.Form.intg(cell_neg, (Wpot  * u) * v);
  mat_neg_u_v = FEPack.pdes.Form.intg(cell_neg, u * v);
  fprintf('(-) half-guide FE matrices computation done.\n');

  % Interior domain
  fprintf('Interior domain FE matrices computation\n');
  mat_int_gradu_gradv = FEPack.pdes.Form.intg(cell_int, (Tmat' * grad2(u))   * (Tmat' * grad2(v)));
  mat_int_gradu_vec1v = FEPack.pdes.Form.intg(cell_int, (Tmat' * grad2(u))   * (Te1   * v));
  mat_int_vec1u_gradv = FEPack.pdes.Form.intg(cell_int, (Te1   * u) * (Tmat' * grad2(v)));
  mat_int_vec1u_vec1v = FEPack.pdes.Form.intg(cell_int, (Te1   * u) * (Te1   * v));
  mat_int_trPot_u_v   = FEPack.pdes.Form.intg(cell_int, (trPot * u) * v);
  mat_int_u_v = FEPack.pdes.Form.intg(cell_int, u * v);
  fprintf('Interior domain FE matrices computation done.\n');

  ecs_int = (((u|edge0xInt) - (u|edge1xInt)) == 0.0);
  ecs_int.applyEcs;
  PPint = ecs_int.P;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  %% Edge states
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  val = zeros(numKpar, numEgy);

  if (parallel.use)

    if ~isfield(parallel, 'numcores')
      error('parallel.use has been set to true; you have to specify the number of cores through parallel.numcores. Use the command feature(''numcores'') to see the number of available cores.');
    end

    % % Parallel version
    % pool = gcp('nocreate'); % Get current parallel pool
    % if isempty(pool)
    %   parpool('local', numCores); % Start a parallel pool if none exists
    % end
    % fprintf('Using parfor (parallel execution)\n');


  else

    for idK = 1:numKpar
      
      % Current parallel quasi-momentum
      kx = kpars(idK);

      for idE = 1:numEgy

        fprintf('(E) %d sur %d ; (K) %d sur %d\n', idE, numEgy, idK, numKpar);

        % Current energy
        Egy = Egies(idE)
        opts.omega = sqrt(Egy);
        
        % (+) half-guide
        % ======= FE matrix
        AApos = mat_pos_gradu_gradv...
              +  1i * kx  * mat_pos_vec1u_gradv...
              -  1i * kx  * mat_pos_gradu_vec1v...
              + (kx * kx) * mat_pos_vec1u_vec1v...
              +         mat_pos_Vpot_u_v...
              + delta * mat_pos_Wpot_u_v...
              - Egy   * mat_pos_u_v;

        % ======= Robin boundary condition
        % BCstruct_pos.BCu = 1i*Egy;
        % BCstruct_pos.BCdu = -1;
        BCstruct_pos.BCu = 1;
        BCstruct_pos.BCdu = 0;

        % ======= Half-guide problem
        [~, ~, ~, ~, ~, BCstruct_pos, Lambda_pos] = PeriodicHalfGuideBVP(mesh_pos, +1, 2, AApos, BCstruct_pos, 1, opts);

        % (-) half-guide
        % ======= FE matrix
        AAneg = mat_neg_gradu_gradv...
              +  1i * kx  * mat_neg_vec1u_gradv...
              -  1i * kx  * mat_neg_gradu_vec1v...
              + (kx * kx) * mat_neg_vec1u_vec1v...
              +         mat_neg_Vpot_u_v...
              - delta * mat_neg_Wpot_u_v...
              - Egy   * mat_neg_u_v;

        % ======= Robin boundary condition
        BCstruct_neg.BCu = -1i*Egy;
        BCstruct_neg.BCdu = 1;
        % BCstruct_neg.BCu = 1;
        % BCstruct_neg.BCdu = 0;

        % ======= Half-guide problem
        [~, ~, ~, ~, ~, BCstruct_neg, Lambda_neg] = PeriodicHalfGuideBVP(mesh_neg, -1, 2, AAneg, BCstruct_neg, 1, opts);

        % Interior problem
        AAint = mat_int_gradu_gradv +...
              +  1i * kx  * mat_int_vec1u_gradv...
              -  1i * kx  * mat_int_gradu_vec1v...
              + (kx * kx) * mat_int_vec1u_vec1v...
              + mat_int_trPot_u_v...
              - Egy * mat_int_u_v;

        % ======= Coefficients associated to boundary condition
        Nbpos = size(Lambda_pos, 1);
        Nbneg = size(Lambda_neg, 1);

        BCu_int_pos = (BCstruct_pos.BCu' * eye(Nbpos) + BCstruct_pos.BCdu * Lambda_pos)...
                    \ (Lambda_pos * BCstruct_pos.BCu - BCstruct_pos.BCdu' * eye(Nbpos));
        BCu_int_neg = (BCstruct_neg.BCu' * eye(Nbneg) + BCstruct_neg.BCdu * Lambda_neg)...
                    \ (Lambda_neg * BCstruct_neg.BCu - BCstruct_neg.BCdu' * eye(Nbneg));

        SSpos = FEPack.pdes.Form.intg_TU_V(edgeYmaxInt, BCu_int_pos, BCstruct_pos.spB0, 'projection');
        SSneg = FEPack.pdes.Form.intg_TU_V(edgeYminInt, BCu_int_neg, BCstruct_neg.spB0, 'projection');

        AAint0 = PPint * (AAint - SSpos - SSneg) * PPint.';
        
        val(idK, idE) = max(max(abs(AAint0 \ eye(size(AAint0)))));
        % val(idK, idE) = det(AAint0);
        % val(idK, idE) = eigs(AAint0, 1, 'smallestabs');

        % error();
      end

    end

  end 

end