function quasi2DHalfGuideDirichlet_par(orientation, meshXY, meshLineZ, mu3D, rho3D, BCstruct, opts)

  %% Initialization
  cutslope = opts.cutvec(2) / opts.cutvec(1);
  DeltaS = meshLineZ.points(2, 1) - meshLineZ.points(1, 1);
  N_s = meshLineZ.numPoints;
  N_XY = meshXY.numPoints;
  nomdossier = opts.nomdossier;

  % Basis functions on faces 
  spBX = BCstruct.spBX; NbX = spBX.numBasis; % X = cst
  spBY = BCstruct.spBY; NbY = spBY.numBasis; % Y = cst

  % Integrands for FE matrices
  u = FEPack.pdes.PDEObject; v = dual(u);
  gradu_gradv = @(muco) (muco * grad2(u)) * grad2(v);
  gradu_vecYv = @(muco) (muco * grad2(u)) * ([0; opts.cutvec(1)] * v);
  vecYu_gradv = @(muco) (muco * ([0; opts.cutvec(1)] * u)) * grad2(v);
  vecYu_vecYv = @(muco) (muco * ([0; opts.cutvec(1)] * u)) * ([0; opts.cutvec(1)] * v);
  u_v = @(rhoco) ((rhoco*u)*v);

  % Edges
  edge0x = meshXY.domain('xmin'); N0x = edge0x.numPoints;
  edge1x = meshXY.domain('xmax'); N1x = edge1x.numPoints;
  edge0y = meshXY.domain('ymin'); N0y = edge0y.numPoints;
  edge1y = meshXY.domain('ymax'); N1y = edge1y.numPoints;

  %% FE matrices
  sVec = meshLineZ.points;
  cellXY = meshXY.domain('volumic');
  cutvec = opts.cutvec;
  suffix = opts.suffix;
  omega = opts.omega;

  % parpool
  parfor idS = 1:N_s-1
    fprintf('Calcul des matrices éléments finis.\n');
    sVar = sVec(idS, 1);
    mu2Ds  = @(x)  mu3D([x(:, 1), cutvec(1)*x(:, 2), cutvec(2)*x(:, 2) + sVar]); %#ok
    rho2Ds = @(x) rho3D([x(:, 1), cutvec(1)*x(:, 2), cutvec(2)*x(:, 2) + sVar]); %#ok

    FEmat = [];
    tic;
    FEmat.mat_gradu_gradv = FEPack.pdes.Form.intg(cellXY, gradu_gradv(mu2Ds));
    FEmat.mat_gradu_vecYv = FEPack.pdes.Form.intg(cellXY, gradu_vecYv(mu2Ds));
    FEmat.mat_vecYu_gradv = FEPack.pdes.Form.intg(cellXY, vecYu_gradv(mu2Ds));
    FEmat.mat_vecYu_vecYv = FEPack.pdes.Form.intg(cellXY, vecYu_vecYv(mu2Ds));
    FEmat.mat_u_v         = FEPack.pdes.Form.intg(cellXY,        u_v(rho2Ds));
    toc;

    parsave([nomdossier, 'FEmat_', suffix, '_', num2str(idS)], FEmat, true);
  end
  % delete(gcp('nocreate'))

  %% Compute Floquet-Bloch transform
  % Homogeneous Dirichlet boundary conditions
  ecs = ((u|edge0x) == 0.0) & ((u|edge1x) == 0.0) &...
        ((u|edge0y) == 0.0) & ((u|edge1y) == 0.0);

  % Surfacic rhs
  B0x = sparse(edge0x.IdPoints(1:end), (1:N0x), 1.0, N_XY, N0x);
  B1x = sparse(edge1x.IdPoints(1:end), (1:N1x), 1.0, N_XY, N1x);
  B0y = sparse(edge0y.IdPoints(2:end-1), (1:N0y-2), 1.0, N_XY, N0y-2);
  B1y = sparse(edge1y.IdPoints(2:end-1), (1:N1y-2), 1.0, N_XY, N1y-2);

  ecs.applyEcs;
  ecs.b = [B0x, B1x, B0y, B1y];
  sidenames = {'0x', '1x', '0y', '1y'};

  % Define shear map
  shearmap.fun0x = @(s) spBX.evaluateBasisFunctions([    opts.cutvec(1) * meshXY.points(edge0x.IdPoints, 2),...
                                                     mod(opts.cutvec(2) * meshXY.points(edge0x.IdPoints, 2) + s, 1), zeros(N0x, 1)]);
  %
  shearmap.fun1x = @(s) spBX.evaluateBasisFunctions([    opts.cutvec(1) * meshXY.points(edge1x.IdPoints, 2),...
                                                     mod(opts.cutvec(2) * meshXY.points(edge1x.IdPoints, 2) + s, 1), zeros(N1x, 1)]);
  %
  shearmap.fun0y = @(s) spBY.evaluateBasisFunctions([mod(           s * ones(N0y-2, 1), 1), meshXY.points(edge0y.IdPoints(2:end-1), 1), zeros(N0y-2, 1)]);
  shearmap.fun1y = @(s) spBY.evaluateBasisFunctions([mod(cutslope + s * ones(N0y-2, 1), 1), meshXY.points(edge1y.IdPoints(2:end-1), 1), zeros(N1y-2, 1)]);
  
  
  for idFB = 1:opts.numFloquetPoints
    
    %% Local cell problems and DtN operators
    fprintf('%d sur %d\n', idFB, opts.numFloquetPoints);
    fprintf('<strong>Problèmes de cellule locaux et opérateurs DtN.</strong>\n');
    
    % Initialize the Face DtN operators
    T0x0x = zeros(NbX, NbX); T0x1x = T0x0x; T1x0x = T0x0x; T1x1x = T0x0x; 
    T0x0y = zeros(NbY, NbX); T0x1y = T0x0y; T1x0y = T0x0y; T1x1y = T0x0y; 
    T0y0x = zeros(NbX, NbY); T0y1x = T0y0x; T1y0x = T0y0x; T1y1x = T0y0x; 
    T0y0y = zeros(NbY, NbY); T0y1y = T0y0y; T1y0y = T0y0y; T1y1y = T0y0y;
    
    % Floquet variable
    FloquetVar = opts.FloquetPoints(idFB);
    
    parfor idS = 1:N_s-1
      
      %% Local cell problems
      % The FE matrix is a linear combination of elementary matrices
      FEmat = load([nomdossier, 'FEmat_', suffix, '_', num2str(idS), '.mat']);

      AA = 1i * FloquetVar * FEmat.mat_vecYu_gradv + FEmat.mat_gradu_gradv...
         - 1i * FloquetVar * FEmat.mat_gradu_vecYv + FloquetVar * FloquetVar * FEmat.mat_vecYu_vecYv...
         - (omega^2) * FEmat.mat_u_v;
      
      AA0 =  ecs.P * AA * ecs.P'; %#ok
      LL0 = -ecs.P * AA * ecs.b;
      
      % Solve the linear system
      % tic;
      Ecell0 = ecs.b + ecs.P' * (AA0 \ LL0);
      % tpscpu = toc;
      % fprintf('Inversion système : %0.3e secondes\n', tpscpu);

      % Local cell solutions
      solcell = [];
      solcell.E0x = Ecell0(:, 1:N0x);   Ecell0(:, 1:N0x)   = [];
      solcell.E1x = Ecell0(:, 1:N1x);   Ecell0(:, 1:N1x)   = [];
      solcell.E0y = Ecell0(:, 1:N0y-2); Ecell0(:, 1:N0y-2) = [];
      solcell.E1y = Ecell0(:, 1:N1y-2);

      % Save local cell solutions
      parsave([nomdossier, 'local_cell_sol_', suffix, '_Floquet_', num2str(idFB), '_S_', num2str(idS)], solcell, true);

      %% Local auxiliary DtN operators
      basisfun = [];
      sVar = sVec(idS, 1);
      
      for idI = 1:4
        nameE = ['E', sidenames{idI}]; %#ok
        nameFun = ['fun', sidenames{idI}];

        basisfun.(nameFun) = solcell.(nameE) * shearmap.(nameFun)(sVar); %#ok
      end

      % tic;
      A_m_fun0x = AA * basisfun.fun0x; 
      A_m_fun1x = AA * basisfun.fun1x; 
      A_m_fun0y = AA * basisfun.fun0y; 
      A_m_fun1y = AA * basisfun.fun1y; 

      T0x0x = T0x0x + basisfun.fun0x' * A_m_fun0x; 
      T1x0x = T1x0x + basisfun.fun0x' * A_m_fun1x; 
      T0x1x = T0x1x + basisfun.fun1x' * A_m_fun0x; 
      T1x1x = T1x1x + basisfun.fun1x' * A_m_fun1x;
      
      T0x0y = T0x0y + basisfun.fun0y' * A_m_fun0x; 
      T1x0y = T1x0y + basisfun.fun0y' * A_m_fun1x; 
      T0x1y = T0x1y + basisfun.fun1y' * A_m_fun0x; 
      T1x1y = T1x1y + basisfun.fun1y' * A_m_fun1x;
      
      T0y0x = T0y0x + basisfun.fun0x' * A_m_fun0y; 
      T1y0x = T1y0x + basisfun.fun0x' * A_m_fun1y; 
      T0y1x = T0y1x + basisfun.fun1x' * A_m_fun0y; 
      T1y1x = T1y1x + basisfun.fun1x' * A_m_fun1y;
      
      T0y0y = T0y0y + basisfun.fun0y' * A_m_fun0y; 
      T1y0y = T1y0y + basisfun.fun0y' * A_m_fun1y; 
      T0y1y = T0y1y + basisfun.fun1y' * A_m_fun0y; 
      T1y1y = T1y1y + basisfun.fun1y' * A_m_fun1y;
      
    end
    system(['cp ', nomdossier, 'local_cell_sol_', suffix, '_Floquet_', num2str(idFB), '_S_1.mat ', nomdossier, 'local_cell_sol_', suffix, '_Floquet_', num2str(idFB), '_S_', num2str(N_s), '.mat']);

    %% Compute DtD operator
    fprintf('<strong>Calcul opérateur DtD.</strong>\n');

    % Divide DtN operators by the quadrature weight
    T0x0x = T0x0x * DeltaS / opts.cutvec(1);
    T1x0x = T1x0x * DeltaS / opts.cutvec(1);
    T0x1x = T0x1x * DeltaS / opts.cutvec(1);
    T1x1x = T1x1x * DeltaS / opts.cutvec(1);
    T0x0y = T0x0y * DeltaS / opts.cutvec(1);
    T1x0y = T1x0y * DeltaS / opts.cutvec(1);
    T0x1y = T0x1y * DeltaS / opts.cutvec(1);
    T1x1y = T1x1y * DeltaS / opts.cutvec(1);
    T0y0x = T0y0x * DeltaS / opts.cutvec(1);
    T1y0x = T1y0x * DeltaS / opts.cutvec(1);
    T0y1x = T0y1x * DeltaS / opts.cutvec(1);
    T1y1x = T1y1x * DeltaS / opts.cutvec(1);
    T0y0y = T0y0y * DeltaS / opts.cutvec(1);
    T1y0y = T1y0y * DeltaS / opts.cutvec(1);
    T0y1y = T0y1y * DeltaS / opts.cutvec(1);
    T1y1y = T1y1y * DeltaS / opts.cutvec(1);

    % Solve linear system
    DtD = (T1y0y + T0y0y + T1y1y + T0y1y) \ [-(T0x0y + T0x1y), -(T1x0y + T1x1y)];
    solguide.DtD0 = DtD(:,     1:NbX);
    solguide.DtD1 = DtD(:, NbX+1:end);

    %% Traces and normal traces of the local cell solutions
    % Ekl is the trace of Ek on Sigmal
    E00 = speye(NbX);
    E10 = sparse(NbX, NbX);
    E01 = sparse(NbX, NbX);
    E11 = speye(NbX);

    % Fkl is the normal trace of Ek on Sigmal. It is evaluated weakly and then projected
    fprintf('<strong>Opérateurs DtN locaux.</strong>\n');
    F00 = spBX.invmassmat * ((T0y0x + T1y0x) * solguide.DtD0 + T0x0x);
    F10 = spBX.invmassmat * ((T0y0x + T1y0x) * solguide.DtD1 + T1x0x);
    F01 = spBX.invmassmat * ((T0y1x + T1y1x) * solguide.DtD0 + T0x1x);
    F11 = spBX.invmassmat * ((T0y1x + T1y1x) * solguide.DtD1 + T1x1x);
    
    %% Solve the Riccati equation
    fprintf('<strong>Equation de Riccati.</strong>\n');
    
    flux = @(V) modesFlux(V, E00, E10, F00, F10, orientation, spB0, opts.omega);
    riccatiOpts.tol = 1.0e-2;
    if isfield(opts, 'suffix')
      riccatiOpts.suffix = suffix;
    else
      riccatiOpts.suffix = '';
    end

    % Solve the linearized eigenvalue problem associated to the Riccati equation
    tic;
    [solguide.R, solguide.D] = propagationOperators([E01, E11;  orientation*F01,  orientation*F11], ...
                                                    [E00, E10; -orientation*F00, -orientation*F10], flux, riccatiOpts);
    tpscpu = toc;
    fprintf('Résolution équation de Riccati : %0.3e secondes\n', tpscpu);
    
    %% Compute the DtN operator
    U0 = E00 + E10 * solguide.D;
    dU0 = F00 + F10 * solguide.D;

    solguide.Lambda = -BCstruct.BCu' * dU0 + BCstruct.BCdu' * U0;

    %% Save outputs
    save([nomdossier, 'half_guide_sol_', suffix, '_Floquet_', num2str(idFB)], '-struct', 'solguide', '-v7.3');
  end

end

%% Auxiliary function: The flux function
function flux = modesFlux(V, E00, E10, F00, F10, orientation, spB, omega)
  
  Nb = size(E00, 1);
  Psi  = E00 * V(1:Nb, :) + E10 * V(Nb+1:end, :);
  dPsi = F00 * V(1:Nb, :) + F10 * V(Nb+1:end, :);

  flux = -orientation * imag(diag(Psi' * spB.massmat * dPsi) / omega);

end
