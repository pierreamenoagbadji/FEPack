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
  shearmap.fun1y = @(s) spBY.evaluateBasisFunctions([mod(cutslope + s * ones(N1y-2, 1), 1), meshXY.points(edge1y.IdPoints(2:end-1), 1), zeros(N1y-2, 1)]);

  % shearmap2.fun0y = @(s) spBY.evaluateBasisFunctions([mod(           s * ones(N0y, 1), 1), meshXY.points(edge0y.IdPoints, 1), zeros(N0y, 1)]);
  % shearmap2.fun1y = @(s) spBY.evaluateBasisFunctions([mod(cutslope + s * ones(N0y, 1), 1), meshXY.points(edge1y.IdPoints, 1), zeros(N1y, 1)]);

  % Quadrature for DtN operators
  Nquad = N_s * 8;
  quadpoints = linspace(0, 1, Nquad + 1).';
  quadpoints(end) = [];
  structLoc = meshLineZ.domain('volumic').locateInDomain([quadpoints, zeros(Nquad, 2)]);
  elts = meshLineZ.domain('volumic').elements(structLoc.elements, :);
  coos = structLoc.barycoos;
  quadwgt = quadpoints(2) - quadpoints(1);
  
  for idFB = 1:opts.numFloquetPoints
    
    %% Local cell problems and DtN operators
    fprintf('%d sur %d\n', idFB, opts.numFloquetPoints);
    fprintf('<strong>Problèmes de cellule locaux et opérateurs DtN.</strong>\n');
    
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

      % DtN operators
      for idI = 1:4
        nameEsolI = ['E', sidenames{idI}];
        AA_m_solI = AA * solcell.(nameEsolI);

        for idJ = 1:4
          nameEsolJ = ['E', sidenames{idJ}];
          nameEdgeT = ['edgeT', sidenames{idI}, sidenames{idJ}];
          solcell.(nameEdgeT) = solcell.(nameEsolJ)' * AA_m_solI;
        end
      end

      % Save local cell solutions
      parsave([nomdossier, 'local_cell_sol_', suffix, '_Floquet_', num2str(idFB), '_S_', num2str(idS)], solcell, true);
    end
    system(['cp ', nomdossier, 'local_cell_sol_', suffix, '_Floquet_', num2str(idFB), '_S_1.mat ', nomdossier, 'local_cell_sol_', suffix, '_Floquet_', num2str(idFB), '_S_', num2str(N_s), '.mat']);
    
    % Initialize the Face DtN operators
    T0x0x = zeros(NbX, NbX); T0x1x = T0x0x; T1x0x = T0x0x; T1x1x = T0x0x; 
    T0x0y = zeros(NbY, NbX); T0x1y = T0x0y; T1x0y = T0x0y; T1x1y = T0x0y; 
    T0y0x = zeros(NbX, NbY); T0y1x = T0y0x; T1y0x = T0y0x; T1y1x = T0y0x; 
    T0y0y = zeros(NbY, NbY); T0y1y = T0y0y; T1y0y = T0y0y; T1y1y = T0y0y;
    % M0x0x = zeros(NbX, NbX); M1x1x = zeros(NbX, NbX);
    % M0y0y = zeros(NbY, NbY); M1y1y = zeros(NbY, NbY);
    % M_edge_0x0x = FEPack.pdes.Form.intg(meshXY.domain('xmin'), u * v); M_edge_0x0x = M_edge_0x0x(meshXY.domain('xmin').IdPoints, meshXY.domain('xmin').IdPoints);
    % M_edge_1x1x = FEPack.pdes.Form.intg(meshXY.domain('xmax'), u * v); M_edge_1x1x = M_edge_1x1x(meshXY.domain('xmax').IdPoints, meshXY.domain('xmax').IdPoints);
    % M_edge_0y0y = FEPack.pdes.Form.intg(meshXY.domain('ymin'), u * v); M_edge_0y0y = M_edge_0y0y(meshXY.domain('ymin').IdPoints, meshXY.domain('ymin').IdPoints);
    % M_edge_1y1y = FEPack.pdes.Form.intg(meshXY.domain('ymax'), u * v); M_edge_1y1y = M_edge_1y1y(meshXY.domain('ymax').IdPoints, meshXY.domain('ymax').IdPoints);
    
    parfor idQuad = 1:Nquad

      %% Local auxiliary DtN operators
      % basisfun = [];
      sVar = quadpoints(idQuad);
      solcell1 = load([nomdossier, 'local_cell_sol_', suffix, '_Floquet_', num2str(idFB), '_S_', num2str(elts(idQuad, 1))]);
      solcell2 = load([nomdossier, 'local_cell_sol_', suffix, '_Floquet_', num2str(idFB), '_S_', num2str(elts(idQuad, 2))]);

      % % tic;
      % A_m_fun0x = AA * basisfun.fun0x; 
      % A_m_fun1x = AA * basisfun.fun1x; 
      % A_m_fun0y = AA * basisfun.fun0y; 
      % A_m_fun1y = AA * basisfun.fun1y; 

      % T0x0x = T0x0x + shearmap.fun0x(sVar)' 

      T0x0x = T0x0x + shearmap.fun0x(sVar)' * (coos(idQuad, 1) * solcell1.edgeT0x0x + coos(idQuad, 2) * solcell2.edgeT0x0x) * shearmap.fun0x(sVar); 
      T1x0x = T1x0x + shearmap.fun0x(sVar)' * (coos(idQuad, 1) * solcell1.edgeT1x0x + coos(idQuad, 2) * solcell2.edgeT1x0x) * shearmap.fun1x(sVar); 
      T0x1x = T0x1x + shearmap.fun1x(sVar)' * (coos(idQuad, 1) * solcell1.edgeT0x1x + coos(idQuad, 2) * solcell2.edgeT0x1x) * shearmap.fun0x(sVar); 
      T1x1x = T1x1x + shearmap.fun1x(sVar)' * (coos(idQuad, 1) * solcell1.edgeT1x1x + coos(idQuad, 2) * solcell2.edgeT1x1x) * shearmap.fun1x(sVar);
      
      T0x0y = T0x0y + shearmap.fun0y(sVar)' * (coos(idQuad, 1) * solcell1.edgeT0x0y + coos(idQuad, 2) * solcell2.edgeT0x0y) * shearmap.fun0x(sVar); 
      T1x0y = T1x0y + shearmap.fun0y(sVar)' * (coos(idQuad, 1) * solcell1.edgeT1x0y + coos(idQuad, 2) * solcell2.edgeT1x0y) * shearmap.fun1x(sVar); 
      T0x1y = T0x1y + shearmap.fun1y(sVar)' * (coos(idQuad, 1) * solcell1.edgeT0x1y + coos(idQuad, 2) * solcell2.edgeT0x1y) * shearmap.fun0x(sVar); 
      T1x1y = T1x1y + shearmap.fun1y(sVar)' * (coos(idQuad, 1) * solcell1.edgeT1x1y + coos(idQuad, 2) * solcell2.edgeT1x1y) * shearmap.fun1x(sVar);
      
      T0y0x = T0y0x + shearmap.fun0x(sVar)' * (coos(idQuad, 1) * solcell1.edgeT0y0x + coos(idQuad, 2) * solcell2.edgeT0y0x) * shearmap.fun0y(sVar); 
      T1y0x = T1y0x + shearmap.fun0x(sVar)' * (coos(idQuad, 1) * solcell1.edgeT1y0x + coos(idQuad, 2) * solcell2.edgeT1y0x) * shearmap.fun1y(sVar); 
      T0y1x = T0y1x + shearmap.fun1x(sVar)' * (coos(idQuad, 1) * solcell1.edgeT0y1x + coos(idQuad, 2) * solcell2.edgeT0y1x) * shearmap.fun0y(sVar); 
      T1y1x = T1y1x + shearmap.fun1x(sVar)' * (coos(idQuad, 1) * solcell1.edgeT1y1x + coos(idQuad, 2) * solcell2.edgeT1y1x) * shearmap.fun1y(sVar);
      
      T0y0y = T0y0y + shearmap.fun0y(sVar)' * (coos(idQuad, 1) * solcell1.edgeT0y0y + coos(idQuad, 2) * solcell2.edgeT0y0y) * shearmap.fun0y(sVar); 
      T1y0y = T1y0y + shearmap.fun0y(sVar)' * (coos(idQuad, 1) * solcell1.edgeT1y0y + coos(idQuad, 2) * solcell2.edgeT1y0y) * shearmap.fun1y(sVar); 
      T0y1y = T0y1y + shearmap.fun1y(sVar)' * (coos(idQuad, 1) * solcell1.edgeT0y1y + coos(idQuad, 2) * solcell2.edgeT0y1y) * shearmap.fun0y(sVar); 
      T1y1y = T1y1y + shearmap.fun1y(sVar)' * (coos(idQuad, 1) * solcell1.edgeT1y1y + coos(idQuad, 2) * solcell2.edgeT1y1y) * shearmap.fun1y(sVar);

      % M0x0x = M0x0x + shearmap.fun0x(sVar)' * M_edge_0x0x * shearmap.fun0x(sVar); 
      % M1x1x = M1x1x + shearmap.fun0x(sVar)' * M_edge_1x1x * shearmap.fun1x(sVar);

      % M0y0y = M0y0y + shearmap2.fun0y(sVar)' * M_edge_0y0y * shearmap2.fun0y(sVar); 
      % M1y1y = M1y1y + shearmap2.fun1y(sVar)' * M_edge_1y1y * shearmap2.fun1y(sVar); 
      
    end

    % %% Compute DtD operator
    % fprintf('<strong>Calcul opérateur DtD.</strong>\n');

    % Divide DtN operators by the quadrature weight
    T0x0x = T0x0x * quadwgt / opts.cutvec(1);
    T1x0x = T1x0x * quadwgt / opts.cutvec(1);
    T0x1x = T0x1x * quadwgt / opts.cutvec(1);
    T1x1x = T1x1x * quadwgt / opts.cutvec(1);
    T0x0y = T0x0y * quadwgt / opts.cutvec(1);
    T1x0y = T1x0y * quadwgt / opts.cutvec(1);
    T0x1y = T0x1y * quadwgt / opts.cutvec(1);
    T1x1y = T1x1y * quadwgt / opts.cutvec(1);
    T0y0x = T0y0x * quadwgt / opts.cutvec(1);
    T1y0x = T1y0x * quadwgt / opts.cutvec(1);
    T0y1x = T0y1x * quadwgt / opts.cutvec(1);
    T1y1x = T1y1x * quadwgt / opts.cutvec(1);
    T0y0y = T0y0y * quadwgt / opts.cutvec(1);
    T1y0y = T1y0y * quadwgt / opts.cutvec(1);
    T0y1y = T0y1y * quadwgt / opts.cutvec(1);
    T1y1y = T1y1y * quadwgt / opts.cutvec(1);

    % M0x0x = spBX.invmassmat * M0x0x * quadwgt / opts.cutvec(1);
    % M0y0y = spBY.invmassmat * M0y0y * quadwgt / opts.cutvec(1);
    % M1x1x = spBX.invmassmat * M1x1x * quadwgt / opts.cutvec(1);
    % M1y1y = spBY.invmassmat * M1y1y * quadwgt / opts.cutvec(1);

    % fprintf('M0x0x: %d et %d\n', norm(M0x0x - eye(size(M0x0x))), norm(M0x0x - diag(diag(M0x0x))));
    % fprintf('M0y0y: %d et %d\n', norm(M0y0y - eye(size(M0y0y))), norm(M0y0y - diag(diag(M0y0y))));
    % fprintf('M1x1x: %d et %d\n', norm(M1x1x - eye(size(M1x1x))), norm(M1x1x - diag(diag(M1x1x))));
    % fprintf('M1y1y: %d et %d\n', norm(M1y1y - eye(size(M1y1y))), norm(M1y1y - diag(diag(M1y1y))));
    % % error('Stop');
    % % disp(real(M0x0x));

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
    
    %% Cas homogène
    % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Cas homogène
    % FourierIdsX = spBX.FourierIds.X(:).'; dx = length(FourierIdsX);
    % FourierIdsY = spBX.FourierIds.Y(:).'; dy = length(FourierIdsY);
    % FourierIdsZ = spBX.FourierIds.Z(:).'; dz = length(FourierIdsZ);

    % [Ix, Iy, ~] = ind2sub([dx, dy, dz], (1:spBX.numBasis)');
    % r_xi_fun = FloquetVar*opts.cutvec(1) + 2*pi*(FourierIdsX(Ix) * opts.cutvec(1) + FourierIdsY(Iy) * opts.cutvec(2)).';
    % l_xi_fun = sqrt(r_xi_fun.^2 - opts.omega^2);

    % mZXpts = spBY.domain.mesh.points;

    % solguide.DtD0 = spBY.FE_to_spectral * (sinh((1 - orientation*mZXpts(:, 2)) * l_xi_fun.') .* exp(2i*pi*mZXpts(:, 1) * FourierIdsY(Iy)) ./ (ones(size(mZXpts, 1), 1) * sinh(l_xi_fun.')));
    % %
    % % norm(solguide.DtD0)
    % %
    % solguide.DtD1 = spBY.FE_to_spectral * (sinh(orientation*mZXpts(:, 2) * l_xi_fun.') .* exp(2i*pi*mZXpts(:, 1) * FourierIdsY(Iy)) ./ (ones(size(mZXpts, 1), 1) * sinh(l_xi_fun.')));


    % F00 = diag(l_xi_fun./tanh(l_xi_fun));
    % F01 = -diag(l_xi_fun./sinh(l_xi_fun));
    % F10 = F01;
    % F11 = F00;

    % solguide.R = diag(exp(-l_xi_fun));
    % solguide.D = solguide.R;
    % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
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
