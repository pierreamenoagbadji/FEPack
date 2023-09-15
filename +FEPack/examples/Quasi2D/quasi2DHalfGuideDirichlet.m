function quasi2DHalfGuideDirichlet(orientation, meshYZ, meshZX, meshXY, omega, mu3D, rho3D, BCstruct, opts)

  %% Initialization
  cutslope = opts.cutvec(2) / opts.cutvec(1);
  sLine = meshYZ.domain('xmin');
  svec = meshYZ.points(sLine.IdPoints, 2);
  DeltaS = svec(2) - svec(1);
  N_s = numel(svec);
  N_XY = meshXY.numPoints;

  % Basis functions on faces 
  % faceX = meshYZ.domain('volumic'); % X = cst
  % faceY = meshZX.domain('volumic'); % Y = cst
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
  for idS = 1:N_s-1
    fprintf('Calcul des matrices éléments finis.\n');
    sVar = svec(idS);
    mu2Ds  = @(x)  mu3D([x(:, 1), opts.cutvec(1)*x(:, 2), opts.cutvec(2)*x(:, 2) + sVar]);
    rho2Ds = @(x) rho3D([x(:, 1), opts.cutvec(1)*x(:, 2), opts.cutvec(2)*x(:, 2) + sVar]);

    tic;
    FEmat.mat_gradu_gradv = FEPack.pdes.Form.intg(meshXY.domain('volumic'), gradu_gradv(mu2Ds));
    FEmat.mat_gradu_vecYv = FEPack.pdes.Form.intg(meshXY.domain('volumic'), gradu_vecYv(mu2Ds));
    FEmat.mat_vecYu_gradv = FEPack.pdes.Form.intg(meshXY.domain('volumic'), vecYu_gradv(mu2Ds));
    FEmat.mat_vecYu_vecYv = FEPack.pdes.Form.intg(meshXY.domain('volumic'), vecYu_vecYv(mu2Ds));
    FEmat.mat_u_v         = FEPack.pdes.Form.intg(meshXY.domain('volumic'),        u_v(rho2Ds));
    toc;

    save(['outputs/FEmat_', opts.suffix, '_', num2str(idS)], '-struct', 'FEmat', '-v7.3');
  end
  
  %% Compute Floquet-Bloch transform
  % Floquet points
  FloquetPoints = linspace(-pi/opts.period, pi/opts.period, opts.numFloquetPoints);

  % Homogeneous Dirichlet boundary conditions
  ecs = ((u|edge0x) == 0.0) & ((u|edge1x) == 0.0) &...
        ((u|edge0y) == 0.0) & ((u|edge1y) == 0.0);

  % Surfacic rhs
  B0x = sparse(edge0x.IdPoints(1:end), (1:N0x), 1.0, N_XY, N0x);
  B1x = sparse(edge1x.IdPoints(1:end), (1:N1x), 1.0, N_XY, N1x);
  B0y = sparse(edge0y.IdPoints(1:end), (1:N0y), 1.0, N_XY, N0y);
  B1y = sparse(edge1y.IdPoints(1:end), (1:N1y), 1.0, N_XY, N1y);

  ecs.applyEcs;
  ecs.b = [B0x, B1x, B0y, B1y];
  sidenames = {'0x', '1x', '0y', '1y'};

  % Define shear map
  % shearmapX = @(Y, s) spBX.evaluateBasisFunctions([opts.cutvec(1) * Y, opts.cutvec(2) * Y + s]);
  % shearmap0y = @(X, s) spBY.evaluateBasisFunctions([s, X]);
  % shearmap1y = @(X, s) spBY.evaluateBasisFunctions([cutslope + s, X]);

  % shearmapX(meshXY.points(edge0x.IdPoints, 1), sVar - cutslope);
  % shearmapX(meshXY.points(edge1x.IdPoints, 1), sVar + cutslope);
  % shearmap0y(meshXY.points(edge0y.IdPoints, 1), sVar - cutslope);
  % shearmap0y(meshXY.points(edge1y.IdPoints, 1), sVar + cutslope);
  
  % shearmap.fun0x = @(s) spBX.evaluateBasisFunctions([zeros(N0x, 1), opts.cutvec(1) * meshXY.points(edge0x.IdPoints, 1),...
  %                                                               mod(opts.cutvec(2) * meshXY.points(edge0x.IdPoints, 1) + s, 1)]);
  % %
  % shearmap.fun1x = @(s) spBX.evaluateBasisFunctions([ones(N0x,  1), opts.cutvec(1) * meshXY.points(edge1x.IdPoints, 1),...
  %                                                               mod(opts.cutvec(2) * meshXY.points(edge1x.IdPoints, 1) + s, 1)]);
  %
  shearmap.fun0x = @(s) spBX.evaluateBasisFunctions([    opts.cutvec(1) * meshXY.points(edge0x.IdPoints, 2),...
                                                     mod(opts.cutvec(2) * meshXY.points(edge0x.IdPoints, 2) + s, 1), zeros(N0x, 1)]);
  %
  shearmap.fun1x = @(s) spBX.evaluateBasisFunctions([    opts.cutvec(1) * meshXY.points(edge1x.IdPoints, 2),...
                                                     mod(opts.cutvec(2) * meshXY.points(edge1x.IdPoints, 2) + s, 1), zeros(N0x, 1)]);
  %
  shearmap.fun0y = @(s) spBY.evaluateBasisFunctions([mod(           s * ones(N0y, 1), 1), meshXY.points(edge0y.IdPoints, 1), zeros(N0y, 1)]);
  shearmap.fun1y = @(s) spBY.evaluateBasisFunctions([mod(cutslope + s * ones(N0y, 1), 1), meshXY.points(edge1y.IdPoints, 1), zeros(N0y, 1)]);
  
  
  for idFB = 1:opts.numFloquetPoints
    %% Local cell problems
    fprintf('<strong>Problèmes de cellule locaux.</strong>\n');
    fprintf('%d sur %d\n', idFB, opts.numFloquetPoints);
    
    FloquetVar = FloquetPoints(idFB);
    
    for idS = 1:N_s-1
      % The FE matrix is a linear combination of elementary matrices
      FEmat = load(['FEmat_', opts.suffix, '_', num2str(idS), '.mat']);

      AA = 1i * FloquetVar * FEmat.mat_vecYu_gradv + FEmat.mat_gradu_gradv...
         - 1i * FloquetVar * FEmat.mat_gradu_vecYv + FloquetVar * FloquetVar * FEmat.mat_vecYu_vecYv...
         - (omega^2) * FEmat.mat_u_v;
      
      AA0 =  ecs.P * AA * ecs.P';
      LL0 = -ecs.P * AA * ecs.b;
      
      % Solve the linear system
      tic;
      Ecell0 = ecs.b + ecs.P' * (AA0 \ LL0);
      tpscpu = toc;
      fprintf('Inversion système : %0.3e secondes\n', tpscpu);

      % Local cell solutions
      solcell.E0x = Ecell0(:, 1:N0x); Ecell0(:, 1:N0x) = [];
      solcell.E1x = Ecell0(:, 1:N1x); Ecell0(:, 1:N1x) = [];
      solcell.E0y = Ecell0(:, 1:N0y); Ecell0(:, 1:N0y) = [];
      solcell.E1y = Ecell0(:, 1:N1y);

      % local edge DtN operators
      tic;
      for idI = 1:4
        for idJ = 1:4
          nameEi  = ['E', sidenames{idI}];
          nameEj  = ['E', sidenames{idJ}];
          nameTij = ['edgeT', sidenames{idI}, sidenames{idJ}];
          solcell.(nameTij) = solcell.(nameEj)' * (AA * solcell.(nameEi));
        end
      end
      tpscpu = toc;
      fprintf('Calcul DtN : %0.3e secondes\n', tpscpu);

      % Save local cell data
      save(['outputs/local_cell_sol_', opts.suffix, '_s_', num2str(idS), '_Floquet_', num2str(idFB)], '-struct', 'solcell', '-v7.3');
    end
    system(['cp outputs/local_cell_sol_', opts.suffix, '_s_1_Floquet_', num2str(idFB), '.mat outputs/local_cell_sol_', opts.suffix, '_s_', num2str(N_s), '_Floquet_', num2str(idFB), '.mat']);

    %% Auxiliary face DtN operators
    % Initialize the operators
    aux.T0x0x = zeros(NbX, NbX); aux.T0x1x = aux.T0x0x; aux.T1x0x = aux.T0x0x; aux.T1x1x = aux.T0x0x; 
    aux.T0x0y = zeros(NbY, NbX); aux.T0x1y = aux.T0x0y; aux.T1x0y = aux.T0x0y; aux.T1x1y = aux.T0x0y; 
    aux.T0y0x = zeros(NbX, NbY); aux.T0y1x = aux.T0y0x; aux.T1y0x = aux.T0y0x; aux.T1y1x = aux.T0y0x; 
    aux.T0y0y = zeros(NbY, NbY); aux.T0y1y = aux.T0y0y; aux.T1y0y = aux.T0y0y; aux.T1y1y = aux.T0y0y;

    % Quadrature rule
    for idS = 1:N_s-1
      sVar = svec(idS);
      solcell = load(['outputs/local_cell_sol_', opts.suffix, '_s_', num2str(idS), '_Floquet_', num2str(idFB)]);

      % solcell.edgeT0x1x_m_slope = interpolateOperator([0, mod(sVar - cutslope, 1), 0], 'local_cell_sol', 'edgeT0x1x', sLine);
      % solcell.edgeT1x1x_m_slope = interpolateOperator([0, mod(sVar - cutslope, 1), 0], 'local_cell_sol', 'edgeT1x1x', sLine);

      solcell.edgeT0x1y_m_slope = interpolateOperator([0, mod(sVar - cutslope, 1), 0], 'local_cell_sol', 'edgeT0x1y', sLine);
      solcell.edgeT1x1y_m_slope = interpolateOperator([0, mod(sVar - cutslope, 1), 0], 'local_cell_sol', 'edgeT1x1y', sLine);

      % solcell.edgeT0y1x_m_slope = interpolateOperator([0, mod(sVar - cutslope, 1), 0], 'local_cell_sol', 'edgeT0y1x', sLine);
      % solcell.edgeT1y1x_m_slope = interpolateOperator([0, mod(sVar - cutslope, 1), 0], 'local_cell_sol', 'edgeT1y1x', sLine);

      solcell.edgeT0y1y_m_slope = interpolateOperator([0, mod(sVar - cutslope, 1), 0], 'local_cell_sol', 'edgeT0y1y', sLine);
      solcell.edgeT1y1y_m_slope = interpolateOperator([0, mod(sVar - cutslope, 1), 0], 'local_cell_sol', 'edgeT1y1y', sLine);

      % aux.T0x0x = aux.T0x0x + shearmap.fun0x(sVar)' * solcell.edgeT0x0x         * shearmap.fun0x(sVar);
      % aux.T0x1x = aux.T0x1x + shearmap.fun1x(sVar)' * solcell.edgeT0x1x_m_slope * shearmap.fun0x(sVar - cutslope);
      % aux.T1x0x = aux.T1x0x + shearmap.fun0x(sVar)' * solcell.edgeT1x0x         * shearmap.fun1x(sVar + cutslope);
      % aux.T1x1x = aux.T1x1x + shearmap.fun1x(sVar)' * solcell.edgeT1x1x_m_slope * shearmap.fun1x(sVar);

      % aux.T0x0y = aux.T0x0y + shearmap.fun0y(sVar)' * solcell.edgeT0x0y         * shearmap.fun0x(sVar);
      % aux.T0x1y = aux.T0x1y + shearmap.fun1y(sVar)' * solcell.edgeT0x1y_m_slope * shearmap.fun0x(sVar - cutslope);
      % aux.T1x0y = aux.T1x0y + shearmap.fun0y(sVar)' * solcell.edgeT1x0y         * shearmap.fun1x(sVar + cutslope);
      % aux.T1x1y = aux.T1x1y + shearmap.fun1y(sVar)' * solcell.edgeT1x1y_m_slope * shearmap.fun1x(sVar);

      % aux.T0y0x = aux.T0y0x + shearmap.fun0x(sVar)' * solcell.edgeT0y0x         * shearmap.fun0y(sVar);
      % aux.T0y1x = aux.T0y1x + shearmap.fun1x(sVar)' * solcell.edgeT0y1x_m_slope * shearmap.fun0y(sVar - cutslope);
      % aux.T1y0x = aux.T1y0x + shearmap.fun0x(sVar)' * solcell.edgeT1y0x         * shearmap.fun1y(sVar + cutslope);
      % aux.T1y1x = aux.T1y1x + shearmap.fun1x(sVar)' * solcell.edgeT1y1x_m_slope * shearmap.fun1y(sVar);

      % aux.T0y0y = aux.T0y0y + shearmap.fun0y(sVar)' * solcell.edgeT0y0y         * shearmap.fun0y(sVar);
      % aux.T0y1y = aux.T0y1y + shearmap.fun1y(sVar)' * solcell.edgeT0y1y_m_slope * shearmap.fun0y(sVar - cutslope);
      % aux.T1y0y = aux.T1y0y + shearmap.fun0y(sVar)' * solcell.edgeT1y0y         * shearmap.fun1y(sVar + cutslope);
      % aux.T1y1y = aux.T1y1y + shearmap.fun1y(sVar)' * solcell.edgeT1y1y_m_slope * shearmap.fun1y(sVar);
      aux.T0x0x = aux.T0x0x + shearmap.fun0x(sVar)' * solcell.edgeT0x0x * shearmap.fun0x(sVar);
      aux.T0x1x = aux.T0x1x + shearmap.fun1x(sVar)' * solcell.edgeT0x1x * shearmap.fun0x(sVar);
      aux.T1x0x = aux.T1x0x + shearmap.fun0x(sVar)' * solcell.edgeT1x0x * shearmap.fun1x(sVar);
      aux.T1x1x = aux.T1x1x + shearmap.fun1x(sVar)' * solcell.edgeT1x1x * shearmap.fun1x(sVar);

      aux.T0x0y = aux.T0x0y + shearmap.fun0y(sVar)' * solcell.edgeT0x0y         * shearmap.fun0x(sVar);
      aux.T0x1y = aux.T0x1y + shearmap.fun1y(sVar)' * solcell.edgeT0x1y_m_slope * shearmap.fun0x(sVar - cutslope);
      aux.T1x0y = aux.T1x0y + shearmap.fun0y(sVar)' * solcell.edgeT1x0y         * shearmap.fun1x(sVar);
      aux.T1x1y = aux.T1x1y + shearmap.fun1y(sVar)' * solcell.edgeT1x1y_m_slope * shearmap.fun1x(sVar - cutslope);

      aux.T0y0x = aux.T0y0x + shearmap.fun0x(sVar)' * solcell.edgeT0y0x * shearmap.fun0y(sVar);
      aux.T0y1x = aux.T0y1x + shearmap.fun1x(sVar)' * solcell.edgeT0y1x * shearmap.fun0y(sVar);
      aux.T1y0x = aux.T1y0x + shearmap.fun0x(sVar)' * solcell.edgeT1y0x * shearmap.fun1y(sVar + cutslope);
      aux.T1y1x = aux.T1y1x + shearmap.fun1x(sVar)' * solcell.edgeT1y1x * shearmap.fun1y(sVar + cutslope);

      aux.T0y0y = aux.T0y0y + shearmap.fun0y(sVar)' * solcell.edgeT0y0y         * shearmap.fun0y(sVar);
      aux.T0y1y = aux.T0y1y + shearmap.fun1y(sVar)' * solcell.edgeT0y1y_m_slope * shearmap.fun0y(sVar - cutslope);
      aux.T1y0y = aux.T1y0y + shearmap.fun0y(sVar)' * solcell.edgeT1y0y         * shearmap.fun1y(sVar + cutslope);
      aux.T1y1y = aux.T1y1y + shearmap.fun1y(sVar)' * solcell.edgeT1y1y_m_slope * shearmap.fun1y(sVar);

      % % Shift some of the edge operators
      % for idI = 1:4
      %   for idJ = [2,4] % '1x' and '1y'
      %     nameTij = ['edgeT', sidenames{idI}, sidenames{idJ}];
      %     solcell.(nameTij) = interpolateOperator([0, mod(sVar - cutslope, 1), 0], 'local_cell_sol', nameTij, sLine);
      %   end
      % end

      % % Compute the auxiliary face operators
      % for idI = 1:4
      %   namefunI = ['fun', sidenames{idI}];

      %   for idJ = 1:4
      %     nameEdgeT = ['edgeT', sidenames{idI}, sidenames{idJ}];
      %     nameAuxT  = ['T',  sidenames{idI}, sidenames{idJ}];
      %     namefunJ = ['fun', sidenames{idJ}];
      %     shiftpas = (str2double(sidenames{idI}(1)) - str2double(sidenames{idJ}(1))) * cutslope;

      %     aux.(nameAuxT) = aux.(nameAuxT) + shearmap.(namefunJ)(sVar)' * solcell.(nameEdgeT) * shearmap.(namefunI)(sVar + shiftpas);
      %   end
      % end
    end

    % Divide by the quadrature weight
    for idI = 1:4
      for idJ = 1:4
        nameAuxT  = ['T',  sidenames{idI}, sidenames{idJ}];
        aux.(nameAuxT) = aux.(nameAuxT) * DeltaS / opts.cutvec(1);
      end
    end

    %% Compute DtD operator
    fprintf('<strong>Calcul opérateur DtD.</strong>\n');
    spy(aux.T1y0y + aux.T0y0y + aux.T0y1y + aux.T1y1y);
    DtD = (aux.T1y0y + aux.T0y0y + aux.T0y1y + aux.T1y1y) \ [aux.T0x0y + aux.T0x1y, aux.T1x0y + aux.T1x1y];
    DtD0 = DtD(:,     1:NbX);
    DtD1 = DtD(:, NbX+1:end);

    %% Traces and normal traces of the local cell solutions
    % Ekl is the trace of Ek on Sigmal
    E00 = speye(NbX);
    E10 = sparse(NbX, NbX);
    E01 = sparse(NbX, NbX);
    E11 = speye(NbX);

    % Fkl is the normal trace of Ek on Sigmal. It is evaluated weakly and then projected
    fprintf('<strong>Opérateurs DtN locaux.</strong>\n');
    F00 = spBX.invmassmat * ((aux.T0y0x + aux.T1y0x) * DtD0 + aux.T0x0x);
    F10 = spBX.invmassmat * ((aux.T0y0x + aux.T1y0x) * DtD1 + aux.T1x0x);
    F01 = spBX.invmassmat * ((aux.T0y1x + aux.T1y1x) * DtD0 + aux.T0x1x);
    F11 = spBX.invmassmat * ((aux.T0y1x + aux.T1y1x) * DtD1 + aux.T1x1x);
    
    %% Solve the Riccati equation
    fprintf('<strong>Equation de Riccati.</strong>\n');
    
    flux = @(V) modesFlux(V, E00, E10, F00, F10, orientation, spB0, opts.omega);
    riccatiOpts.tol = 1.0e-2;
    if isfield(opts, 'suffix')
      riccatiOpts.suffix = opts.suffix;
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
    save(['outputs/half_guide_sol_', opts.suffix, '_Floquet_', num2str(idFB)], '-struct', 'solguide', '-v7.3');
  end

end

%% Auxiliary function: The flux function
function flux = modesFlux(V, E00, E10, F00, F10, orientation, spB, omega)
  
  Nb = size(E00, 1);
  Psi  = E00 * V(1:Nb, :) + E10 * V(Nb+1:end, :);
  dPsi = F00 * V(1:Nb, :) + F10 * V(Nb+1:end, :);

  flux = -orientation * imag(diag(Psi' * spB.massmat * dPsi) / omega);

end


%% Auxiliary function: interpolateOperator()
function opres = interpolateOperator(P, prefix, op_name, dom)
  structLoc = dom.locateInDomain(P);
  elts = dom.elements(structLoc.elements, :);
  coos = structLoc.barycoos;
  opres = cell(size(P, 1), 1);

  for idI = 1:size(P, 1)
    % Initialization
    solcell = load([prefix, '_', num2str(find(dom.IdPoints == elts(idI, 1))), '.mat']);
    op = solcell.(op_name);
    opres{idI} = coos(idI, 1) * op;

    % Update
    for idJ = 2:size(coos, 2)
      solcell = load([prefix, '_', num2str(find(dom.IdPoints == elts(idI, idJ))), '.mat']);
      op = solcell.(op_name);
      opres{idI} = opres{idI} + coos(idI, idJ) * op;
    end
  end

  if (size(P, 1) == 1)
    opres = opres{1};
  end
end