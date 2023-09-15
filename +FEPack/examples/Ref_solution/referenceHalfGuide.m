function referenceHalfGuide(meshcell, orientation, omega, mu_coeff, rho_coeff, numCellsZ, sizeCellZ, Zorigin, basis_function, suffix)

  fprintf('1. Problèmes de cellule locaux\n');
  Sigma0x = meshcell.domain('xmin'); N0x = Sigma0x.numPoints;
  Sigma1x = meshcell.domain('xmax'); N1x = Sigma1x.numPoints;
  Sigma0z = meshcell.domain('ymin'); N0z = Sigma0z.numPoints;
  Sigma1z = meshcell.domain('ymax'); N1z = Sigma1z.numPoints;
  domCell = meshcell.domain('volumic');
  N = meshcell.numPoints;
  u = FEPack.pdes.PDEObject; v = dual(u);

  % Homogeneous Dirichlet boundary conditions
  ecs = ((u|Sigma0x) == 0.0) & ((u|Sigma1x) == 0.0) &...
        ((u|Sigma0z) == 0.0) & ((u|Sigma1z) == 0.0);

  % Surfacic rhs
  B0x = sparse(Sigma0x.IdPoints, (1:N0x), 1.0, N, N0x);
  B1x = sparse(Sigma1x.IdPoints, (1:N1x), 1.0, N, N1x);
  B0z = sparse(Sigma0z.IdPoints(2:end-1), (1:N0z-2), 1.0, N, N0z-2);
  B1z = sparse(Sigma1z.IdPoints(2:end-1), (1:N1z-2), 1.0, N, N1z-2);

  ecs.applyEcs;
  ecs.b = [B0x, B1x, B0z, B1z];
  sidenames = {'0x', '1x', '0z', '1z'};

  % ------------------------- %
  % local local cell problems %
  % ------------------------- %
  for idZ = 1:numCellsZ
    fprintf('%3d sur %3d\n', idZ, numCellsZ);
    solcell = struct;
    
    % Compute FE matrices
    cellZorigin = sizeCellZ * (idZ - 1) + Zorigin;

    mu_cell  = @(x)  mu_coeff([x(:, 1), x(:, 2) + cellZorigin]);
    rho_cell = @(x) rho_coeff([x(:, 1), x(:, 2) + cellZorigin]);

    MM = FEPack.pdes.Form.intg(domCell, rho_cell * u * v);
    KK = FEPack.pdes.Form.intg(domCell, (mu_cell * grad2(u)) * grad2(v));
    AA = KK - (omega^2) * MM;

    % Elimination
    AA0 =  ecs.P * AA * ecs.P';
    LL0 = -ecs.P * AA * ecs.b;

    % Solve the linear system
    Ecell0 = ecs.b + ecs.P' * (AA0 \ LL0);
    
    % Local cell solutions
    solcell.E0x = Ecell0(:, 1:N0x  ); Ecell0(:, 1:N0x  ) = [];
    solcell.E1x = Ecell0(:, 1:N1x  ); Ecell0(:, 1:N1x  ) = [];
    solcell.E0z = Ecell0(:, 1:N0z-2); Ecell0(:, 1:N0z-2) = [];
    solcell.E1z = Ecell0(:, 1:N1z-2);

    % DtN operators
    for idI = 1:4
      for idJ = 1:4
        nameEi  = ['E', sidenames{idI}];
        nameEj  = ['E', sidenames{idJ}];
        nameTij = ['T', sidenames{idI}, sidenames{idJ}];
        solcell.(nameTij) = solcell.(nameEj)' * (AA * solcell.(nameEi));
      end
    end

    % Save output
    save(['outputs/local_cell_sol_' suffix, '_', int2str(idZ)], '-struct', 'solcell');
  end

  if (false)
    figure;
    set(groot,'defaultAxesTickLabelInterpreter','latex');
    set(groot,'defaulttextinterpreter','latex');
    set(groot,'defaultLegendInterpreter','latex');

    for idZ = 1:numCellsZ
      SC = load(['outputs/local_cell_sol_' suffix, '_', int2str(idZ)]);

      for idI = 1:size(SC.E0z, 2)
        trisurf(meshcell.triangles, meshcell.points(:, 1),...
                                    meshcell.points(:, 2), full(real(SC.E0z(:, idI))));
        view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
        set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);

        pause;
      end
    end
  end

  % ------------- %
  % Jacobi system %
  % ------------- %
  fprintf('2. Système de Jacobi\n');
  Njaco = (N0z - 2) * (numCellsZ - 1);
  AAjaco = sparse(Njaco, Njaco);

  % Jacobi matrix
  fprintf('\tMatrice de Jacobi\n');

  for idZ = 1:numCellsZ-1 % Main diagonal of Jacobi matrix
    SC1 = load(['outputs/local_cell_sol_' suffix, '_', int2str(idZ)]);
    SC2 = load(['outputs/local_cell_sol_' suffix, '_', int2str(idZ+1)]);
    
    II = (1+(idZ-1)*(N0z-2):idZ*(N0z-2))' * ones(1, N0z-2);
    JJ = II';
    AAjaco = AAjaco + sparse(II(:), JJ(:), SC2.T0z0z(:) + SC1.T1z1z(:), Njaco, Njaco);
  end

  for idZ = 1:numCellsZ-2 % Lower and upper diagonals of Jacobi matrix
    SC = load(['outputs/local_cell_sol_' suffix, '_', int2str(idZ+1)]);
    
    II      = (1+(idZ-1)*(N0z-2):idZ*(N0z-2))' * ones(1, N0z-2);
    IIplus1 = (1+idZ*(N0z-2):(idZ+1)*(N0z-2))' * ones(1, N0z-2);
    JJ      = II';
    JJplus1 = IIplus1';
    AAjaco  = AAjaco + sparse(IIplus1(:), JJ(:), SC.T0z1z(:), Njaco, Njaco)...
                     + sparse(II(:), JJplus1(:), SC.T1z0z(:), Njaco, Njaco);
  end

  %% Right-hand sides
  fprintf('\tSeconds membres\n');
  numPointsXcst = numCellsZ * (N0x - 1) + 1;  % Number of nodes on x = cst
  
  if strcmp(basis_function.name, 'Lagrange')

    B0 = speye(numPointsXcst, numPointsXcst);
    % B0(:, 1) = [];
    % B0(:, end) = [];
    B1 = B0;

  elseif strcmp(basis_function.name, 'Fourier')

    Sigma0xPoints = sort(meshcell.points(Sigma0x.IdPoints, 2));
    pointsXcst = Sigma0xPoints * ones(1, numCellsZ)...
               + sizeCellZ * ones(Sigma0x.numPoints, 1) * (0:numCellsZ-1) + Zorigin;
    
    pointsXcst(1, :) = [];
    pointsXcst = [Sigma0xPoints(1) + Zorigin; pointsXcst(:)];

    B0 = exp(2i * pi * pointsXcst * basis_function.FourierIds(:)' / (sizeCellZ * numCellsZ));
    B1 = B0;
    
  else

    error(['Paramètre ', basis_function.name, ' non reconnu.']);

  end

  numBasisX = size(B0, 2);
  LL0 = sparse(Njaco, numBasisX);% LL0b = LL0;
  LL1 = sparse(Njaco, numBasisX);% LL1b = LL1;

  for idZ = 1:numCellsZ-1
    SC1 = load(['outputs/local_cell_sol_' suffix, '_', int2str(idZ)]);
    SC2 = load(['outputs/local_cell_sol_' suffix, '_', int2str(idZ+1)]);

    idBasis = 1+(idZ-1)*(N0x-1):1+idZ*(N0x-1);

    II = (1+(idZ-1)*(N0z-2):idZ*(N0z-2))' * ones(1, size(B0, 2));
    JJ = ones(N0z-2, 1) * (1:size(B0, 2));

    VV0 = -SC1.T0x1z * B0(idBasis, :) - SC2.T0x0z * B0(idBasis+N0x-1, :);
    VV1 = -SC1.T1x1z * B1(idBasis, :) - SC2.T1x0z * B1(idBasis+N0x-1, :);

    LL0 = LL0 + sparse(II(:), JJ(:), VV0(:), Njaco, size(B0, 2));
    LL1 = LL1 + sparse(II(:), JJ(:), VV1(:), Njaco, size(B1, 2));
  end

  %% Solve system
  fprintf('\tRésolution du système\n');
  traceE  = AAjaco \ [LL0, LL1];

  traceE0 = cell(numCellsZ-1, 1);
  traceE1 = cell(numCellsZ-1, 1);
  for idZ = 1:numCellsZ-1
    idJaco = 1+(idZ-1)*(N0z-2):idZ*(N0z-2);
    traceE0{idZ} = traceE(idJaco, 1:size(LL0, 2));
    traceE1{idZ} = traceE(idJaco, 1+size(LL0, 2):end);
  end
  traceE0 = [{sparse(N0z-2, size(LL0, 2))}; traceE0; {sparse(N0z-2, size(LL0, 2))}];
  traceE1 = [{sparse(N1z-2, size(LL1, 2))}; traceE1; {sparse(N1z-2, size(LL1, 2))}];

  % Plot local cell solutions
  if (false)
    figure;
    set(groot,'defaultAxesTickLabelInterpreter','latex');
    set(groot,'defaulttextinterpreter','latex');
    set(groot,'defaultLegendInterpreter','latex');

    E0 = cell(numCellsZ, 1);
    E1 = cell(numCellsZ, 1);

    for idZ = 1:numCellsZ
      SC = load(['outputs/local_cell_sol_' suffix, '_', int2str(idZ)]);
      idBasis = 1+(idZ-1)*(N0x-1):1+idZ*(N0x-1);

      E0{idZ} = SC.E0x * B0(idBasis, :) + SC.E0z * traceE0{idZ} + SC.E1z * traceE0{idZ+1};
      E1{idZ} = SC.E1x * B1(idBasis, :) + SC.E0z * traceE1{idZ} + SC.E1z * traceE1{idZ+1};
    end

    for idI = 1:size(B1, 2) % (length(basis_function.FourierIds)+1)/2
      for idZ = 1:numCellsZ
        cellZorigin = sizeCellZ * (idZ - 1) + Zorigin;
        trisurf(meshcell.triangles, meshcell.points(:, 1),...
                                    meshcell.points(:, 2)+cellZorigin, full(real(E1{idZ}(:, idI))));
        hold on;
        view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
        set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);
        % caxis([-0.2, 1]);
      end
      hold off;
      pause;
    end
  end

  % --------------------------------------------------- %
  % Local DtN operators - Riccati - global DtN operator %
  % --------------------------------------------------- %
  % Compute mass matrix associated to boundary
  MMx = sparse(numBasisX, numBasisX);
  MM = FEPack.pdes.Form.intg(Sigma0x, u * v);
  MM = MM(Sigma0x.IdPoints, Sigma0x.IdPoints);

  for idZ = 1:numCellsZ
    idBasis = (1+(idZ-1)*(N0x-1):1+idZ*(N0x-1))';
    
    MMx = MMx + B0(idBasis, :)' * MM * B0(idBasis, :);
  end
  
  % Trace of local cell solutions
  E00 = speye(numBasisX, numBasisX); % MMx;
  E10 = sparse(numBasisX, numBasisX);
  E01 = sparse(numBasisX, numBasisX);
  E11 = speye(numBasisX, numBasisX); % MMx;

  % Normal trace of local cell solutions
  F00 = zeros(numBasisX, numBasisX);
  F01 = zeros(numBasisX, numBasisX);
  F10 = zeros(numBasisX, numBasisX);
  F11 = zeros(numBasisX, numBasisX);

  for idZ = 1:numCellsZ
    SC = load(['outputs/local_cell_sol_' suffix, '_', int2str(idZ)]);
    idBasis = 1+(idZ-1)*(N0x-1):1+idZ*(N0x-1);

    F00 = F00 + B0(idBasis, :)' * (SC.T0x0x * B0(idBasis, :) + SC.T0z0x * traceE0{idZ} + SC.T1z0x * traceE0{idZ+1});
    F01 = F01 + B0(idBasis, :)' * (SC.T0x1x * B0(idBasis, :) + SC.T0z1x * traceE0{idZ} + SC.T1z1x * traceE0{idZ+1});
    F10 = F10 + B1(idBasis, :)' * (SC.T1x0x * B1(idBasis, :) + SC.T0z0x * traceE1{idZ} + SC.T1z0x * traceE1{idZ+1});
    F11 = F11 + B1(idBasis, :)' * (SC.T1x1x * B1(idBasis, :) + SC.T0z1x * traceE1{idZ} + SC.T1z1x * traceE1{idZ+1});
  end

  % Riccati equation
  fprintf('3. Equation de Riccati\n');
  [Pop, Sop] = propagationOperators([E01, E11;  orientation*F01,  orientation*F11],...
                                    [E00, E10; -orientation*F00, -orientation*F10],...
                                    @(V, E00, E10, F00, F10) -orientation * imag(diag(...
                                          (E00 * V(1:Nb, :) + E10 * V(Nb+1:end, :))' *...
                                          (F00 * V(1:Nb, :) + F10 * V(Nb+1:end, :))) / opts.omega)); %#ok

  % Global DtN operator
  Lambda = F00 + F10 * Sop; %#ok

  % Save output
  save(['outputs/half_guide_solution_', suffix], 'traceE0', 'traceE1', 'Pop', 'Sop', 'Lambda', 'B0', 'B1');
end