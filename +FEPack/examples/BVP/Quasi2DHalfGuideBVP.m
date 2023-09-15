 function [U, E0, E1, R, D, newBCstruct, Lambda] = Quasi2DHalfGuideBVP(meshDir, meshTrs, meshAux, volBilinearIntg, BCstruct, numCells, opts)
  
  % % ************************* %
   % Preliminary verifications %
   % ************************* %
  if ~isfield(opts, 'computeSol')
    opts.computeSol = true;
  end

  if ~isfield(opts, 'solBasis')
    opts.solBasis = false;
  end

  if ~isfield(opts, 'verbose')
    opts.verbose = 1;
  end

  %% % ************** %
  %  % Initialization %
  %  % ************** %
  Ndir = meshDir.numPoints;
  Spoints = meshTrs.points(meshTrs.domain('xmin').IdPoints, 2);
  Ns = numel(Spoints);
  DeltaS = Spoints(2) - Spoints(1);

  if isa(volBilinearIntg, 'FEPack.pdes.Form')
    % Compute the FE matrix if not done already
    AA = FEPack.pdes.Form.intg(meshDir.domain('volumic'), volBilinearIntg);
  else
    AA = volBilinearIntg;
  end

  % % Domains
  % Xeq1 = meshDir.domains{1};
  % Xeq0 = meshDir.domains{2};
  % Yeq1 = meshDir.domains{3};
  % Yeq0 = meshDir.domains{4};

  % % Boundary basis functions
  % spBX0 = BCstruct.spBX0;
  % spBX1 = BCstruct.spBX1;
  % NbX = spBX0.numBasis; % Number of basis functions

  % spBY0 = BCstruct.spBY0;
  % spBY1 = BCstruct.spBY1;
  % NbY = spBY0.numBasis; % Number of basis functions

  %% % ***************************** %
  %  % Solve the local cell problems %
  %  % ***************************** %
  if (opts.verbose)
    fprintf('1. Resolution des problemes de cellule locaux\n');
  end

  %%%%%%%%%%%%%%%%%%
  % Dirichlet only %
  %%%%%%%%%%%%%%%%%%
  ecs = FEPack.pdes.EssentialConditions;
  for idom = 1:4
    ecs = ecs & assignEcs( (u|meshDir.domains{idom}) == 0.0 );
  end
  ecs.applyEcs;

   % Add the surfacic contributions
  Ndir2D = 2*meshDir.domain('xmin') + 2*meshDir.domain('ymin') - 4;
  ecs.b = sparse(Ndir, Ndir2D);
  NbXdir = meshDir.domain('xmin').numPoints;
  NbYdir = meshDir.domain('ymin').numPoints;

  FX0 = cell(Ns, 1);
  FX1 = cell(Ns, 1);
  FY0 = cell(Ns, 1);
  FY1 = cell(Ns, 1);
  tX0Y0 = cell(Ns, 1);
  tX0Y1 = cell(Ns, 1);
  tX1Y0 = cell(Ns, 1);
  tX1Y1 = cell(Ns, 1);
  tY0Y0 = cell(Ns, 1);
  tY0Y1 = cell(Ns, 1);
  tY1Y0 = cell(Ns, 1);
  tY1Y1 = cell(Ns, 1);
  
  for idTrs = 1:Ns-1

    AA0 =  ecs.P * AA{idTrs} * ecs.P';
    LL0 = -ecs.P * AA{idTrs} * ecs.b;

    UU0 = AA0 \ LL0;
    Fsol = ecs.b + ecs.P' * UU0;
    
    FX0{idTrs} = Fsol(1:NbX);
    FX1{idTrs} = Fsol(NbX:2*NbX-1);
    FY0{idTrs} = Fsol(2*NbX-1:2*NbX+NbY-2);
    FY1{idTrs} = Fsol(2*NbX+NbY-2:end);

    tX0Y0{idTrs} = FY0' * AA{idTrs} * FX0;
    tX0Y1{idTrs} = FY1' * AA{idTrs} * FX0;
    tX1Y0{idTrs} = FY0' * AA{idTrs} * FX1;
    tX1Y1{idTrs} = FY1' * AA{idTrs} * FX1;

    tY0Y0{idTrs} = FY0' * AA{idTrs} * FY0;
    tY0Y1{idTrs} = FY1' * AA{idTrs} * FY0;
    tY1Y0{idTrs} = FY0' * AA{idTrs} * FY1;
    tY1Y1{idTrs} = FY1' * AA{idTrs} * FY1;

  end

  FX0{Ns} = FX0{1};
  FX1{Ns} = FX1{1};
  FY0{Ns} = FY0{1};
  FY1{Ns} = FY1{1};
  tX0Y0{Ns} = tX0Y0{1};
  tX0Y1{Ns} = tX0Y1{1};
  tX1Y0{Ns} = tX1Y0{1};
  tX1Y1{Ns} = tX1Y1{1};
  tY0Y0{Ns} = tY0Y0{1};
  tY0Y1{Ns} = tY0Y1{1};
  tY1Y0{Ns} = tY1Y0{1};
  tY1Y1{Ns} = tY1Y1{1};

  % Compute t01_{s-delta} and t11_{s-delta}
  delta = opts.cutvec(2) / opts.cutvec(1);
  deltavec = zeros(Ndir, 1) * [0, 0, delta];

  tX0Y1_m_delta = cell(Ns, 1);
  tX1Y1_m_delta = cell(Ns, 1);
  tY0Y1_m_delta = cell(Ns, 1);
  tY1Y1_m_delta = cell(Ns, 1);
  
  structLoc = meshTrs.domain('xmin').locateInDomain(Spoints - delta);
  elts = meshTrs.domain('xmin').elements(structLoc.elements, :)
  coos = structLoc.barycoos;

  for idS = 1:Ns 

    tX0Y1_m_delta{idS} = coos(idS, 1) * tX0Y1{elts(idS, 1)} + coos(idS, 2) * tX0Y1{elts(idS, 2)};
    tX1Y1_m_delta{idS} = coos(idS, 1) * tX1Y1{elts(idS, 1)} + coos(idS, 2) * tX1Y1{elts(idS, 2)};

    tY0Y1_m_delta{idS} = coos(idS, 1) * tY0Y1{elts(idS, 1)} + coos(idS, 2) * tY0Y1{elts(idS, 2)};
    tY1Y1_m_delta{idS} = coos(idS, 1) * tY1Y1{elts(idS, 1)} + coos(idS, 2) * tY1Y1{elts(idS, 2)};

  end


  % %%%%%%%%%%%%%%%%%%%% %
  % Compute DtD operator %
  % %%%%%%%%%%%%%%%%%%%% %
  Ntrs = BCstruct.spBtrs.numBasis;
  Naux = BCstruct.spBaux.numBasis;

  TX0Y0 = zeros(Ntrs, Naux);
  TX0Y1 = zeros(Ntrs, Naux);
  TX1Y0 = zeros(Ntrs, Naux);
  TX1Y1 = zeros(Ntrs, Naux);

  TY0Y0 = zeros(Naux, Naux);
  TY0Y1 = zeros(Naux, Naux);
  TY1Y0 = zeros(Naux, Naux);
  TY1Y1 = zeros(Naux, Naux);
  

  for idS = 1:Ns
    
    x3Dtrs = [zeros(Ndir, 1), opts.cutvec(1)*meshDir.points(:, 2),...
                              opts.cutvec(2)*meshDir.points(:, 2) + Spoints(idS)];
    phisTrs = BCstruct.spBtrs.phis(x3Dtrs, 1:Ntrs);
    phisTrs_p_delta = BCstruct.spBtrs.phis(x3Dtrs + deltavec, 1:Ntrs);
    phisTrs_m_delta = BCstruct.spBtrs.phis(x3Dtrs - deltavec, 1:Ntrs);

    x3D_Y0 = [meshDir.points(:, 1), zeros(Ndir, 1), Spoints(idS) * ones(Ndir, 1)];
    x3D_Y1 = [meshDir.points(:, 1), ones(Ndir, 1), deltavec + Spoints(idS)];
    phisAux_Y0 = BCstruct.spBaux.phis(x3D_Y0, 1:Naux);
    phisAux_Y1 = BCstruct.spBaux.phis(x3D_Y1, 1:Naux);
    
    phisAux_Y0_m_delta = BCstruct.spBaux.phis(x3D_Y0 - deltavec, 1:Naux);
    phisAux_Y1_p_delta = BCstruct.spBaux.phis(x3D_Y1 + deltavec, 1:Naux);

    TX0Y0 = TX0Y0 + DeltaS * phisAux_Y0' * tX0Y0{idS} * phisTrs;
    TX0Y1 = TX0Y1 + DeltaS * phisAux_Y1' * tX0Y1_m_delta{idS} * phisTrs_m_delta;
    TX1Y0 = TX1Y0 + DeltaS * phisAux_Y0' * tX1Y0{idS} * phisTrs_p_delta;
    TX1Y1 = TX1Y1 + DeltaS * phisAux_Y1' * tX1Y1_m_delta{idS} * phisTrs;
    
    TY0Y0 = TY0Y0 + DeltaS * phisAux_Y0' * tY0Y0{idS} * phisAux_Y0;
    TY0Y1 = TY0Y1 + DeltaS * phisAux_Y1' * tY0Y1_m_delta{idS} * phisAux_Y0_m_delta;
    TY1Y0 = TY1Y0 + DeltaS * phisAux_Y0' * tY1Y0{idS} * phisAux_Y1_p_delta;
    TY1Y1 = TY1Y1 + DeltaS * phisAux_Y1' * tY1Y1_m_delta{idS} * phisAux_Y1;

  end

  DtD = (TY0Y0 + TY0Y1 + TY1Y0 + TY1Y1) \ [TX0Y0 + TX0Y1, TX1Y0 + TX1Y1]; % Naux x Ntrs

  % ************************ %
  % Deduce the DtN operators %
  % ************************ %
  % T00 = 


  %% % ************************** %
  %  % Solve the Riccati equation %
  %  % ************************** %
  if (opts.verbose)
    fprintf('3. Resolution du système de Riccati\n');
  end

  % Solve the linearized eigenvalue problem associated to the Riccati equation
  flux = @(V) modesFlux(V, E00, E10, F00, F10, orientation, spB0, opts.omega);
  riccatiOpts.tol = 1.0e-2;
  if isfield(opts, 'suffix')
    riccatiOpts.suffix = opts.suffix;
  else
    riccatiOpts.suffix = '';
  end
  [R, D] = propagationOperators([E01, E11;  orientation*F01,  orientation*F11], ...
                                [E00, E10; -orientation*F00, -orientation*F10], flux, riccatiOpts);

  %% % ********************************* %
  %  % Compute the solution cell by cell %
  %  % ********************************* %
  if (opts.verbose)
    fprintf('4. Construction de la solution\n');
  end

  if (opts.computeSol)
    if (opts.solBasis)

      % Compute the solution for any basis function
      U = cell(numCells, 1);
      R0 = eye(Nb);
      R1 = D;

      for idCell = 0:numCells-1
        % Compute the solution in the current cell
        U{idCell + 1} = E0 * R0 + E1 * R1;

        % Update the coefficients
        R0 = R * R0;
        R1 = D * R0;
      end

    else

      % Compute the solution for one boundary data
      U = zeros(N, numCells);
      R0Phi = spB0.FE_to_spectral * BCstruct.phi(mesh.points(Sigma0.IdPoints, :));
      R1Phi = D * R0Phi;

      for idCell = 0:numCells-1
        % Compute the solution in the current cell
        U(:, idCell + 1) = E0 * R0Phi + E1 * R1Phi;

        % Update the coefficients
        R0Phi = R * R0Phi;
        R1Phi = D * R0Phi;
      end

    end
  else
    U = [];
  end
  %% % *********************************************** %
  %  % The transmission coefficient and the derivative %
  %  % *********************************************** %
  if (opts.verbose)
    fprintf('5. Calcul de l''operateur de transmission global\n');
  end

  % Compute the trace and the normal trace of the half-guide solution
  U0 = E00 + E10 * D;
  dU0 = F00 + F10 * D;

  Lambda = -BCstruct.BCu' * dU0 + BCstruct.BCdu' * U0;
  newBCstruct = BCstruct;
end

%% % ******** %
%  % Appendix %
%  % ******** %
% The flux function
function flux = modesFlux(V, E00, E10, F00, F10, orientation, spB, omega)
  
  Nb = size(E00, 1);
  Psi  = E00 * V(1:Nb, :) + E10 * V(Nb+1:end, :);
  dPsi = F00 * V(1:Nb, :) + F10 * V(Nb+1:end, :);

  flux = -orientation * imag(diag(Psi' * spB.massmat * dPsi) / omega);

end
