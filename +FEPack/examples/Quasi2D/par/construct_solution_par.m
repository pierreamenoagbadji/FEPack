function construct_solution_par(orientation, meshXY, meshLineZ, BCstruct, opts)

  % Initialization
  sLine = meshLineZ.domain('volumic');
  cutslope = opts.cutvec(2) / opts.cutvec(1);
  spBX = BCstruct.spBX;
  spBY = BCstruct.spBY; 

  % Edges
  edge0x = meshXY.domain('xmin'); N0x = edge0x.numPoints;
  edge1x = meshXY.domain('xmax'); N1x = edge1x.numPoints;
  edge0y = meshXY.domain('ymin'); N0y = edge0y.numPoints;
  edge1y = meshXY.domain('ymax'); N1y = edge1y.numPoints;

  % Define shear map
  shearmap.fun0x = @(s) spBX.evaluateBasisFunctions([    opts.cutvec(1) * meshXY.points(edge0x.IdPoints, 2),...
                                                     mod(opts.cutvec(2) * meshXY.points(edge0x.IdPoints, 2) + s, 1), zeros(N0x, 1)]);
  %
  shearmap.fun1x = @(s) spBX.evaluateBasisFunctions([    opts.cutvec(1) * meshXY.points(edge1x.IdPoints, 2),...
                                                     mod(opts.cutvec(2) * meshXY.points(edge1x.IdPoints, 2) + s, 1), zeros(N1x, 1)]);
  %
  shearmap.fun0y = @(s) spBY.evaluateBasisFunctions([mod(           s * ones(N0y-2, 1), 1), meshXY.points(edge0y.IdPoints(2:end-1), 1), zeros(N0y-2, 1)]);
  shearmap.fun1y = @(s) spBY.evaluateBasisFunctions([mod(cutslope + s * ones(N1y-2, 1), 1), meshXY.points(edge1y.IdPoints(2:end-1), 1), zeros(N1y-2, 1)]);
  
  W = (2*pi/opts.period) ./ (opts.numFloquetPoints - 1);

  %% Traces of FB transforms
  suffix = opts.suffix;
  numCellsSemiInfinite = opts.numCellsSemiInfinite;
  numCellsInfinite = opts.numCellsInfinite;
  nomdossier = opts.nomdossier;

  parfor idFB = 1:opts.numFloquetPoints
    soltrace = load([nomdossier, 'sol_trace_Floquet_', num2str(idFB)]);
    solguide = load([nomdossier, 'half_guide_sol_', suffix, '_Floquet_', num2str(idFB)]);
    
    soltrace.R0phi = zeros(size(soltrace.vec, 1), numCellsSemiInfinite);
    soltrace.R0phi(:, 1) = soltrace.vec;

    for idX = 2:numCellsSemiInfinite
      soltrace.R0phi(:, idX) = solguide.R * soltrace.R0phi(:, idX-1);
    end
    soltrace.R1phi = solguide.D * soltrace.R0phi;

    % save([nomdossier, 'sol_trace_Floquet_', num2str(idFB)], 'R0phi', 'R1phi', '-append');
    parsave([nomdossier, 'sol_trace_Floquet_', num2str(idFB)], soltrace, true);
  end

  %% Compute inverse FB transform
  idsY = (1:2*numCellsInfinite)';
  sVars = mod((idsY - numCellsInfinite - 1) * cutslope, 1);  % Transverse variables
  % sVars = mod((idsY - numCellsInfinite - 1), 1);  % Transverse variables
  
  % Locate transverse variable in mesh
  structLoc = sLine.locateInDomain([sVars, zeros(2*numCellsInfinite, 2)]);
  elts = sLine.elements(structLoc.elements, :);
  coos = structLoc.barycoos;

  FloquetPoints = opts.FloquetPoints;
  ptsY = meshXY.points(:, 2);
  cutvec = opts.cutvec;
  period = opts.period;
  
  for idY = 1:2*numCellsInfinite
    fprintf('%d sur %d\n', idY, 2*numCellsInfinite);
    sVar = sVars(idY);
    Usol = zeros(meshXY.numPoints, opts.numCellsSemiInfinite);
    coosY = coos(idY, :);
    eltsY = elts(idY, :);
   
    for idFB = 1:opts.numFloquetPoints
      FloquetVar = FloquetPoints(idFB);
      solguide = load([nomdossier, 'half_guide_sol_', suffix, '_Floquet_', num2str(idFB)]);
      soltrace = load([nomdossier, 'sol_trace_Floquet_', num2str(idFB)]);
      solcell1 = load([nomdossier, 'local_cell_sol_', suffix, '_Floquet_', num2str(idFB), '_S_', num2str(eltsY(1))]); %#ok
      solcell2 = load([nomdossier, 'local_cell_sol_', suffix, '_Floquet_', num2str(idFB), '_S_', num2str(eltsY(2))]);


      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      FourierIdsX = spBX.FourierIds.X(:).'; dx = length(FourierIdsX);
      FourierIdsY = spBX.FourierIds.Y(:).'; dy = length(FourierIdsY);
      FourierIdsZ = spBX.FourierIds.Z(:).'; dz = length(FourierIdsZ);

      [Ix, Iy, ~] = ind2sub([dx, dy, dz], (1:spBX.numBasis)');
      r_xi_fun = FloquetVar*opts.cutvec(1) + 2*pi*(FourierIdsX(Ix) * opts.cutvec(1) + FourierIdsY(Iy) * opts.cutvec(2)).';
      l_xi_fun = sqrt(r_xi_fun.^2 - opts.omega^2);

      E0 = sinh((1 - orientation*meshXY.points(:, 1)) * l_xi_fun.') .* exp(2i*pi* meshXY.points(:, 2) * (FourierIdsX(Ix) * opts.cutvec(1) + FourierIdsY(Iy) * opts.cutvec(2))) .* exp(2i*pi*FourierIdsY(Iy)*sVar) ./ (ones(meshXY.numPoints, 1) * sinh(l_xi_fun.'));
      %
      E1 = sinh(orientation*meshXY.points(:, 1) * l_xi_fun.') .* exp(2i*pi * meshXY.points(:, 2) * (FourierIdsX(Ix) * opts.cutvec(1) + FourierIdsY(Iy) * opts.cutvec(2))) .* exp(2i*pi*FourierIdsY(Iy)*sVar) ./ (ones(meshXY.numPoints, 1) * sinh(l_xi_fun.'));
      %
      TFB_U = E0 * soltrace.R0phi + E1 * soltrace.R1phi;
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % e0x = coosY(1) * solcell1.E0x + coosY(2) * solcell2.E0x; %#ok
      % e1x = coosY(1) * solcell1.E1x + coosY(2) * solcell2.E1x;
      % e0y = coosY(1) * solcell1.E0y + coosY(2) * solcell2.E0y;
      % e1y = coosY(1) * solcell1.E1y + coosY(2) * solcell2.E1y;
      
      % TFB_U = e0x * shearmap.fun0x(sVar) * soltrace.R0phi... % X direction, E0
      %       + e1x * shearmap.fun1x(sVar) * soltrace.R1phi... % X direction, E1
      %       + e0y * shearmap.fun0y(sVar) * solguide.DtD0 * soltrace.R0phi...    % Y direction, E0
      %       + e1y * shearmap.fun1y(sVar) * solguide.DtD0 * soltrace.R0phi...    % Y direction, E0
      %       + e0y * shearmap.fun0y(sVar) * solguide.DtD1 * soltrace.R1phi...    % Y direction, E1
      %       + e1y * shearmap.fun1y(sVar) * solguide.DtD1 * soltrace.R1phi; %#ok % Y direction, E1
      % % TFB_U = e0x * ones(size(e0x, 2), size(soltrace.R0phi, 2)) % shearmap.fun0x(sVar) * soltrace.R0phi... % X direction, E0
      % %       + e1x * ones(size(e1x, 2), size(soltrace.R1phi, 2));% shearmap.fun1x(sVar) * soltrace.R1phi;%... % X direction, E1
      % %       % + e0y * shearmap.fun0y(sVar) * solguide.DtD0 * soltrace.R0phi...    % Y direction, E0
      % %       % + e1y * shearmap.fun1y(sVar) * solguide.DtD0 * soltrace.R0phi...    % Y direction, E0
      % %       % + e0y * shearmap.fun0y(sVar) * solguide.DtD1 * soltrace.R1phi...    % Y direction, E1
      % %       % + e1y * shearmap.fun1y(sVar) * solguide.DtD1 * soltrace.R1phi; %#ok % Y direction, E1
      %
      exp_floquet = exp(1i * FloquetVar * (cutvec(1) * ptsY + idY - numCellsInfinite - 1)); %#ok
      Usol = Usol + TFB_U .* exp_floquet;
    end

    Usol = Usol * W * sqrt(period / (2*pi)); %#ok

    save([nomdossier, 'solution_', suffix, '_Y_', num2str(idY)], 'Usol', '-v7.3');

  end

end
