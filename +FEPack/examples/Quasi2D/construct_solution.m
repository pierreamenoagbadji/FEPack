function construct_solution(meshXY, meshLineZ, BCstruct, opts)

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

  % Shear maps
  shearmap.fun0x = @(s) spBX.evaluateBasisFunctions([    opts.cutvec(1) * meshXY.points(edge0x.IdPoints, 2),...
                                                     mod(opts.cutvec(2) * meshXY.points(edge0x.IdPoints, 2) + s, 1), zeros(N0x, 1)]);
  %
  shearmap.fun1x = @(s) spBX.evaluateBasisFunctions([    opts.cutvec(1) * meshXY.points(edge1x.IdPoints, 2),...
                                                     mod(opts.cutvec(2) * meshXY.points(edge1x.IdPoints, 2) + s, 1), zeros(N1x, 1)]);
  %
  shearmap.fun0y = @(s) spBY.evaluateBasisFunctions([mod(           s * ones(N0y-2, 1), 1), meshXY.points(edge0y.IdPoints(2:end-1), 1), zeros(N0y-2, 1)]);
  shearmap.fun1y = @(s) spBY.evaluateBasisFunctions([mod(cutslope + s * ones(N0y-2, 1), 1), meshXY.points(edge1y.IdPoints(2:end-1), 1), zeros(N1y-2, 1)]);
  W = (2*pi/opts.period) ./ (opts.numFloquetPoints - 1);

  %% Traces of FB transforms
  for idFB = 1:opts.numFloquetPoints
    soltrace = load(['outputs/sol_trace_Floquet_', num2str(idFB)]);
    solguide = load(['outputs/half_guide_sol_', opts.suffix, '_Floquet_', num2str(idFB)]);
    
    R0phi = zeros(size(soltrace.soltrace, 1), opts.numCellsSemiInfinite);
    R0phi(:, 1) = soltrace.soltrace;

    for idX = 2:opts.numCellsSemiInfinite
      R0phi(:, idX) = solguide.R * R0phi(:, idX-1);
    end
    R1phi = solguide.D * R0phi; %#ok
    
    save(['outputs/sol_trace_Floquet_', num2str(idFB)], 'R0phi', 'R1phi', '-append');
  end

  %% Compute inverse FB transform
  idsY = (1:2*opts.numCellsInfinite)';
  sVars = mod((idsY - opts.numCellsInfinite - 1) * cutslope, 1);  % Transverse variables
  
  % Locate transverse variable in mesh
  structLoc = sLine.locateInDomain([sVars, zeros(2*opts.numCellsInfinite, 2)]);
  elts = sLine.elements(structLoc.elements, :);
  coos = structLoc.barycoos;
  
  for idY = 1:2*opts.numCellsInfinite
    fprintf('%d sur %d\n', idY, 2*opts.numCellsInfinite);
    sVar = sVars(idY);
    Usol = zeros(meshXY.numPoints, opts.numCellsSemiInfinite);
   
    for idFB = 1:opts.numFloquetPoints
      FloquetVar = opts.FloquetPoints(idFB);
      solguide = load(['outputs/half_guide_sol_', opts.suffix, '_Floquet_', num2str(idFB)]);
      soltrace = load(['outputs/sol_trace_Floquet_', num2str(idFB)]);
      solcell1 = load(['outputs/local_cell_sol_', opts.suffix, '_Floquet_', num2str(idFB), '_S_', num2str(elts(idY, 1))]);
      solcell2 = load(['outputs/local_cell_sol_', opts.suffix, '_Floquet_', num2str(idFB), '_S_', num2str(elts(idY, 2))]);

      e0x = coos(idY, 1) * solcell1.E0x + coos(idY, 2) * solcell2.E0x;
      e1x = coos(idY, 1) * solcell1.E1x + coos(idY, 2) * solcell2.E1x;
      e0y = coos(idY, 1) * solcell1.E0y + coos(idY, 2) * solcell2.E0y;
      e1y = coos(idY, 1) * solcell1.E1y + coos(idY, 2) * solcell2.E1y;
      
      TFB_U = e0x * shearmap.fun0x(sVar) * soltrace.R0phi... % X direction, E0
            + e1x * shearmap.fun1x(sVar) * soltrace.R1phi... % X direction, E1
            + e0y * shearmap.fun0y(sVar) * solguide.DtD0 * soltrace.R0phi...  % Y direction, E0
            + e1y * shearmap.fun1y(sVar) * solguide.DtD0 * soltrace.R0phi...  % Y direction, E0
            + e0y * shearmap.fun0y(sVar) * solguide.DtD1 * soltrace.R1phi...  % Y direction, E1
            + e1y * shearmap.fun1y(sVar) * solguide.DtD1 * soltrace.R1phi;    % Y direction, E1
      %
      exp_floquet = exp(1i * FloquetVar * (opts.cutvec(1) * meshXY.points(:, 2) + idY - opts.numCellsInfinite - 1));
      Usol = Usol + TFB_U .* exp_floquet;
    end

    Usol = Usol * W * sqrt(opts.period / (2*pi)); %#ok

    save(['outputs/solution_', opts.suffix, '_Y_', num2str(idY)], 'Usol', '-v7.3');

  end

end
