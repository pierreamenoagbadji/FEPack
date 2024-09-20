function construct_solution(meshXY, meshYZ, BCstruct, opts)

  % Initialization
  sLine = meshYZ.domain('xmin');
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

  %% Compute solution
  for idY = 1:2*opts.numCellsInfinite
    fprintf('%d sur %d\n', idY, 2*opts.numCellsInfinite);
    sVar = (idY - opts.numCellsInfinite - 1) * cutslope;
    Usol = zeros(meshXY.numPoints, opts.numCellsSemiInfinite);


    for idFB = 1:opts.numFloquetPoints
      FloquetVar = opts.FloquetPoints(idFB);
      solguide = load(['outputs/half_guide_sol_', opts.suffix, '_Floquet_', num2str(idFB)]);
      soltrace = load(['outputs/sol_trace_Floquet_', num2str(idFB)]);

      e0x = interpolateObj([0, mod(sVar, 1), 0], ['outputs/local_cell_sol_', opts.suffix, '_Floquet_', num2str(idFB), '_S'], 'E0x', sLine);
      e1x = interpolateObj([0, mod(sVar, 1), 0], ['outputs/local_cell_sol_', opts.suffix, '_Floquet_', num2str(idFB), '_S'], 'E1x', sLine);
      e0y = interpolateObj([0, mod(sVar, 1), 0], ['outputs/local_cell_sol_', opts.suffix, '_Floquet_', num2str(idFB), '_S'], 'E0y', sLine);
      e1y = interpolateObj([0, mod(sVar, 1), 0], ['outputs/local_cell_sol_', opts.suffix, '_Floquet_', num2str(idFB), '_S'], 'E1y', sLine);

      % solcell = load(['outputs/local_cell_sol_', opts.suffix, '_Floquet_', num2str(idFB), '_S_1.mat']); 
      % e0x = solcell.E0x;
      % e1x = solcell.E1x;
      % e0y = solcell.E0y;
      % e1y = solcell.E1y;
      
      R0phi = zeros(size(soltrace.soltrace, 1), opts.numCellsSemiInfinite);
      R0phi(:, 1) = soltrace.soltrace;
      for idX = 2:opts.numCellsSemiInfinite
        R0phi(:, idX) = solguide.R * R0phi(:, idX-1);
      end
      R1phi = solguide.D * R0phi;

      TFB_U = e0x * shearmap.fun0x(sVar) * R0phi... % X direction, E0
            + e1x * shearmap.fun1x(sVar) * R1phi... % X direction, E1
            + e0y * shearmap.fun0y(sVar) * solguide.DtD0 * R0phi...  % Y direction, E0
            + e1y * shearmap.fun1y(sVar) * solguide.DtD0 * R0phi...  % Y direction, E0
            + e0y * shearmap.fun0y(sVar) * solguide.DtD1 * R1phi...  % Y direction, E1
            + e1y * shearmap.fun1y(sVar) * solguide.DtD1 * R1phi;    % Y direction, E1
      %
      exp_floquet = exp(1i * FloquetVar * (opts.cutvec(1) * meshXY.points(:, 2) + idY - opts.numCellsInfinite - 1));

      Usol = Usol + TFB_U .* exp_floquet;
    end

    Usol = Usol * W * sqrt(opts.period / (2*pi)); %#ok

    save(['outputs/solution_', opts.suffix, '_Y_', num2str(idY)], 'Usol', '-v7.3');

  end

end
