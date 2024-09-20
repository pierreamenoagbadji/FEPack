function compute_guide_solution_script(idZmin, idZmax, cheminDonnees)

  load([cheminDonnees, '/inputs.mat']);

  % Restrict idZmax to the number of cells
  idZmax = min(idZmax, numCellsZ);

  for idZ = idZmin:idZmax

    fprintf('%3d sur %3d\n', idZ, numCellsZ);
  
    % Solution in the positive guide
    R0phi = solphi;
    R1phi = solGuidePos.Sop * solphi;
    SC = load([cheminDonnees, '/local_cell_sol_pos_', int2str(idZ)]);
    idBasis = 1+(idZ-1)*(N0x-1):1+idZ*(N0x-1);

    for idX = 1:numCellsXpos
      idX

      Usol = (SC.E0x * solGuidePos.B0(idBasis, :) +...
              SC.E0z * solGuidePos.traceE0{idZ} +...
              SC.E1z * solGuidePos.traceE0{idZ+1}) * R0phi...
              ...
           + (SC.E1x * solGuidePos.B1(idBasis, :) +...
              SC.E0z * solGuidePos.traceE1{idZ} +...
              SC.E1z * solGuidePos.traceE1{idZ+1}) * R1phi;
      
      % Save output
      save([cheminDonnees, '/guide_sol_pos_', int2str(idX), '_', int2str(idZ)], 'Usol', '-v7.3');

      % Update the coefficients
      R0phi = solGuidePos.Pop * R0phi;
      R1phi = solGuidePos.Sop * R0phi;
    end

    % Solution in the negative guide
    R0phi = solphi;
    R1phi = solGuideNeg.Sop * solphi;
    SC = load([cheminDonnees, '/local_cell_sol_neg_', int2str(idZ)]);
    idBasis = 1+(idZ-1)*(N0x-1):1+idZ*(N0x-1);

    for idX = 1:numCellsXneg
      idX
      
      Usol = (SC.E0x * solGuideNeg.B0(idBasis, :) +...
              SC.E0z * solGuideNeg.traceE0{idZ} +...
              SC.E1z * solGuideNeg.traceE0{idZ+1}) * R0phi...
              ...
           + (SC.E1x * solGuideNeg.B1(idBasis, :) +...
              SC.E0z * solGuideNeg.traceE1{idZ} +...
              SC.E1z * solGuideNeg.traceE1{idZ+1}) * R1phi;
      
      % Save outputs
      save([cheminDonnees, '/guide_sol_neg_', int2str(idX), '_', int2str(idZ)], 'Usol', '-v7.3');

      % Update the coefficients
      R0phi = solGuideNeg.Pop * R0phi;
      R1phi = solGuideNeg.Sop * R0phi;
    end


  end

end