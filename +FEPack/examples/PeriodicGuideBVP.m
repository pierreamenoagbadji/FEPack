function U = PeriodicGuideBVP(infiniteDirection, volBilinearIntg,...
                              mesh_negative, BoundaryStruct_negative, numCells_negative,...
                              mesh_positive, BoundaryStruct_positive, numCells_positive,...
                              opts)

  % Solve the problem in the positive half-guide
  [~, U0pos, dU0pos] = PeriodicHalfGuideBVP(mesh_positive, +1, infiniteDirection, volBilinearIntg, BoundaryStruct_positive, numCells_positive, opts);

  % Solve the problem in the negative half-guide
  [~, U0neg, dU0neg] = PeriodicHalfGuideBVP(mesh_negative, -1, infiniteDirection, volBilinearIntg, BoundaryStruct_negative, numCells_negative, opts);


end
