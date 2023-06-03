function solve_TFB_waveguide(idFB)

  load('mat');

  FloquetVar = FloquetPoints(idFB);

  % The FE matrix is a linear combination of the elementary pieces
  AApos = mat_gradu_gradv_pos + 1i * FloquetVar * mat_vec1u_gradv_pos - 1i * FloquetVar * mat_gradu_vec1v_pos + FloquetVar * FloquetVar * mat_vec1u_vec1v_pos + mat_u_v_pos;
  AAneg = mat_gradu_gradv_neg + 1i * FloquetVar * mat_vec1u_gradv_neg - 1i * FloquetVar * mat_gradu_vec1v_neg + FloquetVar * FloquetVar * mat_vec1u_vec1v_neg + mat_u_v_neg;

  % The Floquet-Bloch transform of the boundary data
  jumpData_FB = @(x) BlochTransform(x, FloquetVar, G3D, infiniteDirection);

  % Compute the Floquet-Bloch transform of the solution
  TFBU = PeriodicGuideJumpBVP(semiInfiniteDirection,...
                              AApos, mesh3Dpos, BCstruct_pos, numCellsSemiInfinite_pos,...
                              AAneg, mesh3Dneg, BCstruct_neg, numCellsSemiInfinite_neg,...
                              jumpData_FB, opts);
  
  % Save the Floquet-Bloch transform of the solution
  save(['outputs/TFBU_', int2str(idFB), '.mat'], '-struct', 'TFBU');
  
end