function solve_TFB_waveguide_gen(idBEGIN, idEND, numNodes, cheminDonnees)

  % Load inputs
  load([cheminDonnees, '/inputs_', int2str(numNodes), '.mat']);

  % Restrict idEND to the number of Floquet points
  idEND = min(idEND, numFloquetPoints);
  
  % Solve the waveguide problems
  for idFB = idBEGIN:idEND
    fprintf('%d sur %d\n', idFB, idEND-idBEGIN+1);
    
    FloquetVar = FloquetPoints_pos(idFB);

    % The FE matrix is a linear combination of the elementary pieces
    AApos = mat_gradu_gradv_pos + 1i * FloquetVar * mat_vec1u_gradv_pos - 1i * FloquetVar * mat_gradu_vec1v_pos + FloquetVar * FloquetVar * mat_vec1u_vec1v_pos - (omega^2) * mat_u_v_pos;

    % Compute the Floquet-Bloch transform of the solution
    [~, sol_pos_data.E0, sol_pos_data.E1, sol_pos_data.R, sol_pos_data.D, ~, TFBlambdaPos] = PeriodicHalfGuideBVP(mesh3Dpos, +1, semiInfiniteDirection, AApos, BCstruct_pos, numCellsSemiInfinite_pos, opts);
    
    % Save the Floquet-Bloch transform of the solution
    save([cheminDonnees, '/sol_pos_data_', int2str(idFB), '.mat'], '-struct', 'sol_pos_data', '-v7.3');
    save([cheminDonnees, '/TFBlambdaPos_', int2str(idFB), '.mat'], 'TFBlambdaPos', '-v7.3');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    FloquetVar = FloquetPoints_neg(idFB);

    % The FE matrix is a linear combination of the elementary pieces
    AAneg = mat_gradu_gradv_neg + 1i * FloquetVar * mat_vec1u_gradv_neg - 1i * FloquetVar * mat_gradu_vec1v_neg + FloquetVar * FloquetVar * mat_vec1u_vec1v_neg - (omega^2) * mat_u_v_neg;

    % Compute the Floquet-Bloch transform of the solution
    [~, sol_neg_data.E0, sol_neg_data.E1, sol_neg_data.R, sol_neg_data.D, ~, TFBlambdaNeg] = PeriodicHalfGuideBVP(mesh3Dneg, +1, semiInfiniteDirection, AAneg, BCstruct_neg, numCellsSemiInfinite_neg, opts);
    
    % Save the Floquet-Bloch transform of the solution
    save([cheminDonnees, '/sol_neg_data_', int2str(idFB), '.mat'], '-struct', 'sol_neg_data', '-v7.3');
    save([cheminDonnees, '/TFBlambdaNeg_', int2str(idFB), '.mat'], 'TFBlambdaNeg', '-v7.3');

  end
  
end