function local_cell_script(idZmin, idZmax, cheminDonnees)

  load([cheminDonnees, '/inputs_half_guide.mat']);

  % Restrict idZmax to the number of cells
  idZmax = min(idZmax, numCellsZ);

  for idZ = idZmin:idZmax
    fprintf('%3d sur %3d\n', idZ, numCellsZ);
    solcell = struct;
    
    % Compute FE matrices
    cellZorigin = sizeCellZ * (idZ - 1) + Zorigin;

    mu_cell  = @(x)  mu_coeff([x(:, 1), x(:, 2) + cellZorigin]);
    rho_cell = @(x) rho_coeff([x(:, 1), x(:, 2) + cellZorigin]);

    FEmatrix.MM = FEPack.pdes.Form.intg(domCell, rho_cell * u * v);
    FEmatrix.KK = FEPack.pdes.Form.intg(domCell, (mu_cell * grad2(u)) * grad2(v));
    AA = FEmatrix.KK - (omega^2) * FEmatrix.MM;

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
    save([cheminDonnees, '/local_cell_sol_' suffix, '_', int2str(idZ)], '-struct', 'solcell', '-v7.3');
    % save([cheminDonnees, '/local_cell_mat_' suffix, '_', int2str(idZ)], '-struct', 'FEmatrix');
  end
  
end