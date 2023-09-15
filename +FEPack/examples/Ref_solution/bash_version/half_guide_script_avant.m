function half_guide_script_avant(guide, cheminDonnees)
  load([cheminDonnees, '/inputs.mat']);
  
  if strcmpi(guide, 'pos')

    suffix = 'pos';
    meshcell = mesh_pos;
    orientation = 1;
    omega = opts.omega;
    mu_coeff = mu_pos;
    rho_coeff = rho_pos;

  elseif strcmpi(guide, 'neg')

    suffix = 'neg';
    meshcell = mesh_neg;
    orientation = -1;
    omega = opts.omega;
    mu_coeff = mu_neg;
    rho_coeff = rho_neg;
    
  else

    error(['argument', guide, ' non reconnu.']);

  end

  save([cheminDonnees, '/inputs_half_guide.mat'], 'meshcell', 'orientation', 'omega', 'mu_coeff', 'rho_coeff', 'numCellsZ', 'sizeCellZ', 'Zorigin', 'suffix', 'basis_function', '-v7.3');
  %
  clearvars -except cheminDonnees; % Clear all variables except cheminDonnees
  %
  load([cheminDonnees, '/inputs_half_guide.mat']);

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

  save([cheminDonnees, '/inputs_half_guide.mat'], '-v7.3');
end