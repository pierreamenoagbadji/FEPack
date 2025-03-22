function [Edirac, minEgy, maxEgy, phi1K] = compute_dirac_range(...
                                            numEigvals,...
                                            vecPer1, vecPer2, dualVec1, dualVec2,...
                                            Vpot, Wpot, delta,...
                                            edge_vec1, ~, edge_dual_vec1, edge_dual_vec2,...
                                            numNodesX, numNodesY, numNodesK, is_rational,...
                                            use_parallel, numcores, tol, plot_dispersion_functions)

  if (nargin < 20)
    tol = 1e-1;
  end

  %% Generate mesh
  meshXY = FEPack.meshes.MeshRectangle(1, [0 1], [0 1], numNodesX, numNodesY);
  cellXY = meshXY.domain('volumic');
  edge0x = meshXY.domain('xmin');
  edge1x = meshXY.domain('xmax');
  edge0y = meshXY.domain('ymin');
  edge1y = meshXY.domain('ymax');

  u = FEPack.pdes.PDEObject; 
  v = dual(u);

  Rmat = [vecPer1,   vecPer2];
  Tmat = [dualVec1'; dualVec2'] / (2*pi);
  % Tmat = Rmat \ eye(2);
  Te1 = Tmat' * [1; 0];
  Te2 = Tmat' * [0; 1];
  V_trs_pot = @(x) Vpot((Rmat * x(:, 1:2)')');
  W_trs_pot = @(x) Wpot((Rmat * x(:, 1:2)')');

  %% Finite elements matrices
  fprintf('FE matrices computation\n');
  mat_gradu_gradv = FEPack.pdes.Form.intg(cellXY, (Tmat' * grad2(u))   * (Tmat' * grad2(v)));
  mat_gradu_vec1v = FEPack.pdes.Form.intg(cellXY, (Tmat' * grad2(u))   * (Te1   * v));
  mat_gradu_vec2v = FEPack.pdes.Form.intg(cellXY, (Tmat' * grad2(u))   * (Te2   * v));
  mat_vec1u_gradv = FEPack.pdes.Form.intg(cellXY, (Te1   * u) * (Tmat' * grad2(v)));
  mat_vec2u_gradv = FEPack.pdes.Form.intg(cellXY, (Te2   * u) * (Tmat' * grad2(v)));
  mat_vec1u_vec1v = FEPack.pdes.Form.intg(cellXY, (Te1   * u) * (Te1   * v));
  mat_vec1u_vec2v = FEPack.pdes.Form.intg(cellXY, (Te1   * u) * (Te2   * v));
  mat_vec2u_vec1v = FEPack.pdes.Form.intg(cellXY, (Te2   * u) * (Te1   * v));
  mat_vec2u_vec2v = FEPack.pdes.Form.intg(cellXY, (Te2   * u) * (Te2   * v));
  mat_VtrPot_u_v  = FEPack.pdes.Form.intg(cellXY, (V_trs_pot * u) * v);
  mat_WtrPot_u_v  = FEPack.pdes.Form.intg(cellXY, (W_trs_pot * u) * v);
  mat_u_v = FEPack.pdes.Form.intg(cellXY, u * v);
  fprintf('FE matrices computation done.\n');

  % Periodic boundary conditions
  ecs = (((u|edge0x) - (u|edge1x)) == 0.0) & (((u|edge0y) - (u|edge1y)) == 0.0);
  ecs.applyEcs;
  PP = ecs.P;

  %% Compute dispersion surfaces at K point
  highSymK = (dualVec1 - dualVec2) / 3;
  kx = highSymK' * vecPer1;
  ky = highSymK' * vecPer2;

  % FE matrices
  AA = mat_gradu_gradv +...
     + 1i * kx * mat_vec1u_gradv...
     + 1i * ky * mat_vec2u_gradv...
     - 1i * kx * mat_gradu_vec1v...
     - 1i * ky * mat_gradu_vec2v...
     + (kx * kx) * mat_vec1u_vec1v...
     + (kx * ky) * mat_vec1u_vec2v...
     + (ky * kx) * mat_vec2u_vec1v...
     + (ky * ky) * mat_vec2u_vec2v...
     + mat_VtrPot_u_v;

  BB = mat_u_v;

  % Apply essential conditions
  AA0 = PP * AA * PP';
  BB0 = PP * BB * PP';

  % Solve eigenvalue problem
  [V, D] = eigs(AA0, BB0, numEigvals, 'smallestabs');
  eigvals = diag(D).';

  %% Find Dirac point
  % The dispersion functions with the smallest difference (which should be 0 in theory)
  % are the one we consider as intersecting
  diff_disp = abs(eigvals(2:end) - eigvals(1:end-1));
  [ecartmin, Imin] = min(diff_disp); 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  Imin = 1;     % À SUPPRIMER
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  Edirac = eigvals(Imin);
  disp(ecartmin);

  if (ecartmin > tol/numNodesX)
    warning(['The closest dispersion functions still have a difference of ', num2str(ecartmin, '%0.3e'), '.']);
  end

  if (imag(Edirac) > tol/numNodesX)
    warning(['The Dirac point has an imaginary part: ', num2str(imag(Edirac), '%0.3e'), '.']);
  end

  Edirac = real(Edirac);

  %% Dirac eigenbasis
  % Phi = a_+ phi_+(x, K) + a_- phi_-(x, K) where
  % a_+ and a_- are chosen to minimize the difference
  % Phi(R^* x) - Phi(x), R being the clockwise rotation matrix by 2pi/3
  phi_m = V(:, Imin);
  phi_p = V(:, Imin+1);

  % rotation matrix
  rotmat = [-0.5, 0.5*sqrt(3); -0.5*sqrt(3), -0.5];
  rotmat = [rotmat, [0; 0]; 0 0 0]; % Padding with 0s

  % Compute phi_± (R^* x)
  structLoc = cellXY.locateInDomain(mod(rotmat.' * meshXY.points.', 1).');
  elts = cellXY.elements(structLoc.elements, :);
  coos = structLoc.barycoos;

  phi_m_extended = PP' * phi_m;
  phi_p_extended = PP' * phi_p;

  Rphi_m = PP * diag(coos * phi_m_extended(elts).');
  Rphi_p = PP * diag(coos * phi_p_extended(elts).');

  % The most optimal combination, that is the vector [a- a+]
  % such that Phi(R^* x) - Phi(x) is the smallest with
  % Phi = a_+ phi_+(x, K) + a_- phi_-(x, K), is obtained from 
  % the eigenvector associated to the smallest eigenvalue of
  % the Gram matrix 
  gram.mat = [Rphi_m' * BB0 * Rphi_m, Rphi_p' * BB0 * Rphi_m;...
              Rphi_m' * BB0 * Rphi_p, Rphi_p' * BB0 * Rphi_p];

  [gram.V, gram.D] = eig(gram.mat);

  [gram.minD, gram.Imin] = min(diag(gram.D));

  % phi_1 of the Dirac eigenbasis
  phi1K = gram.V(1, gram.Imin) * phi_m + gram.V(2, gram.Imin) * phi_p;

  % % Compute upsilon
  % tr_W = @(x) Wpot((Rmat * x(:, 1:2)')');
  % mat_W_u_v = FEPack.pdes.Form.intg(cellXY, (tr_W * u) * v);
  % upsilon = det(Tmat) * phi1K' * (PP * mat_W_u_v * PP') * phi1K / (phi1K' * BB0 * phi1K);
  % upsilon = real(upsilon);

  %% Range around Dirac point from spectrum of -Δ + V + δ*W
  mat_gradu_vec2v = FEPack.pdes.Form.intg(cellXY, (Tmat' * grad2(u))   * (Te2   * v));
  mat_vec2u_gradv = FEPack.pdes.Form.intg(cellXY, (Te2   * u) * (Tmat' * grad2(v)));
  mat_vec1u_vec2v = FEPack.pdes.Form.intg(cellXY, (Te1   * u) * (Te2   * v));
  mat_vec2u_vec1v = FEPack.pdes.Form.intg(cellXY, (Te2   * u) * (Te1   * v));
  mat_vec2u_vec2v = FEPack.pdes.Form.intg(cellXY, (Te2   * u) * (Te2   * v));

  if (is_rational)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    lamb = linspace(-0.01, 0.01, numNodesK);
    Kpoints = highSymK + (-0.01*edge_dual_vec1) + edge_dual_vec2 * lamb;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  else

    meshK = FEPack.meshes.MeshRectangle(1, [0 1], [0 1], numNodesK, numNodesK);
    Kpoints = edge_dual_vec1 * meshK.points(:, 1).'...
            + edge_dual_vec2 * meshK.points(:, 2).';

  end

  numKpoints = size(Kpoints, 2);
  disp_pert = zeros(numKpoints, numEigvals);

  if (use_parallel)

    if (nargin < 18)
      error('use_parallel has been set to true; you need to specify the number of cores. Use the command feature(''numcores'') to see the number of available cores.');
    end

    % Parallel version
    pool = gcp('nocreate'); % Get current parallel pool
    if isempty(pool)
      parpool('local', numcores); % Start a parallel pool if none exists
    end
    fprintf('Using parfor (parallel execution)\n');

    parfor idK = 1:numKpoints
      fprintf('%d sur %d\n', idK, numKpoints);
      
      kv = Kpoints(:, idK);
      kx = kv' * vecPer1;
      ky = kv' * vecPer2;

      % FE matrices
      AA = mat_gradu_gradv +...
          + 1i * kx * mat_vec1u_gradv...
          + 1i * ky * mat_vec2u_gradv...
          - 1i * kx * mat_gradu_vec1v...
          - 1i * ky * mat_gradu_vec2v...
          + (kx * kx) * mat_vec1u_vec1v...
          + (kx * ky) * mat_vec1u_vec2v...
          + (ky * kx) * mat_vec2u_vec1v...
          + (ky * ky) * mat_vec2u_vec2v...
          +         mat_VtrPot_u_v...
          + delta * mat_WtrPot_u_v;

      BB = mat_u_v;

      AA0 = PP * AA * PP';
      BB0 = PP * BB * PP';

      [~, D] = eigs(AA0, BB0, numEigvals, 'smallestabs');
      disp_pert(idK, :) = diag(D).';
    end

  else

    for idK = 1:numKpoints
      fprintf('%d sur %d\n', idK, numKpoints);
  
      kv = Kpoints(:, idK);
      kx = kv' * vecPer1;
      ky = kv' * vecPer2;

      % FE matrices
      AA = mat_gradu_gradv +...
          + 1i * kx * mat_vec1u_gradv...
          + 1i * ky * mat_vec2u_gradv...
          - 1i * kx * mat_gradu_vec1v...
          - 1i * ky * mat_gradu_vec2v...
          + (kx * kx) * mat_vec1u_vec1v...
          + (kx * ky) * mat_vec1u_vec2v...
          + (ky * kx) * mat_vec2u_vec1v...
          + (ky * ky) * mat_vec2u_vec2v...
          +         mat_VtrPot_u_v...
          + delta * mat_WtrPot_u_v;

      BB = mat_u_v;

      AA0 = PP * AA * PP';
      BB0 = PP * BB * PP';

      [~, D] = eigs(AA0, BB0, numEigvals, 'smallestabs');
      disp_pert(idK, :) = diag(D).';
    end

  end

  disp_pert = real(disp_pert);
  minEgy = max(disp_pert(:, Imin));
  maxEgy = min(disp_pert(:, Imin+1));

  if (plot_dispersion_functions)

    if (is_rational)

      figure;
      set(groot,'defaultAxesTickLabelInterpreter','latex');
      set(groot,'defaulttextinterpreter','latex');
      set(groot,'defaultLegendInterpreter','latex');

      for idE = 1:numEigvals
        plot(lamb, disp_pert(:, idE), 'b');
        hold on;
      end
      plot(lamb, disp_pert(:, Imin), 'r');
      plot(lamb, disp_pert(:, Imin+1), 'r');

      % plot([-0.5 0.5], [Edirac, Edirac], 'r--');
      plot([-0.5 0.5], [maxEgy, maxEgy], 'r--');
      plot([-0.5 0.5], [minEgy, minEgy], 'r--');
      
      plot([highSymK.' * edge_vec1 highSymK.' * edge_vec1], [min(min(disp_pert)) max(max(disp_pert))], 'r--');
      xlim([-0.5 0.5]);
      ylim([min(min(disp_pert)) max(max(disp_pert))]);
      xlabel('$\lambda$');
      ylabel('Eigenvalues for $V + \delta W$');
    
    else

      warning('Not coded yet.');

    end
    [minEgy, maxEgy]
    pause;

  end
end
