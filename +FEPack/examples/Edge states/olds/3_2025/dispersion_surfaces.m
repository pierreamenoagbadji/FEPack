function eigvals = dispersion_surfaces(Kpoints, numEigvals, ...
                                       vecPer1, vecPer2, dualVec1, dualVec2, ...
                                       potential, ...
                                       numNodesX, numNodesY,...
                                       use_parallel, numCores)
%% eigvals = DISPERSION_SURFACES

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
Tmat = [dualVec1'; dualVec2'] / (2*pi); % Inverse of Rmat
Te1 = Tmat' * [1; 0];
Te2 = Tmat' * [0; 1];
transformed_potential = @(x) potential((Rmat * x(:, 1:2)')');

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
mat_trPot_u_v   = FEPack.pdes.Form.intg(cellXY, (transformed_potential * u) * v);
mat_u_v = FEPack.pdes.Form.intg(cellXY, u * v);
fprintf('FE matrices computation done.\n');

% Periodic boundary conditions
ecs = (((u|edge0x) - (u|edge1x)) == 0.0) & (((u|edge0y) - (u|edge1y)) == 0.0);
ecs.applyEcs;
PP = ecs.P;

%% Compute dispersion surfaces
numK = size(Kpoints, 1);
eigvals = zeros(numK, numEigvals);

if (use_parallel)

  if (nargin < 7)
    error('use_parallel has been set to true; you need to specify the number of cores. Use the command feature(''numcores'') to see the number of available cores.');
  end

  % Parallel version
  pool = gcp('nocreate'); % Get current parallel pool
  if isempty(pool)
    parpool('local', numCores); % Start a parallel pool if none exists
  end
  fprintf('Using parfor (parallel execution)\n');

  parfor idI = 1:numK

    fprintf('%d sur %d\n', idI, numK);
      
    kv = Kpoints(idI, :);
    kx = kv * vecPer1;
    ky = kv * vecPer2;

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
       + mat_trPot_u_v;

    BB = mat_u_v;

    % Apply essential conditions
    AA0 = PP * AA * PP';
    BB0 = PP * BB * PP';

    % Solve eigenvalue problem
    [~, D] = eigs(AA0, BB0, numEigvals, 'smallestabs');
    eigvals(idI, :) = diag(D).';

  end

else

  % Sequential execution
  fprintf('Using for (sequential execution)\n');

  for idI = 1:numK

    fprintf('%d sur %d\n', idI, numK);
      
    kv = Kpoints(idI, :);
    kx = kv * vecPer1;
    ky = kv * vecPer2;

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
      + mat_trPot_u_v;

    BB = mat_u_v;

    % Apply essential conditions
    AA0 = PP * AA * PP';
    BB0 = PP * BB * PP';

    % Solve eigenvalue problem
    [~, D] = eigs(AA0, BB0, numEigvals, 'smallestreal');
    eigvals(idI, :) = diag(D).';

  end
  
end

