function init_script_transmission_gen(numFloquetPoints, numNodes, cheminDonnees)
  
  if (nargin < 1)
    numFloquetPoints = 50;
  end

  if (nargin < 2)
    numNodes = 10;
  end

  import FEPack.*

  %% Problem-related variables
  % //////////////////////////
  opts.omega = 8 + 0.25i;
  omega = opts.omega;
  period = 1;
  Lint = 5;
  opts.verbose = 0;

  % 2D coefficients
  vecperfunpos = [-0.5*sqrt(2), 1];
  vecperfunneg = [-pi, 1];% [-pi, 1];

  mu2Dpos = @(x) 0.5 + perCutoffCircle(x, [1; 0], vecperfunpos, [0.5, 0.5], [-0.2, 0.2]);
  rho2Dpos = @(x) 0.5 + perCutoffCuboid(x, [1; 0], vecperfunpos, [0.5, 0.5], [-0.2, 0.2], [-0.2, 0.2], 0.5, 1);
  mu2Dneg  = @(x) ones(size(x, 1), 1); % 1 + 0.5 * cos(2*pi*(x(:, 1) - x(:, 2)*vecperfunneg(1)/vecperfunneg(2))) .* cos(2*pi*x(:, 2)/vecperfunneg(2));
  rho2Dneg = @(x) ones(size(x, 1), 1); % 1 + 0.25 * sin(2*pi*(x(:, 1) - x(:, 2)*vecperfunneg(1)/vecperfunneg(2))) + 0.25 * sin(2*pi*x(:, 2)/vecperfunneg(2));

  % Cut matrices
  vecperpos = vecperfunpos;
  vecperneg = vecperfunneg;

  cutvecpos = [1/vecperpos(2); -vecperpos(1)/vecperpos(2)];
  cutvecneg = [1/vecperneg(2); -vecperneg(1)/vecperneg(2)];

  opts.cutmat_pos = [[1; 0; 0], [0; cutvecpos]];
  opts.cutmat_neg = [[1; 0; 0], [0; cutvecneg]];

  cutslopePos = cutvecpos(2) / cutvecpos(1);
  cutslopeNeg = cutvecneg(2) / cutvecneg(1);

  % 3D coefficients
  fun3D = @(fun2D, vecper, x) fun2D([x(:, 1) + x(:, 3) + x(:, 2) * vecper(1), x(:, 2) * vecper(2)]);
  mu3Dpos  = @(x) fun3D(mu2Dpos,  vecperpos, x);
  mu3Dneg  = @(x) fun3D(mu2Dneg,  vecperneg, x);
  rho3Dpos = @(x) fun3D(rho2Dpos, vecperpos, x);
  rho3Dneg = @(x) fun3D(rho2Dneg, vecperneg, x);

  % Jump data
  Gint = @(x) 100 * FEPack.tools.cutoff(x(:, 1), -0.5, 0.5);
  Gtrs = @(s) ones(size(s, 1), 1);

  % (semi-)infinite directions and numbers of cells
  semiInfiniteDirection = 1;
  infiniteDirection = 2;
  numCellsSemiInfinite_pos = 6;
  numCellsSemiInfinite_neg = 6;
  numCellsInfinite = 5;
  numFloquetPoints = 20;

  %% Mesh
  pregenerate_mesh = 1;
  struct_mesh = 1;
  numNodes2D = 10;
  numNodes3D = numNodes;

  if pregenerate_mesh

    if (struct_mesh)
      mesh_prefix = 'struct';
    else
      mesh_prefix = 'unstruct';
    end

    % Pick mesh from saved file
    m = load(['pregenMeshes/2D/', mesh_prefix, '_mesh_2D_', int2str(numNodes2D), '_positive.mat']);
    mesh2Dpos = m.mesh;

    m = load(['pregenMeshes/3D/', mesh_prefix, '_mesh_3D_', int2str(numNodes3D), '_positive.mat']);
    mesh3Dpos = m.mesh;

    m = load(['pregenMeshes/2D/', mesh_prefix, '_mesh_2D_', int2str(numNodes2D), '_negative_X.mat']);
    mesh2Dneg = m.mesh;

    m = load(['pregenMeshes/3D/', mesh_prefix, '_mesh_3D_', int2str(numNodes3D), '_negative_X.mat']);
    mesh3Dneg = m.mesh;

  else

    mesh2Dpos = meshes.MeshRectangle(struct_mesh, [0 1], [0 1], numNodes2D, numNodes2D);
    mesh3Dpos = meshes.MeshCuboid(struct_mesh, [0 1], [0 1], [0 1], numNodes3D, numNodes3D, numNodes3D);

    mesh2Dneg = meshes.MeshRectangle(struct_mesh, [0 -1], [0 1], numNodes2D, numNodes2D);
    mesh3Dneg = meshes.MeshCuboid(struct_mesh, [0 -1], [0 1], [0 1], numNodes3D, numNodes3D, numNodes3D);

  end

  % Mesh for the interface
  numNodes_int = numNodes2D/4;
  mesh2Dint = meshes.MeshSegment('uniform', -Lint, Lint, floor(2*Lint*numNodes_int));
  mesh2Dtrs = meshes.MeshSegment('uniform', 0, 1, numNodes2D);

  %% Basis functions associated to the interface
  % ////////////////////////////////////////////
  type_basis_int = 'Lagrange'; % 'sine', 'Lagrange'
  type_basis_trs = 'Fourier'; % 'Fourier'

  if strcmpi(type_basis_int, 'sine')

    phisInt = @(x, n) sqrt(1 / Lint) * sin(pi * (x(:, 1) + Lint) * n / (2 * Lint)) .* (abs(x(:, 1)) <= Lint);
    numBasisInt = floor(mesh2Dint.numPoints / 2);

  elseif strcmpi(type_basis_int, 'Lagrange')

    % error('Il n''y a pas de bug. Je veux juste que tu vérifies le nombre de noeuds (pas trop grand).');
    phisInt = @(x, n) ((x(:, 1) - mesh2Dint.points(n, 1).') ./ (mesh2Dint.points(n+1, 1).' - mesh2Dint.points(n, 1).'))...
                                    .* (x(:, 1) >= mesh2Dint.points(n, 1).' & x(:, 1) < mesh2Dint.points(n+1, 1).')...
                                + ((mesh2Dint.points(n+2, 1).' - x(:, 1)) ./ (mesh2Dint.points(n+2, 1).' - mesh2Dint.points(n+1, 1).'))...
                                    .* (x(:, 1) >= mesh2Dint.points(n+1, 1).' & x(:, 1) < mesh2Dint.points(n+2, 1).');
    numBasisInt = mesh2Dint.numPoints - 2;

  else

    error(['Type de fonction de base', type_basis_int, 'non reconnu.']);

  end
  spBint = spaces.SpectralBasis(mesh2Dint.domain('volumic'), phisInt, numBasisInt);
  spBint.computeBasisMatrices;

  numBasisTrs = floor(numNodes2D/4);
  if strcmpi(type_basis_trs, 'Fourier')

    spBtrs = spaces.FourierBasis(mesh2Dtrs.domain('volumic'), numBasisTrs);

  else

    error(['Type de fonction de base', type_basis_trs, 'non reconnu.']);

  end

  numBasis = spBint.numBasis * spBtrs.numBasis;
  [idInt, idTrs] = ind2sub([spBint.numBasis, spBtrs.numBasis], 1:numBasis);

  phis3Dpos = @(x, idI) spBtrs.phis(x(:, 3) - x(:, 2) * cutslopePos, idTrs(idI)) .* phisInt(x(:, 2) / cutvecpos(1), idInt(idI));
  phis3Dneg = @(x, idI) spBtrs.phis(x(:, 3) - x(:, 2) * cutslopeNeg, idTrs(idI)) .* phisInt(x(:, 2) / cutvecneg(1), idInt(idI));

  %%
  % phis3Dtest = @(x, idI) spBtrs.phis(x(:, 2) - x(:, 1) * cutslopePos, idTrs(idI)) .* phisInt(x(:, 1) / cutvecpos(1), idInt(idI));
  % meshPhi = FEPack.meshes.MeshRectangle(1, [-2 2], [0 1], 40, 40);
  % for idI = 1:numBasis
  %   trisurf(meshPhi.triangles, meshPhi.points(:, 1), meshPhi.points(:, 2), real(phis3Dtest(meshPhi.points, idI)));
  %   view(2);
  %   shading interp
  %   set(gca, 'DataAspectRatio', [1 1 1])
  %   pause;
  % end

  %% Plot coefficients?
  plot_coefficients = false;
  if (plot_coefficients)
    set(groot,'defaultAxesTickLabelInterpreter','latex');
    set(groot,'defaulttextinterpreter','latex');
    set(groot,'defaultLegendInterpreter','latex');

    mu2Dglob =  @(x) (x(:, 1) >= 0) .*  mu2Dpos(x) + (x(:, 1) <  0) .* mu2Dneg(x);
    rho2Dglob = @(x) (x(:, 1) >= 0) .* rho2Dpos(x) + (x(:, 1) <  0) .* rho2Dneg(x);

    for idS = 1:(2*numCellsSemiInfinite_pos)
      for idI = 1:(2*numCellsInfinite)
        X = mesh2Dpos.points(:, 1) + (idS - numCellsSemiInfinite_pos - 1);
        Y = mesh2Dpos.points(:, 2) + (idI - numCellsInfinite - 1);
        figure(1);
        trisurf(mesh2Dpos.triangles, X, Y, mu2Dglob([X, Y]));
        hold on;
        view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
        set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);

        figure(2);
        trisurf(mesh2Dpos.triangles, X, Y, rho2Dglob([X, Y]));
        hold on;
        view(2); shading interp; colorbar('TickLabelInterpreter', 'latex');
        set(gca,'DataAspectRatio',[1 1 1], 'FontSize', 16);

      end
    end

    figure(1);
    xlim([-numCellsSemiInfinite_pos, numCellsSemiInfinite_pos]);
    ylim([-numCellsInfinite, numCellsInfinite]);

    figure(2);
    xlim([-numCellsSemiInfinite_pos, numCellsSemiInfinite_pos]);
    ylim([-numCellsInfinite, numCellsInfinite]);

    figure(3);
    x = [zeros(256, 1), linspace(-2, 2, 256)', zeros(256, 1)];
    plot(x(:, 2), Gint(x));
    
  end

  %% Bilinear and linear forms
  % //////////////////////////
  u = pdes.PDEObject; v = dual(u);
  gradu_gradv = @(cutmat, muco) (muco * (cutmat' * grad3(u))) * (cutmat' * grad3(v));
  gradu_vec1v = @(cutmat, muco) (muco * (cutmat' * grad3(u))) * (cutmat' * [0; 1; 0] * v);
  vec1u_gradv = @(cutmat, muco) (muco * (cutmat' * [0; 1; 0] * u)) * (cutmat' * grad3(v));
  vec1u_vec1v = @(cutmat, muco) (muco * (cutmat' * [0; 1; 0] * u)) * (cutmat' * [0; 1; 0] * v);
  u_v = @(rhoco) ((rhoco*u)*v);

  %% Compute FE elementary matrices
  % ///////////////////////////////
  % Positive side
  tic; mat_gradu_gradv_pos = FEPack.pdes.Form.intg(mesh3Dpos.domain('volumic'), gradu_gradv(opts.cutmat_pos, mu3Dpos)); toc;
  tic; mat_gradu_vec1v_pos = FEPack.pdes.Form.intg(mesh3Dpos.domain('volumic'), gradu_vec1v(opts.cutmat_pos, mu3Dpos)); toc;
  tic; mat_vec1u_gradv_pos = FEPack.pdes.Form.intg(mesh3Dpos.domain('volumic'), vec1u_gradv(opts.cutmat_pos, mu3Dpos)); toc;
  tic; mat_vec1u_vec1v_pos = FEPack.pdes.Form.intg(mesh3Dpos.domain('volumic'), vec1u_vec1v(opts.cutmat_pos, mu3Dpos)); toc;
  tic; mat_u_v_pos         = FEPack.pdes.Form.intg(mesh3Dpos.domain('volumic'),        u_v(rho3Dpos)); toc;

  % Negative side
  tic; mat_gradu_gradv_neg = FEPack.pdes.Form.intg(mesh3Dneg.domain('volumic'), gradu_gradv(opts.cutmat_neg, mu3Dneg)); toc;
  tic; mat_gradu_vec1v_neg = FEPack.pdes.Form.intg(mesh3Dneg.domain('volumic'), gradu_vec1v(opts.cutmat_neg, mu3Dneg)); toc;
  tic; mat_vec1u_gradv_neg = FEPack.pdes.Form.intg(mesh3Dneg.domain('volumic'), vec1u_gradv(opts.cutmat_neg, mu3Dneg)); toc;
  tic; mat_vec1u_vec1v_neg = FEPack.pdes.Form.intg(mesh3Dneg.domain('volumic'), vec1u_vec1v(opts.cutmat_neg, mu3Dneg)); toc;
  tic; mat_u_v_neg         = FEPack.pdes.Form.intg(mesh3Dneg.domain('volumic'),        u_v(rho3Dneg)); toc;

  %% Boundary conditions
  % ////////////////////
  basis_functions = 'Lagrange';
  if strcmpi(basis_functions, 'Lagrange')
    BCstruct_pos.spB0 = FEPack.spaces.PeriodicLagrangeBasis(mesh3Dpos.domain('xmin'));
    BCstruct_pos.spB1 = FEPack.spaces.PeriodicLagrangeBasis(mesh3Dpos.domain('xmax'));
  else
    FourierIds = [0, numNodes3D/4, numNodes3D/4];
    BCstruct_pos.spB0 = spaces.FourierBasis(mesh3Dpos.domain('xmin'), FourierIds);
    BCstruct_pos.spB1 = spaces.FourierBasis(mesh3Dpos.domain('xmax'), FourierIds);
  end
  BCstruct_pos.BCdu = 0.0;
  BCstruct_pos.BCu = 1.0;
  BCstruct_pos.representation = 'projection';

  if strcmpi(basis_functions, 'Lagrange')
    BCstruct_neg.spB0 = FEPack.spaces.PeriodicLagrangeBasis(mesh3Dneg.domain('xmin'));
    BCstruct_neg.spB1 = FEPack.spaces.PeriodicLagrangeBasis(mesh3Dneg.domain('xmax'));
  else
    FourierIds = [0, numNodes3D/4, numNodes3D/4];
    BCstruct_neg.spB0 = spaces.FourierBasis(mesh3Dneg.domain('xmin'), FourierIds);
    BCstruct_neg.spB1 = spaces.FourierBasis(mesh3Dneg.domain('xmax'), FourierIds);
  end
  BCstruct_neg.BCdu = 0.0;
  BCstruct_neg.BCu = 1.0;
  BCstruct_neg.representation = 'projection';

  FloquetPoints_pos = linspace(-pi/period, pi/period, numFloquetPoints).';
  FloquetPoints_neg = linspace(-pi/period, pi/period, numFloquetPoints).';
  opts.computeSol = false;
  
  %% Save the workspace
  save([cheminDonnees, '/inputs_', int2str(numNodes), '.mat'], '-v7.3');
  
end