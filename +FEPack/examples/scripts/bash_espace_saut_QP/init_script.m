function init_script(numFloquetPoints)
  
  if (nargin < 1)
    numFloquetPoints = 50;
  end

  import FEPack.*

  %% Problem-related variables
  opts.omega = 5 + 0.25i;
  period = 1;
  vecper = [-1, 1]; % [-sqrt(2), 1];
  cutvec = [1/vecper(2); -vecper(1)/vecper(2)];
  opts.cutmat = [[1; 0; 0], [0; cutvec]];
  cutslope = cutvec(2) / cutvec(1);
  semiInfiniteDirection = 1;
  infiniteDirection = 2;
  numCellsSemiInfinite_pos = 7;
  numCellsSemiInfinite_neg = 7;
  numCellsInfinite = 6;
  opts.verbose = 0;
  FloquetPoints = linspace(-pi/period, pi/period, numFloquetPoints);

  %% Mesh
  pregenerate_mesh = 1;
  numNodes2D = 10;
  numNodes3D = 10;

  if pregenerate_mesh

    % Pick mesh from saved file
    m = load(['pregenMeshes/2D/unstruct_mesh_2D_', int2str(numNodes2D), '_positive.mat']);
    mesh2Dpos = m.mesh;

    m = load(['pregenMeshes/3D/unstruct_mesh_3D_', int2str(numNodes3D), '_positive.mat']);
    mesh3Dpos = m.mesh;

    m = load(['pregenMeshes/2D/unstruct_mesh_2D_', int2str(numNodes2D), '_negative_X.mat']);
    mesh2Dneg = m.mesh;

    m = load(['pregenMeshes/3D/unstruct_mesh_3D_', int2str(numNodes3D), '_negative_X.mat']);
    mesh3Dneg = m.mesh;

  else

    structmesh = 0; %#ok
    mesh2Dpos = meshes.MeshRectangle(structmesh, [0 1], [0 1], numNodes2D, numNodes2D);
    mesh3Dpos = meshes.MeshCuboid(structmesh, [0 1], [0 1], [0 1], numNodes3D, numNodes3D, numNodes3D);

    mesh2Dneg = meshes.MeshRectangle(structmesh, [0 -1], [0 1], numNodes2D, numNodes2D);
    mesh3Dneg = meshes.MeshCuboid(structmesh, [0 -1], [0 1], [0 1], numNodes3D, numNodes3D, numNodes3D);

  end

  %% Coefficients
  % Cartesian coefficients
  type_mu_pos = 2;  % 1: constant; 2: periodic; 
  type_mu_neg = 1;  % 1: constant; 2: periodic; 
  type_rho_pos = 2; % 1: constant; 2: periodic;  
  type_rho_neg = 1; % 1: constant; 2: periodic;  

  switch (type_mu_pos)
  case 1
    mu2DposCart = @(x) ones(size(x, 1), 1);
  case 2
    mu2DposCart = @(x) 1 + 0.25*cos(2*pi*x(:, 1)) + 0.25*cos(2*pi*x(:, 2));
  end

  switch (type_mu_neg)
  case 1
    mu2DnegCart = @(x) ones(size(x, 1), 1);
  case 2
    mu2DnegCart = @(x) 1 + 0.5*sin(2*pi*x(:, 1)).*sin(2*pi*x(:, 2));
  end

  switch (type_rho_pos)
  case 1
    rho2DposCart = @(x) ones(size(x, 1), 1);
  case 2
    rho2DposCart = @(x) 2 + 0.5*cos(2*pi*x(:, 1)).*sin(2*pi*x(:, 2));
  end

  switch (type_rho_neg)
  case 1
    rho2DnegCart = @(x) ones(size(x, 1), 1);
  case 2
    rho2DnegCart = @(x) 1 + 0.25*sin(2*pi*x(:, 1)) + 0.25*cos(2*pi*x(:, 2));
  end

  % 2D coefficients
  mu2Dpos  = @(x) mu2DposCart( [x(:, 1) + cutvec(2)*x(:, 2), cutvec(1)*x(:, 2), zeros(size(x, 1), 1)]);
  rho2Dpos = @(x) rho2DposCart([x(:, 1) + cutvec(2)*x(:, 2), cutvec(1)*x(:, 2), zeros(size(x, 1), 1)]);
  mu2Dneg  = @(x) mu2DnegCart( [x(:, 1) + cutvec(2)*x(:, 2), cutvec(1)*x(:, 2), zeros(size(x, 1), 1)]);
  rho2Dneg = @(x) rho2DnegCart([x(:, 1) + cutvec(2)*x(:, 2), cutvec(1)*x(:, 2), zeros(size(x, 1), 1)]);

  % 3D coefficients
  fun3D = @(fun2Dcart, x) fun2Dcart([x(:, 1) + x(:, 2), x(:, 3), zeros(size(x, 1), 1)]);
  mu3Dpos = @(x) fun3D(mu2DposCart, x);
  mu3Dneg = @(x) fun3D(mu2Dneg, x);
  rho3Dpos = @(x) fun3D(rho2DposCart, x);
  rho3Dneg = @(x) fun3D(rho2Dneg, x);

  % Jump data
  G = @(x) FEPack.tools.cutoff(x(:, 2), -0.3, 0.3);
  G3D = @(x) G([zeros(size(x, 1), 1), x(:, 1)/cutvec(1), zeros(size(x, 1), 1)]);

  %% Plot coefficients?
  plot_coefficients = false;
  if (plot_coefficients)
    set(groot,'defaultAxesTickLabelInterpreter','latex'); %#ok
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

    figure(3);
    x = [zeros(256, 1), linspace(-2, 2, 256)', zeros(256, 1)];
    plot(x(:, 2), G(x));
  end

  %% Bilinear and linear forms
  u = pdes.PDEObject; v = dual(u);
  gradu_gradv = @(muco) (muco * (opts.cutmat' * grad3(u))) * (opts.cutmat' * grad3(v));
  gradu_vec1v = @(muco) (muco * (opts.cutmat' * grad3(u))) * (opts.cutmat' * [0; 1; 0] * v);
  vec1u_gradv = @(muco) (muco * (opts.cutmat' * [0; 1; 0] * u)) * (opts.cutmat' * grad3(v));
  vec1u_vec1v = @(muco) (muco * (opts.cutmat' * [0; 1; 0] * u)) * (opts.cutmat' * [0; 1; 0] * v);
  u_v = @(rhoco) - (opts.omega^2) * ((rhoco*u)*v);

  %% Compute FE elementary matrices
  % Positive side
  mat_gradu_gradv_pos = FEPack.pdes.Form.intg(mesh3Dpos.domain('volumic'), gradu_gradv(mu3Dpos));
  mat_gradu_vec1v_pos = FEPack.pdes.Form.intg(mesh3Dpos.domain('volumic'), gradu_vec1v(mu3Dpos));
  mat_vec1u_gradv_pos = FEPack.pdes.Form.intg(mesh3Dpos.domain('volumic'), vec1u_gradv(mu3Dpos));
  mat_vec1u_vec1v_pos = FEPack.pdes.Form.intg(mesh3Dpos.domain('volumic'), vec1u_vec1v(mu3Dpos));
  mat_u_v_pos         = FEPack.pdes.Form.intg(mesh3Dpos.domain('volumic'),        u_v(rho3Dpos));

  % Negative side
  mat_gradu_gradv_neg = FEPack.pdes.Form.intg(mesh3Dneg.domain('volumic'), gradu_gradv(mu3Dneg));
  mat_gradu_vec1v_neg = FEPack.pdes.Form.intg(mesh3Dneg.domain('volumic'), gradu_vec1v(mu3Dneg));
  mat_vec1u_gradv_neg = FEPack.pdes.Form.intg(mesh3Dneg.domain('volumic'), vec1u_gradv(mu3Dneg));
  mat_vec1u_vec1v_neg = FEPack.pdes.Form.intg(mesh3Dneg.domain('volumic'), vec1u_vec1v(mu3Dneg));
  mat_u_v_neg         = FEPack.pdes.Form.intg(mesh3Dneg.domain('volumic'),        u_v(rho3Dneg));

  %% Boundary conditions
  basis_functions = 'Lagrange';
  if strcmpi(basis_functions, 'Lagrange')
    BCstruct_pos.spB0 = FEPack.spaces.PeriodicLagrangeBasis(mesh3Dpos.domain('xmin'));
    BCstruct_pos.spB1 = FEPack.spaces.PeriodicLagrangeBasis(mesh3Dpos.domain('xmax'));
  else
    FourierIds = [0, numNodes3D/2, 0];
    BCstruct_pos.spB0 = spaces.FourierBasis(mesh3Dpos.domain('xmin'), FourierIds);
    BCstruct_pos.spB1 = spaces.FourierBasis(mesh3Dpos.domain('xmax'), FourierIds);
  end
  BCstruct_pos.BCdu = 0.0;
  BCstruct_pos.BCu = 1.0;
  BCstruct_pos.representation = 'weak evaluation';

  if strcmpi(basis_functions, 'Lagrange')
    BCstruct_neg.spB0 = FEPack.spaces.PeriodicLagrangeBasis(mesh3Dneg.domain('xmin'));
    BCstruct_neg.spB1 = FEPack.spaces.PeriodicLagrangeBasis(mesh3Dneg.domain('xmax'));
  else
    FourierIds = [0, numNodes3D/2, 0];
    BCstruct_neg.spB0 = spaces.FourierBasis(mesh3Dneg.domain('xmin'), FourierIds);
    BCstruct_neg.spB1 = spaces.FourierBasis(mesh3Dneg.domain('xmax'), FourierIds);
  end
  BCstruct_neg.BCdu = 0.0;
  BCstruct_neg.BCu = 1.0;
  BCstruct_neg.representation = 'weak evaluation';

  %% Save the workspace
  save('inputs.mat');
  
end