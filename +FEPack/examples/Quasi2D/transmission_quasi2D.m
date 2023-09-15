clear; clc;
%%
import FEPack.*
% profile OFF
% profile ON

%% Problem-related variables
omega = 8 + 1i;
opts.period = 1;
opts.verbose = 0;
problem_setting = 'A'; % 'A' or 'B'

if strcmpi(problem_setting, 'A')

  % 2D coefficients
  period_posFun = 1;
  period_negFun = sqrt(2);

  mu2Dpos = @(x) ones(size(x, 1), 1);  % 0.5 + perCutoffCircle(x, [1; 0], [0; period_posFun], [0.5, 0.5], [-0.2, 0.2]);
  rho2Dpos = @(x) ones(size(x, 1), 1);  % 0.5 + perCutoffCuboid(x, [1; 0], [0; period_posFun], [0.5, 0.5], [-0.2, 0.2], [-0.2, 0.2], 0.5, 1);
  mu2Dneg  = @(x) ones(size(x, 1), 1);  % 1 + 0.5 * cos(2*pi*x(:, 1)) .* cos(2*pi*x(:, 2)/period_negFun);
  rho2Dneg = @(x) ones(size(x, 1), 1);  % 1 + 0.25 * sin(2*pi*x(:, 1)) + 0.25 * sin(2*pi*x(:, 2)/period_negFun);

  % Cut vector
  period_pos = period_posFun;
  period_neg = period_negFun;

  opts.cutvec = [period_pos; period_neg];

  % 3D coefficients
  mu3Dpos = @(x) mu2Dpos([x(:, 1), period_pos * x(:, 2)]);
  rho3Dpos = @(x) rho2Dpos([x(:, 1), period_pos * x(:, 2)]);
  mu3Dneg = @(x) mu2Dneg([x(:, 1), period_neg * x(:, 3)]);
  rho3Dneg = @(x) rho2Dneg([x(:, 1), period_neg * x(:, 3)]);
  
else

  % 2D coefficients
  vecperFun = [-sqrt(2), 1]; % [-sqrt(2), 1];
  
  mu2Dpos = @(x) 0.5 + perCutoffCircle(x, [1; 0], vecperFun, [0.5, 0.5], [-0.2, 0.2]);
  rho2Dpos = @(x) 0.5 + perCutoffCuboid(x, [1; 0], vecperFun, [0.5, 0.5], [-0.2, 0.2], [-0.2, 0.2], 0.5, 1);
  mu2Dneg  = @(x) ones(size(x, 1), 1);
  rho2Dneg = @(x) ones(size(x, 1), 1);

  % Cut vector
  vecper = vecperFun;
  opts.cutvec = [1/vecper(2); -vecper(1)/vecper(2)];
  
  % 3D coefficients
  fun3D = @(fun2D, x) fun2D([x(:, 1) + x(:, 3) + x(:, 2) * vecper(1), x(:, 2) * vecper(2)]);
  mu3Dpos = @(x) fun3D(mu2Dpos, x);
  mu3Dneg = @(x) fun3D(mu2Dneg, x);
  rho3Dpos = @(x) fun3D(rho2Dpos, x);
  rho3Dneg = @(x) fun3D(rho2Dneg, x);

end

% Cut matrix and cut slope
opts.cutmat = [[1; 0; 0], [0; opts.cutvec]];
cutslope = opts.cutvec(2) / opts.cutvec(1);

% Jump data
% alpha_G = 3;
% eps_G = 1e-8;
% supp_G = -log(eps_G) / alpha_G;
% G = @(x) exp(-alpha_G * x(:, 2).^2) .* (abs(x(:, 2)) <= supp_G);
G = @(x) FEPack.tools.cutoff(x(:, 2), -0.5, 0.5);
G3D = @(x) G([zeros(size(x, 1), 1), x(:, 2)/opts.cutvec(1), zeros(size(x, 1), 1)]);

% Numbers of cells
opts.numCellsSemiInfinite_pos = 6;
opts.numCellsSemiInfinite_neg = 6;
numCellsInfinite = 5;
opts.numFloquetPoints = 1;

%% Meshes
struct_mesh = 1;

numNodesX = 20;
numNodesY = 20;
numNodesZ = 20;

meshYZ = meshes.MeshRectangle(struct_mesh, [0 1], [0 1], numNodesY, numNodesZ);
meshZX = meshes.MeshRectangle(struct_mesh, [0 1], [0 1], numNodesZ, numNodesX);
meshXYpos = meshes.MeshRectangle(struct_mesh, [0  1], [0 1/opts.cutvec(1)], numNodesX, ceil(numNodesY/opts.cutvec(1)));
meshXYneg = meshes.MeshRectangle(struct_mesh, [0 -1], [0 1/opts.cutvec(1)], numNodesX, ceil(numNodesY/opts.cutvec(1)));

%% Basis functions and boundary conditions
basis_functions = 'Lagrange';
if strcmpi(basis_functions, 'Lagrange')
  BCstruct.spBX = FEPack.spaces.PeriodicLagrangeBasis(meshYZ.domain('volumic')); % X = cst
  BCstruct.spBY = FEPack.spaces.PeriodicLagrangeBasis(meshZX.domain('volumic')); % Y = cst
else
  BCstruct.spBX = spaces.FourierBasis(meshYZ.domain('xmin'), [floor(numNodesY/4), floor(numNodesZ/4)]);
  BCstruct.spBY = spaces.FourierBasis(meshZX.domain('xmax'), [floor(numNodesZ/4), floor(numNodesX/4)]);
end
BCstruct.BCdu = 0.0;
BCstruct.BCu = 1.0;
BCstruct.representation = 'projection';

%% Solve transmission problem
opts.suffix = 'pos';
quasi2DHalfGuideDirichlet(+1, meshYZ, meshZX, meshXYpos, omega, mu3Dpos, rho3Dpos, BCstruct, opts);